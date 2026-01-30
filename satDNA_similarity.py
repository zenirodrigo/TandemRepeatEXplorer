#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
from collections import defaultdict
from typing import List, Tuple, Dict

try:
    import readline
    import glob

    def _complete_path(text, state):
        buf = readline.get_line_buffer()
        token = buf.split()[-1] if buf.split() else ""
        token = os.path.expanduser(token)
        matches = glob.glob(token + "*")
        matches = [m + ("/" if os.path.isdir(m) else "") for m in matches]
        try:
            return matches[state]
        except IndexError:
            return None

    readline.set_completer_delims(" \t\n;")
    readline.parse_and_bind("tab: complete")
    readline.set_completer(_complete_path)
except Exception:
    pass

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
VALID = set("ACGTN")

FASTA_WRAP = 60

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

PROGRESS_EVERY_SEQS = 25
PROGRESS_EVERY_ALNS = 200

MAX_PROOFS_PER_FAMILY = 10
ALIGN_WRAP = 120

MAX_SUPERFAMS_PER_QUERY = 200


def clean_id(h: str) -> str:
    return h.strip().split()[0]


def clean_seq(s: str) -> str:
    s = s.strip().upper().replace("U", "T")
    return "".join(c for c in s if c in VALID)


def read_fasta(path: str) -> List[Tuple[str, str]]:
    rec = []
    h = None
    buf = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    rec.append((clean_id(h), clean_seq("".join(buf))))
                h = line[1:]
                buf = []
            else:
                buf.append(line)
        if h is not None:
            rec.append((clean_id(h), clean_seq("".join(buf))))
    rec = [(i, s) for i, s in rec if s]
    return rec


def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for hid, seq in records:
            out.write(f">{hid}\n")
            for i in range(0, len(seq), FASTA_WRAP):
                out.write(seq[i:i + FASTA_WRAP] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def revcomp_kmer(km: str) -> str:
    return km.translate(DNA_COMP)[::-1]


def canonical_kmer(km: str) -> str:
    rc = revcomp_kmer(km)
    return km if km <= rc else rc


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        self.parent[rb] = ra
        return True


def choose_k_and_min_shared(seq_len: int) -> Tuple[int, int]:
    if seq_len >= 5000:
        return 13, 6
    if seq_len >= 2000:
        return 11, 5
    if seq_len >= 800:
        return 9, 4
    return 6, 3


def circular_kmers_set(seq: str, k: int) -> set:
    """
    RC-aware circular k-mers:
    we index canonical(kmer) = min(kmer, RC(kmer))
    so reverse-complement matches are discoverable in the prefilter.
    """
    n = len(seq)
    if n < k:
        return set()
    s2 = seq + seq
    out = set()
    for i in range(0, n):
        km = s2[i:i + k]
        out.add(canonical_kmer(km))
    return out


def semiglobal_identity_only(A: str, B2: str) -> float:
    n = len(A)
    m = len(B2)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_SCORE
        tb[i][0] = 2

    for j in range(1, m + 1):
        dp[0][j] = 0
        tb[0][j] = 3

    for i in range(1, n + 1):
        ai = A[i - 1]
        row = dp[i]
        row_tb = tb[i]
        prev = dp[i - 1]
        for j in range(1, m + 1):
            bj = B2[j - 1]
            diag = prev[j - 1] + (MATCH_SCORE if ai == bj else MISMATCH_SCORE)
            up = prev[j] + GAP_SCORE
            left = row[j - 1] + GAP_SCORE

            best = diag
            move = 1
            if up > best:
                best = up
                move = 2
            if left > best:
                best = left
                move = 3

            row[j] = best
            row_tb[j] = move

    best_j = max(range(m + 1), key=lambda j: dp[n][j])

    i, j = n, best_j
    matches = 0
    aln_len = 0
    while i > 0:
        move = tb[i][j]
        if move == 1:
            if A[i - 1] == B2[j - 1]:
                matches += 1
            aln_len += 1
            i -= 1
            j -= 1
        elif move == 2:
            aln_len += 1
            i -= 1
        else:
            aln_len += 1
            j -= 1

    return matches / aln_len if aln_len > 0 else 0.0


def semiglobal_alignment(A: str, B2: str) -> Tuple[float, str, str]:
    n = len(A)
    m = len(B2)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_SCORE
        tb[i][0] = 2
    for j in range(1, m + 1):
        dp[0][j] = 0
        tb[0][j] = 3

    for i in range(1, n + 1):
        ai = A[i - 1]
        row = dp[i]
        row_tb = tb[i]
        prev = dp[i - 1]
        for j in range(1, m + 1):
            bj = B2[j - 1]
            diag = prev[j - 1] + (MATCH_SCORE if ai == bj else MISMATCH_SCORE)
            up = prev[j] + GAP_SCORE
            left = row[j - 1] + GAP_SCORE

            best = diag
            move = 1
            if up > best:
                best = up
                move = 2
            if left > best:
                best = left
                move = 3

            row[j] = best
            row_tb[j] = move

    best_j = max(range(m + 1), key=lambda j: dp[n][j])

    i, j = n, best_j
    alnA = []
    alnB = []
    while i > 0:
        move = tb[i][j]
        if move == 1:
            alnA.append(A[i - 1])
            alnB.append(B2[j - 1])
            i -= 1
            j -= 1
        elif move == 2:
            alnA.append(A[i - 1])
            alnB.append("-")
            i -= 1
        else:
            alnA.append("-")
            alnB.append(B2[j - 1])
            j -= 1

    alnA = "".join(reversed(alnA))
    alnB = "".join(reversed(alnB))

    matches = 0
    L = len(alnA)
    for a, b in zip(alnA, alnB):
        if a != "-" and b != "-" and a == b:
            matches += 1
    identity = matches / L if L > 0 else 0.0

    return identity, alnA, alnB


def best_direction_identity(A: str, B: str) -> Tuple[float, str]:
    id_fwd = semiglobal_identity_only(A, B + B)
    rcB = revcomp(B)
    id_rc = semiglobal_identity_only(A, rcB + rcB)
    if id_fwd >= id_rc:
        return id_fwd, "forward"
    return id_rc, "reverse-complement"


def best_direction_alignment(A: str, B: str) -> Tuple[float, str, str, str]:
    id_fwd, a1, b1 = semiglobal_alignment(A, B + B)
    rcB = revcomp(B)
    id_rc, a2, b2 = semiglobal_alignment(A, rcB + rcB)
    if id_fwd >= id_rc:
        return id_fwd, "forward", a1, b1
    return id_rc, "reverse-complement", a2, b2


def reciprocal_identity_only(A: str, B: str) -> Tuple[float, float, str, float, str]:
    id1, rel1 = best_direction_identity(A, B)
    id2, rel2 = best_direction_identity(B, A)
    id_min = id1 if id1 <= id2 else id2
    return id_min, id1, rel1, id2, rel2


def pretty_alignment(alnA: str, alnB: str) -> str:
    mid = []
    for a, b in zip(alnA, alnB):
        if a != "-" and b != "-" and a == b:
            mid.append("|")
        else:
            mid.append(" ")
    mid = "".join(mid)

    out = []
    for i in range(0, len(alnA), ALIGN_WRAP):
        out.append("A: " + alnA[i:i + ALIGN_WRAP])
        out.append("   " + mid[i:i + ALIGN_WRAP])
        out.append("B: " + alnB[i:i + ALIGN_WRAP])
        out.append("")
    return "\n".join(out).rstrip()


def ask_user() -> Tuple[str, float]:
    print("\n=== satDNA Similarity Clustering ===")
    print("This program groups satellite DNA monomers into families based on")
    print("circular sequence similarity, allowing reverse-complement and gaps.\n")
    print("Tip: use TAB to autocomplete file paths.\n")

    fasta = input("Input FASTA file with monomer sequences: ").strip()
    fasta = os.path.expanduser(fasta)

    identity = float(input(
        "Minimum identity threshold (e.g. 0.80 = 80% similarity) [0.80]: "
    ).strip() or "0.80")

    return fasta, identity


def run(fasta: str, identity_threshold: float) -> None:
    t0 = time.time()
    records = read_fasta(fasta)
    if not records:
        raise SystemExit(f"ERROR: No sequences found in {fasta}")

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    lens = [len(s) for s in seqs]
    n = len(seqs)

    min_len = min(lens)
    max_len = max(lens)
    med_len = sorted(lens)[n // 2]

    print(f"\nLoaded {n} sequences.")
    print(f"Length stats (bp): min={min_len} median={med_len} max={max_len}")
    print("Building k-mer index...")

    K, MIN_SHARED = choose_k_and_min_shared(med_len)
    MAX_CANDIDATES_PER_SEQ = 3000

    print(f"Internal prefilter: k={K}, min_shared_signals={MIN_SHARED}, max_candidates_per_seq={MAX_CANDIDATES_PER_SEQ}")
    print("Clustering + detecting superfamilies...")

    uf = UnionFind(n)

    inv: Dict[str, List[int]] = defaultdict(list)
    kmers_by_i: List[set] = [set() for _ in range(n)]

    for i, s in enumerate(seqs):
        kmset = circular_kmers_set(s, K)
        kmers_by_i[i] = kmset
        for km in kmset:
            inv[km].append(i)

    compared_candidates = 0
    alignments_done = 0
    unions = 0

    proof_edges = []
    superfamily_hits = []

    for i in range(n):
        if i % PROGRESS_EVERY_SEQS == 0 and i > 0:
            dt = time.time() - t0
            print(f"Progress: {i}/{n} | alignments={alignments_done} | unions={unions} | elapsed={dt:.1f}s")

        counts = defaultdict(int)
        for km in kmers_by_i[i]:
            for j in inv.get(km, []):
                if j <= i:
                    continue
                counts[j] += 1

        candidates = [(j, c) for j, c in counts.items() if c >= MIN_SHARED]
        candidates.sort(key=lambda x: x[1], reverse=True)
        if len(candidates) > MAX_CANDIDATES_PER_SEQ:
            candidates = candidates[:MAX_CANDIDATES_PER_SEQ]
        superfams_for_i = 0

        for j, shared in candidates:
            compared_candidates += 1

            id_min, id_i_to_j, rel_i_to_j, id_j_to_i, rel_j_to_i = reciprocal_identity_only(seqs[i], seqs[j])
            alignments_done += 1

            if alignments_done % PROGRESS_EVERY_ALNS == 0:
                dt = time.time() - t0
                print(f"  Alignments computed: {alignments_done} | elapsed={dt:.1f}s")

            if id_min + 1e-12 >= identity_threshold:
                if uf.find(i) == uf.find(j):
                    continue

                if uf.union(i, j):
                    unions += 1

                    id_ab, rel_ab, alnA_ab, alnB_ab = best_direction_alignment(seqs[i], seqs[j])
                    id_ba, rel_ba, alnA_ba, alnB_ba = best_direction_alignment(seqs[j], seqs[i])

                    proof_edges.append((
                        i, j,
                        {
                            "id_min": id_min,
                            "A_id": ids[i], "B_id": ids[j],
                            "A_len": lens[i], "B_len": lens[j],
                            "A_to_B": {"identity": id_ab, "relation": rel_ab, "alnA": alnA_ab, "alnB": alnB_ab},
                            "B_to_A": {"identity": id_ba, "relation": rel_ba, "alnA": alnA_ba, "alnB": alnB_ba},
                        }
                    ))
                continue

            one_way_best = id_i_to_j if id_i_to_j >= id_j_to_i else id_j_to_i
            if one_way_best + 1e-12 >= identity_threshold:
                if superfams_for_i >= MAX_SUPERFAMS_PER_QUERY:
                    continue

                if id_i_to_j >= id_j_to_i:
                    query = i
                    target = j
                    rel = rel_i_to_j
                    id_oneway = id_i_to_j
                    direction = "A_in_B"
                else:
                    query = j
                    target = i
                    rel = rel_j_to_i
                    id_oneway = id_j_to_i
                    direction = "B_in_A"

                superfamily_hits.append({
                    "query_id": ids[query],
                    "target_id": ids[target],
                    "query_len": lens[query],
                    "target_len": lens[target],
                    "direction": direction,
                    "best_orientation": rel,
                    "identity_oneway": id_oneway,
                    "identity_reciprocal_min": id_min,
                    "shared_signals": shared,
                })
                superfams_for_i += 1

    families = defaultdict(list)
    for i in range(n):
        families[uf.find(i)].append(i)

    base = os.path.splitext(fasta)[0]
    out_fasta = base + f".id{int(identity_threshold*100)}.family_reps.fasta"
    out_tsv = base + f".id{int(identity_threshold*100)}.families.tsv"
    out_report = base + f".id{int(identity_threshold*100)}.proof.txt"
    out_super = base + f".id{int(identity_threshold*100)}.superfamilies.tsv"

    reps = []
    with open(out_tsv, "w", encoding="utf-8") as tsv:
        tsv.write("family_id\tfamily_size\trep_id\tmember_id\tmember_len\n")
        for root, members in sorted(families.items(), key=lambda x: len(x[1]), reverse=True):
            rep = min(members, key=lambda x: ids[x])
            fam_id = ids[rep]
            reps.append((fam_id, seqs[rep]))
            for m in sorted(members, key=lambda x: ids[x]):
                tsv.write(f"{fam_id}\t{len(members)}\t{ids[rep]}\t{ids[m]}\t{lens[m]}\n")

    write_fasta(out_fasta, reps)

    with open(out_super, "w", encoding="utf-8") as out:
        out.write(
            "query_id\ttarget_id\tquery_len\ttarget_len\tdirection\tbest_orientation\t"
            "identity_oneway\tidentity_reciprocal_min\tshared_signals\n"
        )
        superfamily_hits.sort(key=lambda d: (d["identity_oneway"], d["shared_signals"]), reverse=True)
        for d in superfamily_hits:
            out.write(
                f"{d['query_id']}\t{d['target_id']}\t{d['query_len']}\t{d['target_len']}\t"
                f"{d['direction']}\t{d['best_orientation']}\t{d['identity_oneway']:.4f}\t"
                f"{d['identity_reciprocal_min']:.4f}\t{d['shared_signals']}\n"
            )

    edges_by_root = defaultdict(list)
    for i, j, proof in proof_edges:
        r = uf.find(i)
        if uf.find(j) == r:
            edges_by_root[r].append(proof)

    report = []
    report.append("# satDNA_similarity.py proof report")
    report.append("# Families: reciprocal circular identity (forward/RC) with gaps included.")
    report.append(f"# Identity threshold (family): {identity_threshold}")
    report.append(f"# Prefilter: k={K}, min_shared_signals={MIN_SHARED}, max_candidates_per_seq={MAX_CANDIDATES_PER_SEQ}")
    report.append(f"# Scoring: match={MATCH_SCORE}, mismatch={MISMATCH_SCORE}, gap={GAP_SCORE}")
    report.append("")
    report.append(f"# Stats: sequences={n}, families={len(families)}, candidate_pairs={compared_candidates}, alignments={alignments_done}, unions={unions}")
    report.append(f"# Superfamily links written to: {out_super}")
    report.append("")

    for root, members in sorted(families.items(), key=lambda x: len(x[1]), reverse=True):
        rep = min(members, key=lambda x: ids[x])
        fam_id = ids[rep]
        report.append("=" * 100)
        report.append(f"FAMILY: {fam_id}")
        report.append(f"Family size: {len(members)}")
        report.append("Members (id\tlength):")
        for m in sorted(members, key=lambda x: ids[x]):
            report.append(f"  {ids[m]}\t{lens[m]}")
        report.append("")

        proofs = edges_by_root.get(root, [])
        if not proofs:
            report.append("No stored union-edge proofs for this family (singleton or no union edges stored).")
            report.append("")
            continue

        report.append("Proof alignments (union edges used to build the family):")
        report.append("")

        printed = 0
        for p in proofs:
            report.append(f"[EDGE] {p['A_id']} (len={p['A_len']})  <->  {p['B_id']} (len={p['B_len']})")
            report.append(f"  Reciprocal identity (min of both directions): {p['id_min']:.4f}")
            report.append(f"  A→B best identity = {p['A_to_B']['identity']:.4f} (orientation: {p['A_to_B']['relation']})")
            report.append(f"  B→A best identity = {p['B_to_A']['identity']:.4f} (orientation: {p['B_to_A']['relation']})")
            report.append("")
            report.append("  Alignment A→B:")
            report.append(pretty_alignment(p["A_to_B"]["alnA"], p["A_to_B"]["alnB"]))
            report.append("")
            report.append("  Alignment B→A:")
            report.append(pretty_alignment(p["B_to_A"]["alnA"], p["B_to_A"]["alnB"]))
            report.append("")
            printed += 1
            if printed >= MAX_PROOFS_PER_FAMILY:
                report.append(f"(Stopped after {MAX_PROOFS_PER_FAMILY} proofs for this family.)")
                report.append("")
                break

    with open(out_report, "w", encoding="utf-8") as out:
        out.write("\n".join(report) + "\n")

    dt = time.time() - t0
    print("\nDONE")
    print(f"Total time: {dt:.1f}s")
    print(f"Sequences: {n}")
    print(f"Families:  {len(families)}")
    print(f"Candidate pairs evaluated: {compared_candidates}")
    print(f"Alignments computed:       {alignments_done}")
    print(f"Unions (family links):     {unions}")
    print(f"FASTA reps:  {out_fasta}")
    print(f"TSV:         {out_tsv}")
    print(f"REPORT:      {out_report}")
    print(f"SUPERFAMS:   {out_super}")


if __name__ == "__main__":
    fasta, identity = ask_user()
    run(fasta, identity)
