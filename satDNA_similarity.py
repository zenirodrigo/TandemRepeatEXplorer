#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Dependência opcional para acelerar:
  pip install parasail

Se parasail não estiver instalado, o script usa fallback em Python puro.
"""

import os
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Set

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

try:
    import parasail as _parasail
    _PARASAIL_AVAILABLE = True
except ImportError:
    _parasail = None
    _PARASAIL_AVAILABLE = False

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
MAX_CANDIDATES_PER_SEQ = 3000

if _PARASAIL_AVAILABLE:
    _DNA_MATRIX = _parasail.matrix_create("ACGTN", MATCH_SCORE, MISMATCH_SCORE)
else:
    _DNA_MATRIX = None


def clean_id(h: str) -> str:
    return h.strip().split()[0]


def clean_seq(s: str) -> str:
    s = s.strip().upper().replace("U", "T")
    return "".join(c for c in s if c in VALID)


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records = []
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((clean_id(header), clean_seq("".join(seq_lines))))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            records.append((clean_id(header), clean_seq("".join(seq_lines))))

    return [(seq_id, seq) for seq_id, seq in records if seq]


def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for seq_id, seq in records:
            out.write(f">{seq_id}\n")
            for i in range(0, len(seq), FASTA_WRAP):
                out.write(seq[i:i + FASTA_WRAP] + "\n")


def sanitize_filename(text: str) -> str:
    safe = []
    for char in text:
        if char.isalnum() or char in ("-", "_", "."):
            safe.append(char)
        else:
            safe.append("_")
    cleaned = "".join(safe).strip("._")
    return cleaned if cleaned else "unnamed"


def write_pairwise_alignment_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    """
    Escreve sequências alinhadas em FASTA.

    Diferente de write_fasta(), aqui preservamos '-' porque as sequências já
    são alinhamentos par-a-par.
    """
    with open(path, "w", encoding="utf-8") as out:
        for seq_id, aligned_seq in records:
            out.write(f">{seq_id}\n")
            for i in range(0, len(aligned_seq), FASTA_WRAP):
                out.write(aligned_seq[i:i + FASTA_WRAP] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def revcomp_kmer(kmer: str) -> str:
    return kmer.translate(DNA_COMP)[::-1]


def canonical_kmer(kmer: str) -> str:
    rc = revcomp_kmer(kmer)
    return kmer if kmer <= rc else rc


def choose_k_and_min_shared(seq_len: int) -> Tuple[int, int]:
    if seq_len >= 5000:
        return 13, 6
    if seq_len >= 2000:
        return 11, 5
    if seq_len >= 800:
        return 9, 4
    return 6, 3


def circular_kmers_set(seq: str, k: int) -> Set[str]:
    n = len(seq)
    if n < k:
        return set()

    seq2 = seq + seq
    kmers = set()
    for i in range(n):
        kmers.add(canonical_kmer(seq2[i:i + k]))
    return kmers


def _parasail_identity(A: str, B2: str) -> float:
    real_b_len = len(B2) // 2
    result = _parasail.sg_de_stats(
        A,
        B2,
        abs(GAP_SCORE),
        abs(GAP_SCORE),
        _DNA_MATRIX,
    )
    denom = max(len(A), real_b_len)
    return result.matches / denom if denom > 0 else 0.0


def _parasail_alignment(A: str, B2: str) -> Tuple[float, str, str]:
    real_b_len = len(B2) // 2
    result = _parasail.sg_de_trace(
        A,
        B2,
        abs(GAP_SCORE),
        abs(GAP_SCORE),
        _DNA_MATRIX,
    )
    traceback = result.traceback
    alnA = traceback.query
    alnB = traceback.ref
    matches = sum(1 for a, b in zip(alnA, alnB) if a == b and a != "-")
    denom = max(len(A), real_b_len)
    identity = matches / denom if denom > 0 else 0.0
    return identity, alnA, alnB


def _py_semiglobal_identity_only(A: str, B2: str) -> float:
    n = len(A)
    m = len(B2)
    real_b_len = m // 2

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
            diag = prev[j - 1] + (MATCH_SCORE if ai == B2[j - 1] else MISMATCH_SCORE)
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

    i = n
    j = best_j
    matches = 0

    while i > 0:
        move = tb[i][j]
        if move == 1:
            if A[i - 1] == B2[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif move == 2:
            i -= 1
        else:
            j -= 1

    denom = max(len(A), real_b_len)
    return matches / denom if denom > 0 else 0.0


def _py_semiglobal_alignment(A: str, B2: str) -> Tuple[float, str, str]:
    n = len(A)
    m = len(B2)
    real_b_len = m // 2

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
            diag = prev[j - 1] + (MATCH_SCORE if ai == B2[j - 1] else MISMATCH_SCORE)
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

    i = n
    j = best_j
    alnA = []
    alnB = []
    matches = 0

    while i > 0:
        move = tb[i][j]
        if move == 1:
            a = A[i - 1]
            b = B2[j - 1]
            alnA.append(a)
            alnB.append(b)
            if a == b:
                matches += 1
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

    denom = max(len(A), real_b_len)
    identity = matches / denom if denom > 0 else 0.0
    return identity, alnA, alnB


def semiglobal_identity_only(A: str, B2: str) -> float:
    if _PARASAIL_AVAILABLE:
        return _parasail_identity(A, B2)
    return _py_semiglobal_identity_only(A, B2)


def semiglobal_alignment(A: str, B2: str) -> Tuple[float, str, str]:
    if _PARASAIL_AVAILABLE:
        return _parasail_alignment(A, B2)
    return _py_semiglobal_alignment(A, B2)


def best_direction_identity(A: str, B: str) -> Tuple[float, str]:
    id_fwd = semiglobal_identity_only(A, B + B)
    rcB = revcomp(B)
    id_rc = semiglobal_identity_only(A, rcB + rcB)
    if id_fwd >= id_rc:
        return id_fwd, "forward"
    return id_rc, "reverse-complement"


def best_direction_alignment(A: str, B: str) -> Tuple[float, str, str, str]:
    id_fwd, alnA_fwd, alnB_fwd = semiglobal_alignment(A, B + B)
    rcB = revcomp(B)
    id_rc, alnA_rc, alnB_rc = semiglobal_alignment(A, rcB + rcB)
    if id_fwd >= id_rc:
        return id_fwd, "forward", alnA_fwd, alnB_fwd
    return id_rc, "reverse-complement", alnA_rc, alnB_rc


def _py_semiglobal_alignment_with_ref_coords(A: str, B2: str) -> Tuple[float, str, str, int, int]:

    n = len(A)
    m = len(B2)
    real_b_len = m // 2

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
    ref_end = best_j

    i, j = n, best_j
    alnA = []
    alnB = []
    matches = 0

    while i > 0:
        move = tb[i][j]
        if move == 1:
            a = A[i - 1]
            b = B2[j - 1]
            alnA.append(a)
            alnB.append(b)
            if a == b:
                matches += 1
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

    ref_start = j
    alnA = "".join(reversed(alnA))
    alnB = "".join(reversed(alnB))

    denom = max(len(A), real_b_len)
    identity = matches / denom if denom > 0 else 0.0
    return identity, alnA, alnB, ref_start, ref_end


def best_direction_alignment_with_member_rotation(reference: str, member: str) -> Tuple[float, str, str, str, int, str]:
    
    if not member:
        return 0.0, "forward", reference, "", 0, member

    id_fwd, aln_ref_fwd, aln_mem_fwd, start_fwd, _end_fwd = _py_semiglobal_alignment_with_ref_coords(
        reference,
        member + member,
    )

    rc_member = revcomp(member)
    id_rc, aln_ref_rc, aln_mem_rc, start_rc, _end_rc = _py_semiglobal_alignment_with_ref_coords(
        reference,
        rc_member + rc_member,
    )

    if id_fwd >= id_rc:
        start = start_fwd % len(member)
        rotated = member[start:] + member[:start]
        return id_fwd, "forward", aln_ref_fwd, aln_mem_fwd, start, rotated

    start = start_rc % len(rc_member)
    rotated = rc_member[start:] + rc_member[:start]
    return id_rc, "reverse-complement", aln_ref_rc, aln_mem_rc, start, rotated



def rotate(seq: str, start: int) -> str:
    if not seq:
        return seq
    start %= len(seq)
    return seq[start:] + seq[:start]


def score_pair(a: str, b: str) -> int:
    return MATCH_SCORE if a == b else MISMATCH_SCORE


def nw_global_align(A: str, B: str) -> Tuple[int, str, str]:
    n, m = len(A), len(B)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]  # 1 diag, 2 up, 3 left

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_SCORE
        tb[i][0] = 2
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + GAP_SCORE
        tb[0][j] = 3

    for i in range(1, n + 1):
        ai = A[i - 1]
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + score_pair(ai, B[j - 1])
            up = dp[i - 1][j] + GAP_SCORE
            left = dp[i][j - 1] + GAP_SCORE
            best = diag
            move = 1
            if up > best:
                best = up
                move = 2
            if left > best:
                best = left
                move = 3
            dp[i][j] = best
            tb[i][j] = move

    i, j = n, m
    alnA = []
    alnB = []
    while i > 0 or j > 0:
        move = tb[i][j]
        if i > 0 and j > 0 and move == 1:
            alnA.append(A[i - 1])
            alnB.append(B[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or move == 2):
            alnA.append(A[i - 1])
            alnB.append("-")
            i -= 1
        else:
            alnA.append("-")
            alnB.append(B[j - 1])
            j -= 1

    return dp[n][m], "".join(reversed(alnA)), "".join(reversed(alnB))


def semi_global_ref_to_dimer(reference: str, member_dimer: str) -> Tuple[int, str, str, int, int]:

    n, m = len(reference), len(member_dimer)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_SCORE
        tb[i][0] = 2
    for j in range(1, m + 1):
        dp[0][j] = 0
        tb[0][j] = 3

    for i in range(1, n + 1):
        ai = reference[i - 1]
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + score_pair(ai, member_dimer[j - 1])
            up = dp[i - 1][j] + GAP_SCORE
            left = dp[i][j - 1] + GAP_SCORE
            best = diag
            move = 1
            if up > best:
                best = up
                move = 2
            if left > best:
                best = left
                move = 3
            dp[i][j] = best
            tb[i][j] = move

    best_j = max(range(m + 1), key=lambda j: dp[n][j])
    best_score = dp[n][best_j]

    i, j = n, best_j
    aln_ref = []
    aln_mem = []
    consumed_positions = []

    while i > 0:
        move = tb[i][j]
        if move == 1:
            aln_ref.append(reference[i - 1])
            aln_mem.append(member_dimer[j - 1])
            consumed_positions.append(j - 1)
            i -= 1
            j -= 1
        elif move == 2:
            aln_ref.append(reference[i - 1])
            aln_mem.append("-")
            i -= 1
        else:
            aln_ref.append("-")
            aln_mem.append(member_dimer[j - 1])
            consumed_positions.append(j - 1)
            j -= 1

    aln_ref = "".join(reversed(aln_ref))
    aln_mem = "".join(reversed(aln_mem))

    if consumed_positions:
        start = min(consumed_positions)
        end = max(consumed_positions) + 1
    else:
        start = best_j
        end = best_j

    return best_score, aln_ref, aln_mem, start, end


def identity_from_alignment(alnA: str, alnB: str, denom: int) -> float:
    matches = sum(1 for a, b in zip(alnA, alnB) if a != "-" and b != "-" and a == b)
    return matches / denom if denom else 0.0


def best_frame_against_first(rep_seq: str, member_seq: str) -> Dict[str, object]:

    if not member_seq:
        _, aln_rep, aln_mem = nw_global_align(rep_seq, member_seq)
        return {
            "orientation": "forward",
            "rotation_start": 0,
            "rotated_seq": member_seq,
            "semiglobal_score": 0,
            "global_score": 0,
            "frame_identity": 0.0,
            "global_identity": 0.0,
            "aln_rep": aln_rep,
            "aln_member": aln_mem,
        }

    best = None
    for orientation, oriented in (("forward", member_seq), ("reverse-complement", revcomp(member_seq))):
        dimer = oriented + oriented
        sg_score, sg_rep, sg_mem, start, _end = semi_global_ref_to_dimer(rep_seq, dimer)
        rotation_start = start % len(oriented)
        rotated = rotate(oriented, rotation_start)
        global_score, aln_rep, aln_member = nw_global_align(rep_seq, rotated)
        frame_identity = identity_from_alignment(sg_rep, sg_mem, max(len(rep_seq), len(member_seq)))
        global_identity = identity_from_alignment(aln_rep, aln_member, max(len(rep_seq), len(member_seq)))
        candidate = {
            "orientation": orientation,
            "rotation_start": rotation_start,
            "rotated_seq": rotated,
            "semiglobal_score": sg_score,
            "global_score": global_score,
            "frame_identity": frame_identity,
            "global_identity": global_identity,
            "aln_rep": aln_rep,
            "aln_member": aln_member,
        }
        key = (global_score, global_identity, sg_score, frame_identity)
        if best is None or key > best[0]:
            best = (key, candidate)
    return best[1]


def pairwise_to_reference_slots(aln_ref: str, aln_member: str, rep_len: int):

    insertions = {i: [] for i in range(-1, rep_len)}
    bases = ["-"] * rep_len
    ref_pos = -1

    for a, b in zip(aln_ref, aln_member):
        if a != "-":
            ref_pos += 1
            if ref_pos < rep_len:
                bases[ref_pos] = b
        else:
            insertions[ref_pos].append(b)

    return bases, insertions


def build_reference_anchored_msa(rep_id: str, rep_seq: str, aligned_members: List[Dict[str, object]]) -> List[Tuple[str, str]]:

    rep_len = len(rep_seq)
    member_slots = []
    max_ins = {i: 0 for i in range(-1, rep_len)}

    for d in aligned_members:
        bases, insertions = pairwise_to_reference_slots(d["aln_rep"], d["aln_member"], rep_len)
        for k, v in insertions.items():
            max_ins[k] = max(max_ins[k], len(v))
        member_slots.append((d, bases, insertions))

    rep_aligned = []
    rep_aligned.extend("-" * max_ins[-1])
    for p, base in enumerate(rep_seq):
        rep_aligned.append(base)
        rep_aligned.extend("-" * max_ins[p])

    records = [(f"{rep_id}|role=representative_first_sequence|original_len={rep_len}", "".join(rep_aligned))]

    for d, bases, insertions in member_slots:
        out = []
        ins = insertions[-1]
        out.extend(ins)
        out.extend("-" * (max_ins[-1] - len(ins)))
        for p in range(rep_len):
            out.append(bases[p])
            ins = insertions[p]
            out.extend(ins)
            out.extend("-" * (max_ins[p] - len(ins)))
        records.append((d["header"], "".join(out)))

    return records

def choose_family_medoid(
    members: List[int],
    ids: List[str],
    seqs: List[str],
    cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]],
    alignments_counter: List[int],
) -> int:

    if len(members) == 1:
        return members[0]

    best_member = members[0]
    best_key = (-1.0, -1, "")
    for candidate in members:
        total = 0.0
        comparisons = 0
        for other in members:
            if other == candidate:
                continue
            id_min, _, _, _, _ = get_reciprocal_identity_cached(
                candidate,
                other,
                seqs,
                cache,
                alignments_counter,
            )
            total += id_min
            comparisons += 1
        mean_id = total / comparisons if comparisons else 1.0
        key = (mean_id, len(seqs[candidate]), "".join(chr(255 - ord(c)) for c in ids[candidate]))
        if key > best_key:
            best_key = key
            best_member = candidate
    return best_member


def reciprocal_identity_only(A: str, B: str) -> Tuple[float, float, str, float, str]:
    id1, rel1 = best_direction_identity(A, B)
    id2, rel2 = best_direction_identity(B, A)
    id_min = min(id1, id2)
    return id_min, id1, rel1, id2, rel2


def pretty_alignment(alnA: str, alnB: str) -> str:
    mid = []
    for a, b in zip(alnA, alnB):
        mid.append("|" if a != "-" and b != "-" and a == b else " ")
    mid = "".join(mid)

    out = []
    for i in range(0, len(alnA), ALIGN_WRAP):
        out.append("A: " + alnA[i:i + ALIGN_WRAP])
        out.append("   " + mid[i:i + ALIGN_WRAP])
        out.append("B: " + alnB[i:i + ALIGN_WRAP])
        out.append("")
    return "\n".join(out).rstrip()


def ask_user() -> Tuple[str, float]:
    print("\n=== satDNA Similarity Clustering — complete-linkage corrigido ===")
    print("Agrupa monômeros de satDNA por identidade circular recíproca.")
    print("Complete-linkage real: fusão de famílias somente se TODOS passam contra TODOS.\n")

    if _PARASAIL_AVAILABLE:
        print("parasail detectado: alinhamentos rápidos via SIMD/C.\n")
    else:
        print("AVISO: parasail não detectado. O fallback em Python puro será mais lento.")
        print("Para acelerar: pip install parasail\n")

    fasta = input("Arquivo FASTA de entrada: ").strip()
    fasta = os.path.expanduser(fasta)

    identity = float(input(
        "Threshold (ex: 0.80 = 80%) [0.80]: "
    ).strip() or "0.80")

    return fasta, identity


def pair_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def get_reciprocal_identity_cached(
    a: int,
    b: int,
    seqs: List[str],
    cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]],
    alignments_counter: List[int],
) -> Tuple[float, float, str, float, str]:
    key = pair_key(a, b)
    if key not in cache:
        cache[key] = reciprocal_identity_only(seqs[a], seqs[b])
        alignments_counter[0] += 1
    return cache[key]


def candidate_passes_all_members(
    candidate: int,
    members: List[int],
    seqs: List[str],
    threshold: float,
    cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]],
    alignments_counter: List[int],
) -> bool:
    for member in members:
        if member == candidate:
            continue
        id_min, _, _, _, _ = get_reciprocal_identity_cached(
            candidate,
            member,
            seqs,
            cache,
            alignments_counter,
        )
        if id_min + 1e-12 < threshold:
            return False
    return True


def families_pass_all_vs_all(
    members_a: List[int],
    members_b: List[int],
    seqs: List[str],
    threshold: float,
    cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]],
    alignments_counter: List[int],
) -> bool:

    for a in members_a:
        for b in members_b:
            if a == b:
                continue
            id_min, _, _, _, _ = get_reciprocal_identity_cached(
                a,
                b,
                seqs,
                cache,
                alignments_counter,
            )
            if id_min + 1e-12 < threshold:
                return False
    return True


def make_proof_edge(
    a: int,
    b: int,
    ids: List[str],
    seqs: List[str],
    lens: List[int],
    id_min: float,
) -> dict:
    id_ab, rel_ab, alnA_ab, alnB_ab = best_direction_alignment(seqs[a], seqs[b])
    id_ba, rel_ba, alnA_ba, alnB_ba = best_direction_alignment(seqs[b], seqs[a])
    return {
        "a": a,
        "b": b,
        "id_min": id_min,
        "A_id": ids[a],
        "B_id": ids[b],
        "A_len": lens[a],
        "B_len": lens[b],
        "A_to_B": {
            "identity": id_ab,
            "relation": rel_ab,
            "alnA": alnA_ab,
            "alnB": alnB_ab,
        },
        "B_to_A": {
            "identity": id_ba,
            "relation": rel_ba,
            "alnA": alnA_ba,
            "alnB": alnB_ba,
        },
    }



def write_family_alignment_outputs(
    out_dir: str,
    sorted_families: List[Tuple[int, List[int]]],
    ids: List[str],
    seqs: List[str],
    lens: List[int],
    cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]],
    alignments_counter: List[int],
) -> Tuple[int, int, int, int]:

    os.makedirs(out_dir, exist_ok=True)

    family_files_written = 0
    frame_corrected_files_written = 0
    pairwise_files_written = 0
    anchored_files_written = 0

    manifest_path = os.path.join(out_dir, "MANIFEST.tsv")
    with open(manifest_path, "w", encoding="utf-8") as manifest:
        manifest.write(
            "family_rank\tfamily_id\tfamily_size\trep_id\t"
            "members_fasta\tframe_corrected_monomers_fasta\t"
            "pairwise_to_first_representative_fasta\t"
            "firstseq_reference_anchored_alignment_fasta\talignment_summary_tsv\n"
        )

        for family_rank, (_fid, members) in enumerate(sorted_families, 1):
            # Ignora singletons nesta etapa final de alinhamento.
            if len(members) < 2:
                continue

            # Ordem original do FASTA de entrada. A primeira sequência da família manda.
            ordered_members = sorted(members)
            representative = ordered_members[0]
            rep_id = ids[representative]
            rep_seq = seqs[representative]
            fam_id = rep_id
            safe_fam_id = sanitize_filename(fam_id)
            prefix = f"Family_{family_rank:06d}_{safe_fam_id}"

            members_fasta = os.path.join(out_dir, f"{prefix}.members.fasta")
            frame_corrected_fasta = os.path.join(out_dir, f"{prefix}.frame_corrected_monomers.fasta")
            pairwise_fasta = os.path.join(out_dir, f"{prefix}.pairwise_to_first_representative.fasta")
            anchored_fasta = os.path.join(out_dir, f"{prefix}.firstseq_reference_anchored_alignment.fasta")
            summary_tsv = os.path.join(out_dir, f"{prefix}.alignment_summary.tsv")

            # 1. FASTA bruto da família.
            member_records = [(ids[m], seqs[m]) for m in ordered_members]
            write_fasta(members_fasta, member_records)
            family_files_written += 1

            # 2-4. Frame corrigido, pairwise com gaps, e alinhamento ancorado.
            aligned_members: List[Dict[str, object]] = []
            frame_corrected_records: List[Tuple[str, str]] = [(
                f"{rep_id}|family_rank={family_rank}|role=representative_first_sequence|"
                f"orientation=forward|rotation_start=0|original_len={lens[representative]}",
                rep_seq,
            )]
            pairwise_records: List[Tuple[str, str]] = []

            with open(summary_tsv, "w", encoding="utf-8") as summary:
                summary.write(
                    "family_rank\tfamily_id\trepresentative\tmember\torientation\t"
                    "rotation_start\tframe_identity\tglobal_identity_to_rep\t"
                    "original_len\trotated_len\tnote\n"
                )
                summary.write(
                    f"{family_rank}\t{fam_id}\t{rep_id}\t{rep_id}\tforward\t0\t"
                    f"1.000000\t1.000000\t{lens[representative]}\t{lens[representative]}\t"
                    f"representative_first_sequence\n"
                )

                for member in ordered_members[1:]:
                    d = best_frame_against_first(rep_seq, seqs[member])
                    d["member_id"] = ids[member]
                    d["header"] = (
                        f"{ids[member]}|family_rank={family_rank}|family={fam_id}|"
                        f"representative={rep_id}|role=aligned_to_first_representative|"
                        f"orientation={d['orientation']}|rotation_start={d['rotation_start']}|"
                        f"frame_identity={d['frame_identity']:.6f}|"
                        f"global_identity_to_rep={d['global_identity']:.6f}|"
                        f"original_len={lens[member]}"
                    )
                    aligned_members.append(d)

                    frame_corrected_records.append((
                        f"{ids[member]}|family_rank={family_rank}|family={fam_id}|"
                        f"representative={rep_id}|role=frame_corrected_monomer|"
                        f"orientation={d['orientation']}|rotation_start={d['rotation_start']}|"
                        f"frame_identity={d['frame_identity']:.6f}|"
                        f"global_identity_to_rep={d['global_identity']:.6f}|"
                        f"original_len={lens[member]}",
                        str(d["rotated_seq"]),
                    ))

                    pair_tag = f"{sanitize_filename(rep_id)}__VS__{sanitize_filename(ids[member])}"
                    pairwise_records.append((
                        f"{pair_tag}|role=representative_first_sequence|member={ids[member]}|"
                        f"orientation={d['orientation']}|rotation_start={d['rotation_start']}|"
                        f"global_identity_to_rep={d['global_identity']:.6f}",
                        str(d["aln_rep"]),
                    ))
                    pairwise_records.append((
                        f"{pair_tag}|role=member_aligned_to_first_representative|member={ids[member]}|"
                        f"orientation={d['orientation']}|rotation_start={d['rotation_start']}|"
                        f"global_identity_to_rep={d['global_identity']:.6f}",
                        str(d["aln_member"]),
                    ))

                    summary.write(
                        f"{family_rank}\t{fam_id}\t{rep_id}\t{ids[member]}\t{d['orientation']}\t"
                        f"{d['rotation_start']}\t{d['frame_identity']:.6f}\t"
                        f"{d['global_identity']:.6f}\t{lens[member]}\t{len(str(d['rotated_seq']))}\t"
                        f"member_dimerized_oriented_rotated_and_aligned_to_first_representative\n"
                    )

            anchored_records = build_reference_anchored_msa(rep_id, rep_seq, aligned_members)

            write_fasta(frame_corrected_fasta, frame_corrected_records)
            write_pairwise_alignment_fasta(pairwise_fasta, pairwise_records)
            write_pairwise_alignment_fasta(anchored_fasta, anchored_records)

            frame_corrected_files_written += 1
            pairwise_files_written += 1
            anchored_files_written += 1

            manifest.write(
                f"{family_rank}\t{fam_id}\t{len(members)}\t{rep_id}\t"
                f"{members_fasta}\t{frame_corrected_fasta}\t{pairwise_fasta}\t"
                f"{anchored_fasta}\t{summary_tsv}\n"
            )

    return (
        family_files_written,
        frame_corrected_files_written,
        pairwise_files_written,
        anchored_files_written,
    )

def run(fasta: str, identity_threshold: float) -> None:
    t0 = time.time()

    records = read_fasta(fasta)
    if not records:
        raise SystemExit(f"ERRO: nenhuma sequência encontrada em {fasta}")

    ids = [record[0] for record in records]
    seqs = [record[1] for record in records]
    lens = [len(seq) for seq in seqs]
    n = len(seqs)

    min_len = min(lens)
    max_len = max(lens)
    med_len = sorted(lens)[n // 2]

    print(f"\nSequências carregadas: {n}")
    print(f"Comprimento dos monômeros: min={min_len} bp | mediana={med_len} bp | max={max_len} bp")
    print("Construindo índice de k-mers circulares canônicos...")

    K, MIN_SHARED = choose_k_and_min_shared(med_len)
    print(f"Pré-filtro: k={K}, min_shared_signals={MIN_SHARED}, max_candidates_per_seq={MAX_CANDIDATES_PER_SEQ}")
    print("Clusterizando com complete-linkage real...\n")

    inv: Dict[str, List[int]] = defaultdict(list)
    kmers_by_i: List[Set[str]] = [set() for _ in range(n)]

    for i, seq in enumerate(seqs):
        kmset = circular_kmers_set(seq, K)
        kmers_by_i[i] = kmset
        for kmer in kmset:
            inv[kmer].append(i)

    seq_to_family: List[int] = [-1] * n
    family_members: Dict[int, List[int]] = {}

    identity_cache: Dict[Tuple[int, int], Tuple[float, float, str, float, str]] = {}

    compared_candidates = 0
    alignments_counter = [0]
    unions = 0

    proof_edges: List[dict] = []
    superfamily_hits: List[dict] = []

    for i in range(n):
        if i % PROGRESS_EVERY_SEQS == 0 and i > 0:
            elapsed = time.time() - t0
            singleton_count = sum(1 for x in seq_to_family if x == -1)
            fam_count = len(family_members) + singleton_count
            print(
                f"Progresso: {i}/{n} | alinhamentos={alignments_counter[0]} | "
                f"pares_cache={len(identity_cache)} | famílias_aprox={fam_count} | elapsed={elapsed:.1f}s"
            )

        counts: Dict[int, int] = defaultdict(int)
        for kmer in kmers_by_i[i]:
            for j in inv.get(kmer, []):
                if j >= i:
                    continue
                counts[j] += 1

        candidates = [(j, c) for j, c in counts.items() if c >= MIN_SHARED]
        candidates.sort(key=lambda item: item[1], reverse=True)
        if len(candidates) > MAX_CANDIDATES_PER_SEQ:
            candidates = candidates[:MAX_CANDIDATES_PER_SEQ]

        superfams_for_i = 0

        for j, shared in candidates:
            compared_candidates += 1

            id_min_ij, id_i_to_j, rel_i_to_j, id_j_to_i, rel_j_to_i = get_reciprocal_identity_cached(
                i,
                j,
                seqs,
                identity_cache,
                alignments_counter,
            )

            if alignments_counter[0] % PROGRESS_EVERY_ALNS == 0:
                elapsed = time.time() - t0
                print(f"  Alinhamentos computados: {alignments_counter[0]} | elapsed={elapsed:.1f}s")

            if id_min_ij + 1e-12 >= identity_threshold:
                fid_i = seq_to_family[i]
                fid_j = seq_to_family[j]

                if fid_i == -1 and fid_j == -1:
                    fid_new = min(i, j)
                    family_members[fid_new] = sorted([i, j])
                    seq_to_family[i] = fid_new
                    seq_to_family[j] = fid_new
                    unions += 1
                    proof_edges.append(make_proof_edge(i, j, ids, seqs, lens, id_min_ij))
                    continue

                if fid_i == -1 and fid_j != -1:
                    members_j = family_members[fid_j]
                    ok = candidate_passes_all_members(
                        i,
                        members_j,
                        seqs,
                        identity_threshold,
                        identity_cache,
                        alignments_counter,
                    )
                    if ok:
                        family_members[fid_j].append(i)
                        family_members[fid_j] = sorted(set(family_members[fid_j]))
                        seq_to_family[i] = fid_j
                        unions += 1
                        proof_edges.append(make_proof_edge(i, j, ids, seqs, lens, id_min_ij))
                    continue

                if fid_i != -1 and fid_j == -1:
                    members_i = family_members[fid_i]
                    ok = candidate_passes_all_members(
                        j,
                        members_i,
                        seqs,
                        identity_threshold,
                        identity_cache,
                        alignments_counter,
                    )
                    if ok:
                        family_members[fid_i].append(j)
                        family_members[fid_i] = sorted(set(family_members[fid_i]))
                        seq_to_family[j] = fid_i
                        unions += 1
                        proof_edges.append(make_proof_edge(i, j, ids, seqs, lens, id_min_ij))
                    continue

                if fid_i != -1 and fid_j != -1:
                    if fid_i == fid_j:
                        continue

                    members_i = family_members[fid_i]
                    members_j = family_members[fid_j]

                    can_merge = families_pass_all_vs_all(
                        members_i,
                        members_j,
                        seqs,
                        identity_threshold,
                        identity_cache,
                        alignments_counter,
                    )

                    if can_merge:
                        fid_keep = min(fid_i, fid_j)
                        fid_drop = max(fid_i, fid_j)

                        merged = sorted(set(family_members[fid_keep] + family_members[fid_drop]))
                        family_members[fid_keep] = merged
                        del family_members[fid_drop]

                        for member in merged:
                            seq_to_family[member] = fid_keep

                        unions += 1
                        proof_edges.append(make_proof_edge(i, j, ids, seqs, lens, id_min_ij))
                    continue

            one_way_best = max(id_i_to_j, id_j_to_i)
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
                    "query_i": query,
                    "target_i": target,
                    "query_id": ids[query],
                    "target_id": ids[target],
                    "query_len": lens[query],
                    "target_len": lens[target],
                    "direction": direction,
                    "best_orientation": rel,
                    "identity_oneway": id_oneway,
                    "identity_reciprocal_min": id_min_ij,
                    "shared_signals": shared,
                })
                superfams_for_i += 1

    final_families: Dict[int, List[int]] = dict(family_members)
    for i in range(n):
        if seq_to_family[i] == -1:
            final_families[i] = [i]
            seq_to_family[i] = i

    sorted_families = sorted(final_families.items(), key=lambda item: len(item[1]), reverse=True)

    base = os.path.splitext(fasta)[0]
    pct = int(identity_threshold * 100)

    out_fasta = f"{base}.id{pct}.family_reps.fasta"
    out_tsv = f"{base}.id{pct}.families.tsv"
    out_report = f"{base}.id{pct}.proof.txt"
    out_super = f"{base}.id{pct}.superfamilies.tsv"
    out_super_groups = f"{base}.id{pct}.superfamilies.groups.txt"
    out_family_fastas_dir = f"{base}.id{pct}.family_fastas"

    reps: List[Tuple[str, str]] = []
    member_to_rep: Dict[str, str] = {}
    rep_to_size: Dict[str, int] = {}
    index_to_final_fid: Dict[int, int] = {}

    for fid, members in sorted_families:
        for member in members:
            index_to_final_fid[member] = fid

    with open(out_tsv, "w", encoding="utf-8") as tsv:
        tsv.write("family_id\tfamily_size\trep_id\tmember_id\tmember_len\n")
        for fid, members in sorted_families:
            rep = min(members, key=lambda x: ids[x])
            fam_id = ids[rep]
            reps.append((fam_id, seqs[rep]))
            rep_to_size[fam_id] = len(members)
            for member in sorted(members, key=lambda x: ids[x]):
                member_to_rep[ids[member]] = fam_id
                tsv.write(f"{fam_id}\t{len(members)}\t{ids[rep]}\t{ids[member]}\t{lens[member]}\n")

    write_fasta(out_fasta, reps)

    (
        family_files_written,
        frame_corrected_files_written,
        pairwise_files_written,
        anchored_files_written,
    ) = write_family_alignment_outputs(
        out_family_fastas_dir,
        sorted_families,
        ids,
        seqs,
        lens,
        identity_cache,
        alignments_counter,
    )

    with open(out_super, "w", encoding="utf-8") as out:
        out.write(
            "query_id\ttarget_id\tquery_len\ttarget_len\t"
            "query_family\ttarget_family\tquery_family_size\ttarget_family_size\t"
            "direction\tbest_orientation\tidentity_oneway\tidentity_reciprocal_min\tshared_signals\n"
        )
        superfamily_hits.sort(key=lambda d: (d["identity_oneway"], d["shared_signals"]), reverse=True)
        for hit in superfamily_hits:
            qfam = member_to_rep.get(hit["query_id"], hit["query_id"])
            tfam = member_to_rep.get(hit["target_id"], hit["target_id"])
            out.write(
                f"{hit['query_id']}\t{hit['target_id']}\t{hit['query_len']}\t{hit['target_len']}\t"
                f"{qfam}\t{tfam}\t{rep_to_size.get(qfam, 1)}\t{rep_to_size.get(tfam, 1)}\t"
                f"{hit['direction']}\t{hit['best_orientation']}\t{hit['identity_oneway']:.4f}\t"
                f"{hit['identity_reciprocal_min']:.4f}\t{hit['shared_signals']}\n"
            )

    fam_edges: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    fam_adj: Dict[str, Set[str]] = defaultdict(set)

    for hit in superfamily_hits:
        qfam = member_to_rep.get(hit["query_id"], hit["query_id"])
        tfam = member_to_rep.get(hit["target_id"], hit["target_id"])
        if qfam == tfam:
            continue
        a, b = (qfam, tfam) if qfam <= tfam else (tfam, qfam)
        fam_edges[(a, b)].append(hit)
        fam_adj[a].add(b)
        fam_adj[b].add(a)

    visited = set()
    components = []
    for node in sorted(fam_adj.keys()):
        if node in visited:
            continue
        stack = [node]
        visited.add(node)
        comp = []
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in fam_adj.get(x, set()):
                if y not in visited:
                    visited.add(y)
                    stack.append(y)
        components.append(sorted(comp))

    with open(out_super_groups, "w", encoding="utf-8") as out:
        out.write("# Grupos de superfamílias\n")
        out.write("# Famílias aqui NÃO foram fundidas porque falharam no critério complete-linkage recíproco.\n")
        out.write("# O link indica similaridade one-way >= threshold entre pelo menos um par de monômeros.\n\n")
        out.write(f"# Threshold: {identity_threshold}\n")
        out.write(f"# Links de superfamília no nível de pares: {len(superfamily_hits)}\n")
        out.write(f"# Grupos no nível de famílias: {len(components)}\n\n")

        if not components:
            out.write("Nenhum grupo de superfamília detectado.\n")
        else:
            for group_id, comp in enumerate(sorted(components, key=len, reverse=True), 1):
                out.write("=" * 100 + "\n")
                out.write(f"SuperfamilyGroup-{group_id} (families={len(comp)}):\n")
                for fam in comp:
                    out.write(f"  Family rep: {fam} (size={rep_to_size.get(fam, 1)})\n")
                out.write("\nEvidências, top links por identidade one-way:\n")

                comp_set = set(comp)
                evid_list = []
                for (a, b), evidences in fam_edges.items():
                    if a in comp_set and b in comp_set:
                        best = max(evidences, key=lambda e: (e["identity_oneway"], e["shared_signals"]))
                        evid_list.append((a, b, best))

                evid_list.sort(key=lambda x: (x[2]["identity_oneway"], x[2]["shared_signals"]), reverse=True)

                for printed, (a, b, best) in enumerate(evid_list):
                    if printed >= 200:
                        out.write("  (parado após 200 links)\n")
                        break
                    out.write(
                        f"  {a}  <->  {b} | best one-way={best['identity_oneway'] * 100:.1f}% "
                        f"(dir={best['direction']}, orient={best['best_orientation']}) | "
                        f"ex: {best['query_id']} -> {best['target_id']} | "
                        f"len {best['query_len']} -> {best['target_len']}\n"
                    )
                out.write("\n")

    proofs_by_final_fid: Dict[int, List[dict]] = defaultdict(list)
    for proof in proof_edges:
        final_fid = index_to_final_fid.get(proof["a"], seq_to_family[proof["a"]])
        proofs_by_final_fid[final_fid].append(proof)

    report = []
    report.append("# similarity2_firstseq_anchored.py — relatório de prova")
    report.append("# Método: complete-linkage real.")
    report.append("# Um novo membro entra na família somente se passar o threshold contra TODOS os membros já presentes.")
    report.append("# Duas famílias só fundem se TODOS os membros de uma passarem contra TODOS os membros da outra.")
    report.append(f"# Threshold de identidade: {identity_threshold}")
    report.append(f"# Pré-filtro: k={K}, min_shared_signals={MIN_SHARED}, max_candidates={MAX_CANDIDATES_PER_SEQ}")
    report.append(f"# Scores: match={MATCH_SCORE}, mismatch={MISMATCH_SCORE}, gap={GAP_SCORE}")
    report.append(f"# Backend: {'parasail' if _PARASAIL_AVAILABLE else 'Python puro'}")
    report.append("")
    report.append(
        f"# Stats: sequences={n}, families={len(sorted_families)}, "
        f"candidate_pairs={compared_candidates}, alignments={alignments_counter[0]}, "
        f"cached_pairs={len(identity_cache)}, unions={unions}"
    )
    report.append(f"# Superfamily links: {out_super}")
    report.append(f"# Superfamily groups: {out_super_groups}")
    report.append(f"# Family FASTA/alignment directory: {out_family_fastas_dir}")
    report.append("")

    for fid, members in sorted_families:
        rep = min(members, key=lambda x: ids[x])
        fam_id = ids[rep]
        report.append("=" * 100)
        report.append(f"FAMILY: {fam_id}")
        report.append(f"Family size: {len(members)}")
        report.append("Members (id\tlength):")
        for member in sorted(members, key=lambda x: ids[x]):
            report.append(f"  {ids[member]}\t{lens[member]}")
        report.append("")

        proofs = proofs_by_final_fid.get(fid, [])
        if not proofs:
            report.append("Singleton or family without stored proof edges.")
            report.append("")
            continue

        report.append("Proof alignments, union edges accepted during clustering:")
        report.append("")

        for printed, proof in enumerate(proofs):
            report.append(f"[EDGE] {proof['A_id']} (len={proof['A_len']})  <->  {proof['B_id']} (len={proof['B_len']})")
            report.append(f"  Reciprocal identity min: {proof['id_min']:.4f}")
            report.append(f"  A→B: {proof['A_to_B']['identity']:.4f} ({proof['A_to_B']['relation']})")
            report.append(f"  B→A: {proof['B_to_A']['identity']:.4f} ({proof['B_to_A']['relation']})")
            report.append("")
            report.append("  Alignment A→B:")
            report.append(pretty_alignment(proof["A_to_B"]["alnA"], proof["A_to_B"]["alnB"]))
            report.append("")
            report.append("  Alignment B→A:")
            report.append(pretty_alignment(proof["B_to_A"]["alnA"], proof["B_to_A"]["alnB"]))
            report.append("")

            if printed + 1 >= MAX_PROOFS_PER_FAMILY:
                report.append(f"(Stopped after {MAX_PROOFS_PER_FAMILY} proofs for this family.)")
                report.append("")
                break

    with open(out_report, "w", encoding="utf-8") as out:
        out.write("\n".join(report) + "\n")

    elapsed = time.time() - t0
    print("\nCONCLUÍDO")
    print(f"Tempo total:             {elapsed:.1f}s")
    print(f"Sequências:              {n}")
    print(f"Famílias:                {len(sorted_families)}")
    print(f"Pares candidatos:        {compared_candidates}")
    print(f"Alinhamentos computados: {alignments_counter[0]}")
    print(f"Pares no cache:          {len(identity_cache)}")
    print(f"Uniões aceitas:          {unions}")
    print(f"FASTA reps:              {out_fasta}")
    print(f"TSV:                     {out_tsv}")
    print(f"RELATÓRIO:               {out_report}")
    print(f"SUPERFAMS:               {out_super}")
    print(f"SF_GROUPS:               {out_super_groups}")
    print(f"FAMILY_FASTAS_DIR:       {out_family_fastas_dir}")
    print(f"Family member FASTAs:    {family_files_written}")
    print(f"Frame-corrected FASTAs:  {frame_corrected_files_written}")
    print(f"Pairwise-to-first FASTAs:{pairwise_files_written}")
    print(f"Anchored alignment FASTAs:{anchored_files_written}")


if __name__ == "__main__":
    fasta_path, identity_value = ask_user()
    run(fasta_path, identity_value)

