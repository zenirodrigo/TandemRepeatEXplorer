#!/usr/bin/env python3
# -*- coding: utf-8 -*-

Uso:
----
python3 satDNA_similarity.py \
    --fasta sequences.fasta \
    --out-prefix name_output \
    --id 0.80 \
    --cov 0.80 \
    --match 2 \
    --mismatch -1 \
    --gap -2 \
    --min-len 1


import argparse
import math
import os
import sys
from collections import defaultdict
from itertools import combinations


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def clean_seq(seq: str) -> str:
    seq = seq.strip().replace(" ", "").replace("\r", "").replace("\n", "")
    seq = seq.upper()
    allowed = set("ACGTN")
    seq = "".join([b if b in allowed else "N" for b in seq])
    return seq


def read_fasta(path: str):
    records = []
    header = None
    seq_chunks = []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = clean_seq("".join(seq_chunks))
                    records.append((header, seq))
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    if header is not None:
        seq = clean_seq("".join(seq_chunks))
        records.append((header, seq))

    return records


def write_fasta(records, out_path, width=80):
    with open(out_path, "w", encoding="utf-8") as out:
        for rid, seq in records:
            out.write(f">{rid}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False

        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1
        return True



def smith_waterman_local(query, target, match_score=2, mismatch_score=-1, gap_score=-2):
    """
    Retorna o melhor alinhamento local entre query e target.

    Saída:
    {
        'score': int,
        'aligned_query': str,
        'aligned_target': str,
        'q_start': int,  # 0-based inclusive no original query
        'q_end': int,    # 0-based exclusive no original query
        't_start': int,  # 0-based inclusive no original target
        't_end': int     # 0-based exclusive no original target
    }
    """
    n = len(query)
    m = len(target)

    H = [[0] * (m + 1) for _ in range(n + 1)]
    TB = [[0] * (m + 1) for _ in range(n + 1)]  # 0 stop, 1 diag, 2 up, 3 left

    best_score = 0
    best_pos = (0, 0)

    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            tj = target[j - 1]
            diag = H[i - 1][j - 1] + (match_score if qi == tj else mismatch_score)
            up = H[i - 1][j] + gap_score
            left = H[i][j - 1] + gap_score
            best = max(0, diag, up, left)

            H[i][j] = best
            if best == 0:
                TB[i][j] = 0
            elif best == diag:
                TB[i][j] = 1
            elif best == up:
                TB[i][j] = 2
            else:
                TB[i][j] = 3

            if best > best_score:
                best_score = best
                best_pos = (i, j)

    i, j = best_pos
    q_end = i
    t_end = j

    aln_q = []
    aln_t = []

    while i > 0 and j > 0 and H[i][j] > 0:
        move = TB[i][j]
        if move == 1:
            aln_q.append(query[i - 1])
            aln_t.append(target[j - 1])
            i -= 1
            j -= 1
        elif move == 2:
            aln_q.append(query[i - 1])
            aln_t.append("-")
            i -= 1
        elif move == 3:
            aln_q.append("-")
            aln_t.append(target[j - 1])
            j -= 1
        else:
            break

    q_start = i
    t_start = j

    aln_q.reverse()
    aln_t.reverse()

    return {
        "score": best_score,
        "aligned_query": "".join(aln_q),
        "aligned_target": "".join(aln_t),
        "q_start": q_start,
        "q_end": q_end,
        "t_start": t_start,
        "t_end": t_end,
    }


def compute_alignment_stats(aln_q: str, aln_t: str, len_a: int, len_b: int):
    """
    Calcula métricas do alinhamento usando normalização pelo maior comprimento.
    """
    matches = 0
    mismatches = 0
    gap_cols = 0
    aligned_a = 0
    aligned_b = 0

    for ca, cb in zip(aln_q, aln_t):
        if ca != "-":
            aligned_a += 1
        if cb != "-":
            aligned_b += 1

        if ca == "-" or cb == "-":
            gap_cols += 1
        elif ca == cb:
            matches += 1
        else:
            mismatches += 1

    max_len = max(len_a, len_b)
    min_aligned = min(aligned_a, aligned_b)

    normalized_identity = matches / max_len if max_len > 0 else 0.0
    coverage = min_aligned / max_len if max_len > 0 else 0.0

    aligned_columns = len(aln_q)

    return {
        "matches": matches,
        "mismatches": mismatches,
        "gap_cols": gap_cols,
        "aligned_a": aligned_a,
        "aligned_b": aligned_b,
        "aligned_columns": aligned_columns,
        "normalized_identity": normalized_identity,
        "coverage": coverage,
    }


def pretty_midline(aln_q: str, aln_t: str) -> str:
    chars = []
    for a, b in zip(aln_q, aln_t):
        if a == b and a != "-":
            chars.append("|")
        else:
            chars.append(" ")
    return "".join(chars)



def circular_local_best(query, target, match_score=2, mismatch_score=-1, gap_score=-2):
    """
    Busca o melhor alinhamento local de query contra target circular.
    Testa:
      - target forward
      - target reverse-complement
    usando target duplicado para permitir wrap-around.

    A melhor solução é escolhida por:
      1) maior score
      2) maior normalized_identity
      3) maior coverage

    Restringimos o trecho real usado no target a no máximo len(target),
    para evitar alinhar mais de uma volta inteira.
    """
    candidates = [
        ("forward", target),
        ("reverse-complement", reverse_complement(target)),
    ]

    best = None
    tlen = len(target)

    for orient, tgt in candidates:
        doubled = tgt + tgt
        result = smith_waterman_local(
            query=query,
            target=doubled,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
        )

        aln_q = result["aligned_query"]
        aln_t = result["aligned_target"]

        stats = compute_alignment_stats(aln_q, aln_t, len(query), len(target))

        # checagem do span real sobre o target duplicado
        t_start = result["t_start"]
        t_end = result["t_end"]
        raw_span = t_end - t_start

        # Se o span passa de uma volta, descartamos.
        # Isso evita "colar" mais de uma cópia do target duplicado.
        if raw_span > tlen:
            continue

        item = {
            "orientation": orient,
            "score": result["score"],
            "aligned_query": aln_q,
            "aligned_target": aln_t,
            "q_start": result["q_start"],
            "q_end": result["q_end"],
            "t_start_doubled": t_start,
            "t_end_doubled": t_end,
            "t_start_circular": t_start % tlen if tlen > 0 else 0,
            "t_end_circular": t_end % tlen if tlen > 0 else 0,
            "stats": stats,
        }

        if best is None:
            best = item
        else:
            bstats = best["stats"]
            if (
                item["score"] > best["score"]
                or (item["score"] == best["score"] and item["stats"]["normalized_identity"] > bstats["normalized_identity"])
                or (
                    item["score"] == best["score"]
                    and math.isclose(item["stats"]["normalized_identity"], bstats["normalized_identity"])
                    and item["stats"]["coverage"] > bstats["coverage"]
                )
            ):
                best = item

    if best is None:
        # fallback impossível na prática se houver bases válidas,
        # mas deixamos por segurança
        best = {
            "orientation": "forward",
            "score": 0,
            "aligned_query": "",
            "aligned_target": "",
            "q_start": 0,
            "q_end": 0,
            "t_start_doubled": 0,
            "t_end_doubled": 0,
            "t_start_circular": 0,
            "t_end_circular": 0,
            "stats": {
                "matches": 0,
                "mismatches": 0,
                "gap_cols": 0,
                "aligned_a": 0,
                "aligned_b": 0,
                "aligned_columns": 0,
                "normalized_identity": 0.0,
                "coverage": 0.0,
            },
        }

    return best


def reciprocal_pair_metrics(seq_a, seq_b, match_score=2, mismatch_score=-1, gap_score=-2):
    """
    Calcula:
      A -> B circular
      B -> A circular

    E define:
      reciprocal_identity = min(identity_A_to_B, identity_B_to_A)
      reciprocal_coverage = min(coverage_A_to_B, coverage_B_to_A)
    """
    a_to_b = circular_local_best(
        query=seq_a,
        target=seq_b,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score,
    )

    b_to_a = circular_local_best(
        query=seq_b,
        target=seq_a,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score,
    )

    rid = min(a_to_b["stats"]["normalized_identity"], b_to_a["stats"]["normalized_identity"])
    rcov = min(a_to_b["stats"]["coverage"], b_to_a["stats"]["coverage"])

    return {
        "a_to_b": a_to_b,
        "b_to_a": b_to_a,
        "reciprocal_identity": rid,
        "reciprocal_coverage": rcov,
    }



def choose_family_representative(member_ids, seq_dict):
    """
    Escolhe o representante da família.
    Regra:
      1) maior comprimento
      2) em empate, ordem alfabética do ID
    """
    return sorted(member_ids, key=lambda x: (-len(seq_dict[x]), x))[0]



def write_families_tsv(families, seq_dict, out_tsv):
    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("Family\tMemberID\tLength\n")
        for fam_name, members in families:
            for mid in members:
                out.write(f"{fam_name}\t{mid}\t{len(seq_dict[mid])}\n")


def write_pairwise_tsv(pair_rows, out_tsv):
    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write(
            "SeqA\tLenA\tSeqB\tLenB\t"
            "A_to_B_identity\tA_to_B_coverage\tA_to_B_orientation\t"
            "B_to_A_identity\tB_to_A_coverage\tB_to_A_orientation\t"
            "ReciprocalIdentity\tReciprocalCoverage\tPass\n"
        )
        for r in pair_rows:
            out.write(
                f"{r['id_a']}\t{r['len_a']}\t{r['id_b']}\t{r['len_b']}\t"
                f"{r['a_to_b_identity']:.4f}\t{r['a_to_b_coverage']:.4f}\t{r['a_to_b_orientation']}\t"
                f"{r['b_to_a_identity']:.4f}\t{r['b_to_a_coverage']:.4f}\t{r['b_to_a_orientation']}\t"
                f"{r['reciprocal_identity']:.4f}\t{r['reciprocal_coverage']:.4f}\t"
                f"{'YES' if r['passed'] else 'NO'}\n"
            )


def write_proof_report(
    out_path,
    families,
    seq_dict,
    used_edges,
    pair_metrics_map,
    args,
    total_pairs_evaluated,
    total_unions,
):
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("# satDNA_similarity.py proof report\n")
        out.write("# Families: circular local alignment with normalized identity + coverage filter.\n")
        out.write(f"# Identity threshold (family): {args.id}\n")
        out.write(f"# Coverage threshold (family): {args.cov}\n")
        out.write(f"# Scoring: match={args.match}, mismatch={args.mismatch}, gap={args.gap}\n")
        out.write("\n")
        out.write(
            f"# Stats: sequences={len(seq_dict)}, families={len(families)}, "
            f"pairwise_comparisons={total_pairs_evaluated}, unions={total_unions}\n"
        )
        out.write("\n")

        for fam_name, members in families:
            out.write("=" * 100 + "\n")
            out.write(f"FAMILY: {fam_name}\n")
            out.write(f"Family size: {len(members)}\n")
            out.write("Members (id\tlength):\n")
            for mid in members:
                out.write(f"  {mid}\t{len(seq_dict[mid])}\n")
            out.write("\n")

            fam_edges = []
            mset = set(members)
            for a, b in used_edges:
                if a in mset and b in mset:
                    fam_edges.append((a, b))

            if not fam_edges:
                out.write("No stored union-edge proofs for this family (singleton or no union edges stored).\n\n")
                continue

            out.write("Proof alignments (union edges used to build the family):\n\n")

            for a, b in fam_edges:
                key = tuple(sorted((a, b)))
                pm = pair_metrics_map[key]

                out.write(f"[EDGE] {a} (len={len(seq_dict[a])})  <->  {b} (len={len(seq_dict[b])})\n")
                out.write(f"  Reciprocal normalized identity: {pm['reciprocal_identity']:.4f}\n")
                out.write(f"  Reciprocal coverage:            {pm['reciprocal_coverage']:.4f}\n")
                out.write(
                    f"  A→B identity = {pm['a_to_b']['stats']['normalized_identity']:.4f} "
                    f"(coverage = {pm['a_to_b']['stats']['coverage']:.4f}, "
                    f"orientation: {pm['a_to_b']['orientation']})\n"
                )
                out.write(
                    f"  B→A identity = {pm['b_to_a']['stats']['normalized_identity']:.4f} "
                    f"(coverage = {pm['b_to_a']['stats']['coverage']:.4f}, "
                    f"orientation: {pm['b_to_a']['orientation']})\n"
                )
                out.write("\n")

                out.write("  Alignment A→B:\n")
                aq = pm["a_to_b"]["aligned_query"]
                at = pm["a_to_b"]["aligned_target"]
                mid = pretty_midline(aq, at)
                out.write(f"A: {aq}\n")
                out.write(f"   {mid}\n")
                out.write(f"B: {at}\n\n")

                out.write("  Alignment B→A:\n")
                bq = pm["b_to_a"]["aligned_query"]
                bt = pm["b_to_a"]["aligned_target"]
                mid2 = pretty_midline(bq, bt)
                out.write(f"A: {bq}\n")
                out.write(f"   {mid2}\n")
                out.write(f"B: {bt}\n\n")



def parse_args():
    parser = argparse.ArgumentParser(
        description="Agrupamento de satDNAs por similaridade circular com identidade normalizada e cobertura."
    )
    parser.add_argument("--fasta", required=True, help="FASTA de entrada")
    parser.add_argument("--out-prefix", required=True, help="Prefixo de saída")
    parser.add_argument("--id", type=float, default=0.80, help="Cutoff mínimo de identidade recíproca normalizada")
    parser.add_argument("--cov", type=float, default=0.80, help="Cutoff mínimo de cobertura recíproca")
    parser.add_argument("--match", type=int, default=2, help="Score de match")
    parser.add_argument("--mismatch", type=int, default=-1, help="Penalty de mismatch")
    parser.add_argument("--gap", type=int, default=-2, help="Penalty de gap")
    parser.add_argument("--min-len", type=int, default=1, help="Comprimento mínimo para manter sequência")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.fasta):
        eprint(f"ERRO: arquivo FASTA não encontrado: {args.fasta}")
        sys.exit(1)

    if not (0.0 <= args.id <= 1.0):
        eprint("ERRO: --id deve estar entre 0 e 1.")
        sys.exit(1)

    if not (0.0 <= args.cov <= 1.0):
        eprint("ERRO: --cov deve estar entre 0 e 1.")
        sys.exit(1)

    records = read_fasta(args.fasta)
    if not records:
        eprint("ERRO: FASTA vazio ou inválido.")
        sys.exit(1)

    seq_dict = {}
    for rid, seq in records:
        if len(seq) < args.min_len:
            continue
        if rid in seq_dict:
            eprint(f"ERRO: ID duplicado no FASTA: {rid}")
            sys.exit(1)
        seq_dict[rid] = seq

    if not seq_dict:
        eprint("ERRO: nenhuma sequência restante após filtro --min-len.")
        sys.exit(1)

    ids = sorted(seq_dict.keys())
    uf = UnionFind(ids)

    total_pairs_evaluated = 0
    total_unions = 0

    pair_metrics_map = {}
    used_edges = []
    pair_rows = []

    eprint(f"[INFO] Sequências válidas: {len(ids)}")
    eprint("[INFO] Iniciando comparações par-a-par...")

    for idx, (id_a, id_b) in enumerate(combinations(ids, 2), start=1):
        seq_a = seq_dict[id_a]
        seq_b = seq_dict[id_b]

        total_pairs_evaluated += 1

        pm = reciprocal_pair_metrics(
            seq_a=seq_a,
            seq_b=seq_b,
            match_score=args.match,
            mismatch_score=args.mismatch,
            gap_score=args.gap,
        )

        key = tuple(sorted((id_a, id_b)))
        pair_metrics_map[key] = pm

        passed = (
            pm["reciprocal_identity"] >= args.id
            and pm["reciprocal_coverage"] >= args.cov
        )

        if passed:
            changed = uf.union(id_a, id_b)
            if changed:
                used_edges.append((id_a, id_b))
                total_unions += 1

        pair_rows.append({
            "id_a": id_a,
            "len_a": len(seq_a),
            "id_b": id_b,
            "len_b": len(seq_b),
            "a_to_b_identity": pm["a_to_b"]["stats"]["normalized_identity"],
            "a_to_b_coverage": pm["a_to_b"]["stats"]["coverage"],
            "a_to_b_orientation": pm["a_to_b"]["orientation"],
            "b_to_a_identity": pm["b_to_a"]["stats"]["normalized_identity"],
            "b_to_a_coverage": pm["b_to_a"]["stats"]["coverage"],
            "b_to_a_orientation": pm["b_to_a"]["orientation"],
            "reciprocal_identity": pm["reciprocal_identity"],
            "reciprocal_coverage": pm["reciprocal_coverage"],
            "passed": passed,
        })

        if idx % 100 == 0:
            eprint(f"[INFO] Comparações concluídas: {idx}")

    fam_map = defaultdict(list)
    for sid in ids:
        fam_map[uf.find(sid)].append(sid)

    families = []
    fam_counter = 1

    for root, members in sorted(fam_map.items(), key=lambda x: (-len(x[1]), sorted(x[1])[0])):
        members = sorted(members)
        rep = choose_family_representative(members, seq_dict)
        fam_name = rep
        families.append((fam_name, members))
        fam_counter += 1

    out_families_tsv = args.out_prefix + ".families.tsv"
    out_reps_fasta = args.out_prefix + ".family_reps.fasta"
    out_proof = args.out_prefix + ".proof.txt"
    out_pairwise = args.out_prefix + ".pairwise.tsv"

    write_families_tsv(families, seq_dict, out_families_tsv)

    rep_records = []
    for fam_name, members in families:
        rep = choose_family_representative(members, seq_dict)
        rep_records.append((fam_name, seq_dict[rep]))
    write_fasta(rep_records, out_reps_fasta)

    write_pairwise_tsv(pair_rows, out_pairwise)

    write_proof_report(
        out_path=out_proof,
        families=families,
        seq_dict=seq_dict,
        used_edges=used_edges,
        pair_metrics_map=pair_metrics_map,
        args=args,
        total_pairs_evaluated=total_pairs_evaluated,
        total_unions=total_unions,
    )

    eprint("[INFO] Finalizado.")
    eprint(f"[INFO] Families TSV:      {out_families_tsv}")
    eprint(f"[INFO] Family reps FASTA: {out_reps_fasta}")
    eprint(f"[INFO] Proof report:      {out_proof}")
    eprint(f"[INFO] Pairwise TSV:      {out_pairwise}")


if __name__ == "__main__":
    main()
