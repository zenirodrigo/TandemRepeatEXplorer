#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    """
    Semiglobal A contra B2 retornando também as coordenadas no alvo B2.

    ref_start/ref_end são coordenadas 0-based half-open em B2 para a região do
    alvo usada no alinhamento de A. Isso permite recuperar a rotação circular.
    """
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
    """
    Alinha reference contra member+member e contra revcomp(member)+revcomp(member).

    Retorna:
      identity, orientation, aln_reference, aln_member, rotation_start, rotated_member

    rotated_member é o monômero inteiro do membro, orientado e rotacionado para
    iniciar no frame que melhor alinha contra o primeiro nucleotídeo do representante.
    """
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
    """Rotaciona uma sequência circular para iniciar em `start` (0-based)."""
    if not seq:
        return seq
    start %= len(seq)
    return seq[start:] + seq[:start]


def score_pair(a: str, b: str) -> int:
    return MATCH_SCORE if a == b else MISMATCH_SCORE


def nw_global_align(A: str, B: str) -> Tuple[int, str, str]:
    """Needleman-Wunsch global com gap linear. Usado após corrigir o frame circular."""
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
    """
    Alinha a representante completa contra qualquer sub-região do membro dimerizado.

    A representante é penalizada ponta-a-ponta. As pontas livres do membro dimerizado
    permitem encontrar o melhor frame circular sem forçar alinhamento contra o dímero todo.
    """
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
    """
    Encontra o melhor frame circular do membro contra a primeira sequência da família.

    1. Testa forward e reverse-complement.
    2. Dimeriza temporariamente o membro orientado.
    3. Alinha a representante contra o dímero para obter o frame.
    4. Rotaciona o membro de volta para um único monômero.
    5. Faz alinhamento global representante × monômero rotacionado, com gaps.
    """
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
    """
    Converte um alinhamento par-a-par para slots por coordenada da representante.

    insertions[-1] = inserções antes da primeira base da representante.
    insertions[p] = inserções depois da posição p da representante.
    bases[p] = base/gap do membro alinhada à posição p da representante.
    """
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
    """
    Junta os alinhamentos par-a-par contra a mesma representante em um alinhamento
    ancorado. Pode ficar largo em famílias divergentes, mas preserva a representante
    como sistema de coordenadas comum.
    """
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
    """
    Escolhe o representante/medoide da família.

    Medoide = membro com maior identidade recíproca média contra os demais.
    Desempates: maior comprimento, depois ID lexicograficamente menor.
    """
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
        "Threshold mínimo de identidade (ex: 0.80 = 80%) [0.80]: "
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
    """
    Complete-linkage real para fusão de famílias.

    A fusão só é permitida se cada membro de A passa o threshold contra cada
    membro de B. Não excluímos i/j aqui. O cache impede recalcular pares já
    avaliados.
    """
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
    """
    Gera arquivos por família usando a lógica definida com a Tarja/Rodrigo:

    - A primeira sequência da família no FASTA de entrada é a representante fixa.
    - A representante é tratada como monômero.
    - Cada outro membro é temporariamente dimerizado, alinhado contra a representante,
      orientado/rotacionado para o melhor frame circular e depois alinhado contra a
      representante com gaps.

    Outputs por família com pelo menos 2 membros:
      1. *.members.fasta
         FASTA bruto da família, na ordem original do FASTA de entrada.

      2. *.frame_corrected_monomers.fasta
         Monômeros orientados/rotacionados para o frame da primeira sequência, sem gaps.

      3. *.pairwise_to_first_representative.fasta
         Output principal para inspeção: pares representante × membro, com gaps.

      4. *.firstseq_reference_anchored_alignment.fasta
         Alinhamento ancorado na representante, fundindo os alinhamentos par-a-par.

      5. *.alignment_summary.tsv
         Metadados de orientação, rotação e identidade.

    Famílias singleton são ignoradas nessa etapa para poupar tempo e espaço.
    """
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


def connected_components_from_edges(nodes: List[str], edges: List[Tuple[str, str]]) -> List[List[str]]:
    """
    Componentes conectados simples para a auditoria pós-clustering.

    Importante: estes grupos NÃO substituem as famílias strict complete-linkage.
    Eles servem apenas para apontar representantes de famílias diferentes que ainda
    têm similaridade circular recíproca acima do threshold.
    """
    adj: Dict[str, Set[str]] = {node: set() for node in nodes}
    for a, b in edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    visited: Set[str] = set()
    components: List[List[str]] = []

    for node in sorted(adj):
        if node in visited:
            continue
        stack = [node]
        visited.add(node)
        comp = []
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj.get(x, set()):
                if y not in visited:
                    visited.add(y)
                    stack.append(y)
        if len(comp) > 1:
            components.append(sorted(comp))

    components.sort(key=lambda c: (-len(c), c[0]))
    return components


def write_representative_redundancy_audit(
    base: str,
    pct: int,
    sorted_families: List[Tuple[int, List[int]]],
    ids: List[str],
    seqs: List[str],
    lens: List[int],
    identity_threshold: float,
) -> Tuple[str, str, str, str, int, int]:
    """
    Auditoria pós-clustering no nível de representantes de famílias.

    Por que existe:
      O clustering principal é strict complete-linkage. Portanto, duas famílias
      podem permanecer separadas mesmo quando seus representantes parecem similares,
      caso algum membro de uma família falhe contra algum membro da outra.

      Além disso, o pré-filtro por k-mers pode deixar de enviar certos pares para
      alinhamento durante o clustering principal. Esta auditoria evita esse ponto:
      ela faz all-vs-all direto entre representantes de famílias, sem pré-filtro.

    O que NÃO faz:
      - Não funde famílias automaticamente.
      - Não altera family_reps.fasta nem families.tsv.
      - Não substitui o resultado strict complete-linkage.

    Saídas:
      <base>.idXX.rep_redundancy.tsv
      <base>.idXX.rep_redundancy.groups.tsv
      <base>.idXX.rep_redundancy.groups.txt
      <base>.idXX.rep_redundancy.alignments.fasta
    """
    out_pairs = f"{base}.id{pct}.rep_redundancy.tsv"
    out_groups_tsv = f"{base}.id{pct}.rep_redundancy.groups.tsv"
    out_groups_txt = f"{base}.id{pct}.rep_redundancy.groups.txt"
    out_alignments = f"{base}.id{pct}.rep_redundancy.alignments.fasta"

    # Mantém a mesma regra usada para family_reps.fasta e families.tsv.
    reps_info = []
    for family_rank, (_fid, members) in enumerate(sorted_families, 1):
        rep = min(members, key=lambda x: ids[x])
        reps_info.append({
            "family_rank": family_rank,
            "family_id": ids[rep],
            "rep_index": rep,
            "rep_id": ids[rep],
            "rep_len": lens[rep],
            "family_size": len(members),
        })

    nodes = [d["family_id"] for d in reps_info]
    edges: List[Tuple[str, str]] = []
    hit_rows = []
    alignment_records: List[Tuple[str, str]] = []

    total_pairs = len(reps_info) * (len(reps_info) - 1) // 2
    checked = 0
    t_start = time.time()

    print("\nAuditoria pós-clustering: all-vs-all entre representantes de famílias...")
    print(f"Representantes: {len(reps_info)} | pares a testar: {total_pairs}")

    for i in range(len(reps_info)):
        a = reps_info[i]
        ai = a["rep_index"]
        for j in range(i + 1, len(reps_info)):
            b = reps_info[j]
            bi = b["rep_index"]
            checked += 1

            if checked % 5000 == 0:
                elapsed = time.time() - t_start
                print(f"  Rep-audit: {checked}/{total_pairs} pares | hits={len(hit_rows)} | elapsed={elapsed:.1f}s")

            id_min, id_a_to_b, rel_a_to_b, id_b_to_a, rel_b_to_a = reciprocal_identity_only(seqs[ai], seqs[bi])

            if id_min + 1e-12 >= identity_threshold:
                edges.append((a["family_id"], b["family_id"]))
                hit_rows.append({
                    "family_A_rank": a["family_rank"],
                    "family_B_rank": b["family_rank"],
                    "family_A": a["family_id"],
                    "family_B": b["family_id"],
                    "rep_A": a["rep_id"],
                    "rep_B": b["rep_id"],
                    "rep_A_len": a["rep_len"],
                    "rep_B_len": b["rep_len"],
                    "family_A_size": a["family_size"],
                    "family_B_size": b["family_size"],
                    "reciprocal_identity_min": id_min,
                    "A_to_B_identity": id_a_to_b,
                    "A_to_B_orientation": rel_a_to_b,
                    "B_to_A_identity": id_b_to_a,
                    "B_to_A_orientation": rel_b_to_a,
                })

                id_ab, rel_ab, alnA_ab, alnB_ab = best_direction_alignment(seqs[ai], seqs[bi])
                id_ba, rel_ba, alnB_ba, alnA_ba = best_direction_alignment(seqs[bi], seqs[ai])
                pair_id = f"{sanitize_filename(a['rep_id'])}__VS__{sanitize_filename(b['rep_id'])}"
                alignment_records.append((
                    f"{pair_id}|A_to_B|A={a['rep_id']}|B={b['rep_id']}|identity={id_ab:.6f}|orientation={rel_ab}|role=A",
                    alnA_ab,
                ))
                alignment_records.append((
                    f"{pair_id}|A_to_B|A={a['rep_id']}|B={b['rep_id']}|identity={id_ab:.6f}|orientation={rel_ab}|role=B_aligned_to_A",
                    alnB_ab,
                ))
                alignment_records.append((
                    f"{pair_id}|B_to_A|A={b['rep_id']}|B={a['rep_id']}|identity={id_ba:.6f}|orientation={rel_ba}|role=B",
                    alnB_ba,
                ))
                alignment_records.append((
                    f"{pair_id}|B_to_A|A={b['rep_id']}|B={a['rep_id']}|identity={id_ba:.6f}|orientation={rel_ba}|role=A_aligned_to_B",
                    alnA_ba,
                ))

    hit_rows.sort(key=lambda d: (d["reciprocal_identity_min"], d["A_to_B_identity"], d["B_to_A_identity"]), reverse=True)

    with open(out_pairs, "w", encoding="utf-8") as out:
        out.write(
            "family_A_rank\tfamily_B_rank\tfamily_A\tfamily_B\trep_A\trep_B\t"
            "rep_A_len\trep_B_len\tfamily_A_size\tfamily_B_size\t"
            "reciprocal_identity_min\tA_to_B_identity\tA_to_B_orientation\t"
            "B_to_A_identity\tB_to_A_orientation\n"
        )
        for row in hit_rows:
            out.write(
                f"{row['family_A_rank']}\t{row['family_B_rank']}\t{row['family_A']}\t{row['family_B']}\t"
                f"{row['rep_A']}\t{row['rep_B']}\t{row['rep_A_len']}\t{row['rep_B_len']}\t"
                f"{row['family_A_size']}\t{row['family_B_size']}\t"
                f"{row['reciprocal_identity_min']:.6f}\t{row['A_to_B_identity']:.6f}\t{row['A_to_B_orientation']}\t"
                f"{row['B_to_A_identity']:.6f}\t{row['B_to_A_orientation']}\n"
            )

    components = connected_components_from_edges(nodes, edges)
    comp_index: Dict[str, int] = {}
    for group_id, comp in enumerate(components, 1):
        for fam in comp:
            comp_index[fam] = group_id

    info_by_family = {d["family_id"]: d for d in reps_info}

    with open(out_groups_tsv, "w", encoding="utf-8") as out:
        out.write("rep_redundancy_group\tfamily_id\tfamily_rank\trep_id\trep_len\tfamily_size\n")
        for group_id, comp in enumerate(components, 1):
            for fam in comp:
                d = info_by_family[fam]
                out.write(
                    f"RepRedundancyGroup_{group_id:06d}\t{fam}\t{d['family_rank']}\t"
                    f"{d['rep_id']}\t{d['rep_len']}\t{d['family_size']}\n"
                )

    with open(out_groups_txt, "w", encoding="utf-8") as out:
        out.write("# Representative-level redundancy audit\n")
        out.write("# These groups are NOT automatic family merges.\n")
        out.write("# They indicate family representatives with circular reciprocal identity >= threshold.\n")
        out.write("# Use this file for manual curation or as candidate superfamily/redundancy evidence.\n\n")
        out.write(f"# Threshold: {identity_threshold}\n")
        out.write(f"# Representatives tested: {len(reps_info)}\n")
        out.write(f"# Pairs tested: {total_pairs}\n")
        out.write(f"# Redundant representative pairs: {len(hit_rows)}\n")
        out.write(f"# Redundancy groups: {len(components)}\n\n")

        if not components:
            out.write("No representative-level redundancy groups detected.\n")
        else:
            hits_by_group: Dict[int, List[dict]] = defaultdict(list)
            for row in hit_rows:
                gid = comp_index.get(row["family_A"]) or comp_index.get(row["family_B"])
                if gid is not None:
                    hits_by_group[gid].append(row)

            for group_id, comp in enumerate(components, 1):
                out.write("=" * 100 + "\n")
                out.write(f"RepRedundancyGroup-{group_id} (families={len(comp)}):\n")
                for fam in comp:
                    d = info_by_family[fam]
                    out.write(
                        f"  Family: {fam} | rank={d['family_rank']} | rep={d['rep_id']} | "
                        f"rep_len={d['rep_len']} | family_size={d['family_size']}\n"
                    )
                out.write("\nEvidence pairs:\n")
                for row in hits_by_group.get(group_id, []):
                    out.write(
                        f"  {row['family_A']} <-> {row['family_B']} | "
                        f"min_id={row['reciprocal_identity_min']:.3f} | "
                        f"A_to_B={row['A_to_B_identity']:.3f} ({row['A_to_B_orientation']}) | "
                        f"B_to_A={row['B_to_A_identity']:.3f} ({row['B_to_A_orientation']}) | "
                        f"len={row['rep_A_len']}:{row['rep_B_len']} | "
                        f"family_size={row['family_A_size']}:{row['family_B_size']}\n"
                    )
                out.write("\n")

    write_pairwise_alignment_fasta(out_alignments, alignment_records)

    print(
        f"Rep-audit concluído: hits={len(hit_rows)} | grupos={len(components)} | "
        f"arquivos: {out_pairs}, {out_groups_txt}"
    )

    return out_pairs, out_groups_tsv, out_groups_txt, out_alignments, len(hit_rows), len(components)



def merge_families_by_representative_similarity(
    base: str,
    pct: int,
    sorted_families: List[Tuple[int, List[int]]],
    ids: List[str],
    seqs: List[str],
    lens: List[int],
    identity_threshold: float,
) -> Tuple[List[Tuple[int, List[int]]], str, str, str, str, int, int, int, int]:
    """
    Etapa pós-clustering que tenta resgatar famílias redundantes usando os
    representantes, MAS sem voltar ao problema de transitividade.

    Diferença importante em relação à versão anterior:
      - A versão anterior fazia componentes conectados entre representantes.
        Isso recriava o problema A~B, B~C => A+B+C, mesmo que A não passasse
        com C.
      - Esta versão usa complete-linkage também no nível dos representantes.
        Dois grupos de famílias só são fundidos se TODOS os representantes de
        um grupo passarem o threshold contra TODOS os representantes do outro.

    Regra:
      1. Calcula um representante por família strict inicial usando a mesma
         regra de family_reps.fasta: menor ID lexicográfico dentro da família.
      2. Faz all-vs-all entre representantes, sem pré-filtro por k-mer.
      3. Registra todos os pares de representantes com identidade circular
         recíproca mínima >= threshold.
      4. Tenta fundir grupos por complete-linkage entre representantes.
      5. Os outputs finais do pipeline passam a ser baseados nas famílias após
         esse merge conservador por representantes.

    Saídas adicionais:
      <base>.idXX.rep_merge_pairs.tsv
      <base>.idXX.rep_merge.groups.tsv
      <base>.idXX.rep_merge.groups.txt
      <base>.idXX.rep_merge.alignments.fasta
    """
    out_pairs = f"{base}.id{pct}.rep_merge_pairs.tsv"
    out_groups_tsv = f"{base}.id{pct}.rep_merge.groups.tsv"
    out_groups_txt = f"{base}.id{pct}.rep_merge.groups.txt"
    out_alignments = f"{base}.id{pct}.rep_merge.alignments.fasta"

    reps_info = []
    for family_rank, (fid, members) in enumerate(sorted_families, 1):
        rep = min(members, key=lambda x: ids[x])
        family_id = ids[rep]
        reps_info.append({
            "family_rank": family_rank,
            "strict_fid": fid,
            "family_id": family_id,
            "rep_index": rep,
            "rep_id": ids[rep],
            "rep_len": lens[rep],
            "family_size": len(members),
            "members": list(members),
        })

    info_by_family = {d["family_id"]: d for d in reps_info}
    pair_metrics: Dict[Tuple[str, str], dict] = {}
    hit_rows: List[dict] = []
    alignment_records: List[Tuple[str, str]] = []

    def fam_pair_key(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    total_pairs = len(reps_info) * (len(reps_info) - 1) // 2
    checked = 0
    t_start = time.time()

    print("\nEtapa pós-clustering: auditoria/merge por representantes com complete-linkage...")
    print(f"Representantes: {len(reps_info)} | pares a testar: {total_pairs}")

    for i in range(len(reps_info)):
        a = reps_info[i]
        ai = a["rep_index"]
        for j in range(i + 1, len(reps_info)):
            b = reps_info[j]
            bi = b["rep_index"]
            checked += 1

            if checked % 5000 == 0:
                elapsed = time.time() - t_start
                print(
                    f"  Rep-audit: {checked}/{total_pairs} pares | "
                    f"candidate_hits={len(hit_rows)} | elapsed={elapsed:.1f}s"
                )

            id_min, id_a_to_b, rel_a_to_b, id_b_to_a, rel_b_to_a = reciprocal_identity_only(seqs[ai], seqs[bi])

            row = {
                "family_A_rank": a["family_rank"],
                "family_B_rank": b["family_rank"],
                "family_A": a["family_id"],
                "family_B": b["family_id"],
                "rep_A": a["rep_id"],
                "rep_B": b["rep_id"],
                "rep_A_len": a["rep_len"],
                "rep_B_len": b["rep_len"],
                "family_A_size": a["family_size"],
                "family_B_size": b["family_size"],
                "reciprocal_identity_min": id_min,
                "A_to_B_identity": id_a_to_b,
                "A_to_B_orientation": rel_a_to_b,
                "B_to_A_identity": id_b_to_a,
                "B_to_A_orientation": rel_b_to_a,
                "passes_threshold": id_min + 1e-12 >= identity_threshold,
            }
            pair_metrics[fam_pair_key(a["family_id"], b["family_id"])] = row

            if row["passes_threshold"]:
                hit_rows.append(row)

                id_ab, rel_ab, alnA_ab, alnB_ab = best_direction_alignment(seqs[ai], seqs[bi])
                id_ba, rel_ba, alnB_ba, alnA_ba = best_direction_alignment(seqs[bi], seqs[ai])
                pair_id = f"{sanitize_filename(a['rep_id'])}__VS__{sanitize_filename(b['rep_id'])}"
                alignment_records.append((
                    f"{pair_id}|A_to_B|A={a['rep_id']}|B={b['rep_id']}|identity={id_ab:.6f}|orientation={rel_ab}|role=A",
                    alnA_ab,
                ))
                alignment_records.append((
                    f"{pair_id}|A_to_B|A={a['rep_id']}|B={b['rep_id']}|identity={id_ab:.6f}|orientation={rel_ab}|role=B_aligned_to_A",
                    alnB_ab,
                ))
                alignment_records.append((
                    f"{pair_id}|B_to_A|A={b['rep_id']}|B={a['rep_id']}|identity={id_ba:.6f}|orientation={rel_ba}|role=B",
                    alnB_ba,
                ))
                alignment_records.append((
                    f"{pair_id}|B_to_A|A={b['rep_id']}|B={a['rep_id']}|identity={id_ba:.6f}|orientation={rel_ba}|role=A_aligned_to_B",
                    alnA_ba,
                ))

    hit_rows.sort(
        key=lambda d: (d["reciprocal_identity_min"], d["A_to_B_identity"], d["B_to_A_identity"]),
        reverse=True,
    )

    with open(out_pairs, "w", encoding="utf-8") as out:
        out.write(
            "family_A_rank\tfamily_B_rank\tfamily_A\tfamily_B\trep_A\trep_B\t"
            "rep_A_len\trep_B_len\tfamily_A_size\tfamily_B_size\t"
            "reciprocal_identity_min\tA_to_B_identity\tA_to_B_orientation\t"
            "B_to_A_identity\tB_to_A_orientation\tpasses_threshold\n"
        )
        for row in hit_rows:
            out.write(
                f"{row['family_A_rank']}\t{row['family_B_rank']}\t{row['family_A']}\t{row['family_B']}\t"
                f"{row['rep_A']}\t{row['rep_B']}\t{row['rep_A_len']}\t{row['rep_B_len']}\t"
                f"{row['family_A_size']}\t{row['family_B_size']}\t"
                f"{row['reciprocal_identity_min']:.6f}\t{row['A_to_B_identity']:.6f}\t{row['A_to_B_orientation']}\t"
                f"{row['B_to_A_identity']:.6f}\t{row['B_to_A_orientation']}\t1\n"
            )

    # ------------------------------------------------------------------
    # Merge conservador: complete-linkage entre representantes.
    # ------------------------------------------------------------------
    rep_families = [d["family_id"] for d in reps_info]
    groups: List[List[str]] = [[fam] for fam in rep_families]

    def find_group_index(fam: str) -> int:
        for idx, group in enumerate(groups):
            if fam in group:
                return idx
        raise KeyError(fam)

    def groups_pass_complete_linkage(group_a: List[str], group_b: List[str]) -> bool:
        for fam_a in group_a:
            for fam_b in group_b:
                if fam_a == fam_b:
                    continue
                row = pair_metrics.get(fam_pair_key(fam_a, fam_b))
                if row is None or not row["passes_threshold"]:
                    return False
        return True

    accepted_merge_edges: List[Tuple[str, str]] = []
    rejected_transitive_edges: List[Tuple[str, str]] = []

    # Processa os pares mais fortes primeiro. A diferença é que um par só funde
    # grupos se todos os demais representantes dos dois grupos também passarem.
    for row in hit_rows:
        fam_a = row["family_A"]
        fam_b = row["family_B"]
        ga = find_group_index(fam_a)
        gb = find_group_index(fam_b)
        if ga == gb:
            continue

        group_a = groups[ga]
        group_b = groups[gb]
        if groups_pass_complete_linkage(group_a, group_b):
            accepted_merge_edges.append((fam_a, fam_b))
            merged = sorted(set(group_a + group_b), key=lambda fam: info_by_family[fam]["family_rank"])
            # remove índices em ordem decrescente para não bagunçar a lista
            for idx in sorted([ga, gb], reverse=True):
                del groups[idx]
            groups.append(merged)
        else:
            rejected_transitive_edges.append((fam_a, fam_b))

    components = [g for g in groups if len(g) > 1]
    components.sort(key=lambda comp: (-len(comp), min(info_by_family[f]["family_rank"] for f in comp)))

    comp_index: Dict[str, int] = {}
    for group_id, comp in enumerate(components, 1):
        for fam in comp:
            comp_index[fam] = group_id

    with open(out_groups_tsv, "w", encoding="utf-8") as out:
        out.write(
            "rep_merge_group\told_family_id\told_family_rank\told_rep_id\told_rep_len\t"
            "old_family_size\tnew_family_id\tnew_family_size\n"
        )
        for group_id, comp in enumerate(components, 1):
            merged_members: List[int] = []
            for fam in comp:
                merged_members.extend(info_by_family[fam]["members"])
            merged_members = sorted(set(merged_members))
            new_rep = min(merged_members, key=lambda x: ids[x])
            new_family_id = ids[new_rep]
            new_family_size = len(merged_members)
            for fam in comp:
                d = info_by_family[fam]
                out.write(
                    f"RepMergeGroup_{group_id:06d}\t{fam}\t{d['family_rank']}\t"
                    f"{d['rep_id']}\t{d['rep_len']}\t{d['family_size']}\t"
                    f"{new_family_id}\t{new_family_size}\n"
                )

    with open(out_groups_txt, "w", encoding="utf-8") as out:
        out.write("# Representative-level conservative family merging\n")
        out.write("# Families listed in the same RepMergeGroup were automatically merged.\n")
        out.write("# Criterion: complete-linkage among family representatives.\n")
        out.write("# This avoids transitive chains such as A~B and B~C merging A+B+C when A~C fails.\n\n")
        out.write(f"# Threshold: {identity_threshold}\n")
        out.write(f"# Strict families before representative merge: {len(sorted_families)}\n")
        out.write(f"# Representatives tested: {len(reps_info)}\n")
        out.write(f"# Pairs tested: {total_pairs}\n")
        out.write(f"# Representative pairs passing threshold: {len(hit_rows)}\n")
        out.write(f"# Accepted complete-linkage merge edges: {len(accepted_merge_edges)}\n")
        out.write(f"# Rejected transitive/non-complete edges: {len(rejected_transitive_edges)}\n")
        out.write(f"# Merge groups accepted: {len(components)}\n\n")

        if not components:
            out.write("No representative-level complete-linkage merge groups detected. Final families remain unchanged.\n")
        else:
            hits_by_group: Dict[int, List[dict]] = defaultdict(list)
            for row in hit_rows:
                gid_a = comp_index.get(row["family_A"])
                gid_b = comp_index.get(row["family_B"])
                if gid_a is not None and gid_a == gid_b:
                    hits_by_group[gid_a].append(row)

            for group_id, comp in enumerate(components, 1):
                merged_members: List[int] = []
                for fam in comp:
                    merged_members.extend(info_by_family[fam]["members"])
                merged_members = sorted(set(merged_members))
                new_rep = min(merged_members, key=lambda x: ids[x])
                out.write("=" * 100 + "\n")
                out.write(
                    f"RepMergeGroup-{group_id} | old_families={len(comp)} | "
                    f"new_family_id={ids[new_rep]} | new_family_size={len(merged_members)}\n"
                )
                for fam in comp:
                    d = info_by_family[fam]
                    out.write(
                        f"  Old family: {fam} | rank={d['family_rank']} | rep={d['rep_id']} | "
                        f"rep_len={d['rep_len']} | family_size={d['family_size']}\n"
                    )
                out.write("\nEvidence pairs within accepted group:\n")
                for row in hits_by_group.get(group_id, []):
                    out.write(
                        f"  {row['family_A']} <-> {row['family_B']} | "
                        f"min_id={row['reciprocal_identity_min']:.3f} | "
                        f"A_to_B={row['A_to_B_identity']:.3f} ({row['A_to_B_orientation']}) | "
                        f"B_to_A={row['B_to_A_identity']:.3f} ({row['B_to_A_orientation']}) | "
                        f"len={row['rep_A_len']}:{row['rep_B_len']} | "
                        f"family_size={row['family_A_size']}:{row['family_B_size']}\n"
                    )
                out.write("\n")

    write_pairwise_alignment_fasta(out_alignments, alignment_records)

    # Aplica a fusão conservative-complete-linkage.
    used_families: Set[str] = set()
    merged_groups: Dict[str, List[int]] = {}

    for group_id, comp in enumerate(components, 1):
        merged_members: List[int] = []
        for fam in comp:
            merged_members.extend(info_by_family[fam]["members"])
            used_families.add(fam)
        merged_members = sorted(set(merged_members))
        new_rep = min(merged_members, key=lambda x: ids[x])
        merged_groups[f"merged_{group_id}_{ids[new_rep]}"] = merged_members

    for d in reps_info:
        fam = d["family_id"]
        if fam in used_families:
            continue
        merged_groups[f"single_{fam}"] = sorted(set(d["members"]))

    merged_final_families: Dict[int, List[int]] = {}
    for members in merged_groups.values():
        rep = min(members, key=lambda x: ids[x])
        merged_final_families[rep] = sorted(set(members))

    merged_sorted_families = sorted(merged_final_families.items(), key=lambda item: len(item[1]), reverse=True)

    old_family_count = len(sorted_families)
    new_family_count = len(merged_sorted_families)

    print(
        f"Rep-merge complete-linkage concluído: candidate_pairs={len(hit_rows)} | "
        f"accepted_groups={len(components)} | famílias strict={old_family_count} -> "
        f"famílias finais={new_family_count} | arquivos: {out_pairs}, {out_groups_txt}"
    )

    return (
        merged_sorted_families,
        out_pairs,
        out_groups_tsv,
        out_groups_txt,
        out_alignments,
        len(hit_rows),
        len(components),
        old_family_count,
        new_family_count,
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

    (
        sorted_families,
        out_rep_merge_pairs,
        out_rep_merge_groups_tsv,
        out_rep_merge_groups_txt,
        out_rep_merge_alignments,
        rep_merge_pairs,
        rep_merge_groups,
        strict_family_count_before_rep_merge,
        final_family_count_after_rep_merge,
    ) = merge_families_by_representative_similarity(
        base,
        pct,
        sorted_families,
        ids,
        seqs,
        lens,
        identity_threshold,
    )

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
    report.append(f"# Strict families before representative merge: {strict_family_count_before_rep_merge}")
    report.append(f"# Final families after representative merge: {final_family_count_after_rep_merge}")
    report.append(f"# Representative merge pairs: {out_rep_merge_pairs}")
    report.append(f"# Representative merge groups TSV: {out_rep_merge_groups_tsv}")
    report.append(f"# Representative merge groups TXT: {out_rep_merge_groups_txt}")
    report.append(f"# Representative merge alignments: {out_rep_merge_alignments}")
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
    print(f"Strict families before rep-merge: {strict_family_count_before_rep_merge}")
    print(f"Final families after rep-merge:   {final_family_count_after_rep_merge}")
    print(f"REP_MERGE_PAIRS:         {out_rep_merge_pairs}")
    print(f"REP_MERGE_GROUPS:        {out_rep_merge_groups_txt}")
    print(f"REP_MERGE_ALNS:          {out_rep_merge_alignments}")
    print(f"Rep merge pairs:         {rep_merge_pairs}")
    print(f"Rep merge groups:        {rep_merge_groups}")
    print(f"Family member FASTAs:    {family_files_written}")
    print(f"Frame-corrected FASTAs:  {frame_corrected_files_written}")
    print(f"Pairwise-to-first FASTAs:{pairwise_files_written}")
    print(f"Anchored alignment FASTAs:{anchored_files_written}")


if __name__ == "__main__":
    fasta_path, identity_value = ask_user()
    run(fasta_path, identity_value)

