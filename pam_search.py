#!/usr/bin/env python3
"""PAM-discovery read processing for the in vitro cleavage (PASIViC) assay.

This script processes paired-end NGS reads from a randomized-PAM library to
recover the PAM sequences that remain after Cas9 cleavage/selection. Each
amplicon has the following layout (5' -> 3'):

    [ ...forward adapter... ][ protospacer (20 nt) ][ N6 / N10 randomized PAM ]
    [ downstream constant region ][ ...reverse adapter... ]

The pipeline runs in four stages per sample, writing intermediate files into
numbered directories so each step can be inspected independently:

    1-rawdata/   ->  input FASTQ files (R1 and R2)                 [input]
    2-fastq/     ->  reads oriented and trimmed to the insert      [stage 1]
    3-notarget/  ->  reads kept for PAM analysis                   [stage 2]
    4-pamsort/   ->  reads with the PAM region at an expected offset[stage 3]
    5-pamresult/ ->  extracted PAM sequences, grouped by offset    [stage 4]

Stage 1  Orient each read (R2 is reverse-complemented) and trim it to the
         region located between the sequencing adapters.
Stage 2  Keep reads that do NOT contain the intact protospacer 10-mer but DO
         contain the downstream constant 10-mer, within an expected length
         window. Reads found on the reverse strand are reverse-complemented so
         that all kept reads share one orientation.
Stage 3  Keep reads in which the downstream constant region starts at an
         expected distance from the read start (i.e. the randomized region has
         the expected length), and record them.
Stage 4  Extract the N-base PAM immediately upstream of the downstream constant
         region, writing one file per observed offset.

Only the CONFIG block below needs editing between experiments.
"""

from __future__ import annotations

import time
from pathlib import Path

# ============================== CONFIG ========================================
# Sample indices to process (used in both the FASTQ file names and outputs).
SAMPLE_INDICES = ["5", "6", "7", "8"]

# Paired-end read numbers to combine (R1 and R2).
READ_NUMBERS = [1, 2]

# gRNA-target 10-mer: positions 20..11 of the protospacer, counted from the
# PAM-proximal end (5' -> 3'). Reads containing this (or its reverse complement)
# are excluded in stage 2.
GRNA_TARGET_10MER = "gagccacatt"

# Constant 10-mer located immediately 3' of the randomized PAM (5' -> 3').
DOWNSTREAM_10MER = "ataaggtggt"

# Number of randomized (N) bases in the PAM library (e.g. 6 or 10).
N_COUNT = 10

# Sequencing adapters (5' -> 3').
ADAPTER_F = "CTTCCGATCT"
ADAPTER_R = "AGATCGGAAG"

# Length window for reads kept in stage 2 (length of the trimmed insert).
MIN_INSERT_LEN = 90   # inclusive
MAX_INSERT_LEN = 103  # inclusive

# Input / output directories.
RAW_DIR = Path("1-rawdata")
FASTQ_DIR = Path("2-fastq")
NOTARGET_DIR = Path("3-notarget")
PAMSORT_DIR = Path("4-pamsort")
PAMRESULT_DIR = Path("5-pamresult")
# =============================================================================


_COMPLEMENT = str.maketrans("ATGC", "TACG")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (non-ACGT left as-is)."""
    return seq.translate(_COMPLEMENT)[::-1]


# Pre-compute adapter reverse complements once.
ADAPTER_F_RC = reverse_complement(ADAPTER_F)
ADAPTER_R_RC = reverse_complement(ADAPTER_R)


def read_fastq_sequences(path: Path) -> list[str]:
    """Return the sequence line of every record in a FASTQ file.

    FASTQ stores four lines per record; the sequence is the 2nd line.
    """
    lines = path.read_text().splitlines()
    return [lines[4 * i + 1].strip() for i in range(len(lines) // 4)]


def trim_between_adapters(seq: str, default: int) -> str:
    """Trim ``seq`` to the insert located between the forward/reverse adapters.

    ``default`` is the fall-back coordinate used when an adapter is absent
    (0 for reads kept in the forward orientation, ``len(seq)`` for reads that
    were reverse-complemented). The two branches below therefore differ only in
    this default, matching the original two code paths for R1 and R2.
    """
    # Forward adapter: locate the insert's 5' boundary.
    if ADAPTER_F in seq:
        f_start = seq.find(ADAPTER_F) + len(ADAPTER_F)
        f_rev_end = default
    elif ADAPTER_F_RC in seq:
        f_rev_end = seq.find(ADAPTER_F_RC)
        f_start = default  # NOTE: the original script had a typo here
                           # ("adaprFpo") that left this value unset, so it
                           # carried over from the previous read. It is assigned
                           # correctly here. This only affects the rare reads in
                           # which the forward adapter appears reverse-
                           # complemented; results are otherwise identical.
    else:
        f_start = f_rev_end = default

    # Reverse adapter: locate the insert's 3' boundary.
    if ADAPTER_R in seq:
        r_end = seq.find(ADAPTER_R)
        r_rev_start = default
    elif ADAPTER_R_RC in seq:
        r_rev_start = seq.find(ADAPTER_R_RC) + len(ADAPTER_R_RC)
        r_end = default
    else:
        r_end = r_rev_start = default

    if f_start < r_end:
        return seq[f_start:r_end]
    if r_rev_start < f_rev_end:
        return seq[r_rev_start:f_rev_end]
    return seq


def stage1_orient_and_trim(index: str) -> None:
    """Orient and trim R1/R2 reads, writing inserts to ``2-fastq/``."""
    out_path = FASTQ_DIR / f"{index}_fastq.txt"
    with out_path.open("w") as out:
        for read_num in READ_NUMBERS:
            fastq_path = RAW_DIR / f"{index}_S{index}_L001_R{read_num}_001.fastq"
            for seq in read_fastq_sequences(fastq_path):
                if read_num == 1:
                    oriented, default = seq, 0
                else:
                    oriented = reverse_complement(seq)
                    default = len(oriented)
                out.write(trim_between_adapters(oriented, default) + "\n")


def stage2_select_reads(index: str) -> None:
    """Select reads for PAM analysis, writing them to ``3-notarget/``.

    Keep reads that lack the intact protospacer 10-mer but contain the
    downstream constant 10-mer within the expected length window. Reads on the
    reverse strand are reverse-complemented so all output shares one orientation.
    """
    target = GRNA_TARGET_10MER.upper()
    target_rc = reverse_complement(target)
    downstream = DOWNSTREAM_10MER.upper()
    downstream_rc = reverse_complement(downstream)

    in_path = FASTQ_DIR / f"{index}_fastq.txt"
    out_path = NOTARGET_DIR / f"{index}_fastq.txt"
    with in_path.open() as fin, out_path.open("w") as out:
        for line in fin:
            seq = line.strip()
            if target in seq or target_rc in seq:
                continue
            if not (MIN_INSERT_LEN <= len(seq) <= MAX_INSERT_LEN):
                continue
            if downstream in seq:
                out.write(seq + "\n")
            elif downstream_rc in seq:
                out.write(reverse_complement(seq) + "\n")


def stage3_and_4_extract_pam(index: str) -> None:
    """Sort reads by PAM offset and extract PAM sequences.

    Stage 3 keeps reads in which the downstream region starts at an expected
    distance from the read start (``N_COUNT`` .. ``N_COUNT + 9``) and records
    them in ``4-pamsort/``. Stage 4 extracts the ``N_COUNT``-base PAM lying
    immediately upstream of the downstream region, writing one file per offset
    into ``5-pamresult/``.
    """
    downstream = DOWNSTREAM_10MER.upper()

    # Stage 3: keep reads whose PAM region has an expected offset.
    sorted_reads: list[str] = []
    in_path = NOTARGET_DIR / f"{index}_fastq.txt"
    sort_path = PAMSORT_DIR / f"{index}_fastq.txt"
    combined_path = PAMRESULT_DIR / f"{index}_Length_3-5_fastq.txt"

    with in_path.open() as fin, sort_path.open("w") as sort_out, \
            combined_path.open("w") as combined_out:
        for line in fin:
            seq = line.strip()
            offset = seq.find(downstream)
            if N_COUNT <= offset <= N_COUNT + 9:
                sort_out.write(seq + "\n")
                sorted_reads.append(seq)
                # Convenience file combining the central offsets.
                if N_COUNT + 2 <= offset <= N_COUNT + 4:
                    combined_out.write(seq[offset - N_COUNT:offset] + "\n")

    # Stage 4: one PAM file per observed offset.
    for offset in range(N_COUNT + 1, N_COUNT + 10):
        length_label = offset - N_COUNT
        result_path = PAMRESULT_DIR / f"{index}_Length_{length_label}_fastq.txt"
        with result_path.open("w") as result_out:
            for seq in sorted_reads:
                if seq.find(downstream) == offset:
                    result_out.write(seq[offset - N_COUNT:offset] + "\n")


def main() -> None:
    start = time.time()

    for directory in (FASTQ_DIR, NOTARGET_DIR, PAMSORT_DIR, PAMRESULT_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    print(f"gRNA target 10-mer      : {GRNA_TARGET_10MER.upper()}")
    print(f"  reverse complement    : {reverse_complement(GRNA_TARGET_10MER.upper())}")
    print(f"downstream 10-mer       : {DOWNSTREAM_10MER.upper()}")
    print(f"  reverse complement    : {reverse_complement(DOWNSTREAM_10MER.upper())}")
    print(f"randomized PAM length N  : {N_COUNT}\n")

    for index in SAMPLE_INDICES:
        print(f"processing sample {index} ...")
        stage1_orient_and_trim(index)
        stage2_select_reads(index)
        stage3_and_4_extract_pam(index)

    print(f"\ndone in {time.time() - start:.1f} sec")


if __name__ == "__main__":
    main()
