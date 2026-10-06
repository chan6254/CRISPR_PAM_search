# CRISPR_PAM_search
PAM identification pipeline for RsCas9, a type II-C CRISPR-Cas9 ortholog

# PAM_search
A read-processing pipeline for identifying protospacer-adjacent motifs (PAMs) from
randomized-PAM in vitro cleavage and selection assays.

The pipeline is nuclease-agnostic. It takes paired-end sequencing of a PAM library that has
been challenged with any RNA-guided nuclease and recovers the PAM sequences of the molecules
that were cleaved. The protospacer, the constant region flanking the randomized positions, the number
of randomized positions, the sequencing adapters and the insert-length window are all set in a
single configuration block, so the same script can be applied to a different nuclease, a
different target or a different library design without changing the code.

It was written for, and used in:

> Kang CY, Ye S, Lee EA, An S, *et al.* Characterization of a new functional CRISPR-Cas9
> from an unidentified bacterial strain.

where it was used to determine the NNRRAY PAM of RsCas9, a type II-C CRISPR-Cas9 ortholog
from an environmental *Rhodobacter sphaeroides* isolate.

---

## What the pipeline does

PAM preferences are determined by an in vitro cleavage and selection assay on a
randomized-PAM library. Each amplicon has the layout

```
[ forward adapter ][ protospacer (20 nt) ][ N6 / N10 randomized PAM ][ constant region ][ reverse adapter ]
```

Molecules carrying a PAM that the nuclease recognizes are cleaved and lose the intact
protospacer. The pipeline recovers, from paired-end sequencing of the surviving library, the
PAM sequences associated with cleavage, which are then summarized as a sequence logo.

The analysis runs in four stages per sample, each writing into a numbered directory so that
intermediate output can be inspected:

| Stage | Output directory | What it does |
|---|---|---|
| input | `1-rawdata/` | paired-end FASTQ files (R1, R2) |
| 1 | `2-fastq/` | orients each read (R2 reverse-complemented) and trims to the insert between the sequencing adapters |
| 2 | `3-notarget/` | keeps reads that lack the intact protospacer 10-mer but retain the downstream constant 10-mer, within an expected insert-length window; reverse-strand reads are reverse-complemented so all output shares one orientation |
| 3 | `4-pamsort/` | keeps reads in which the constant region begins at an expected distance from the read start, i.e. the randomized region has the expected length |
| 4 | `5-pamresult/` | extracts the N-base PAM immediately upstream of the constant region, one file per observed offset |

## Usage

```bash
mkdir -p 1-rawdata                 # place FASTQ files here
python3 pam_search.py
```

Requires Python 3.9 or later. No third-party packages.

FASTQ files are expected to be named `{index}_S{index}_L001_R{1,2}_001.fastq`.

## Configuration

Everything that changes between experiments is collected in the `CONFIG` block at the top of
`pam_search.py`:

| Parameter | Meaning |
|---|---|
| `SAMPLE_INDICES` | sample indices to process |
| `GRNA_TARGET_10MER` | protospacer positions 20–11 counted from the PAM-proximal end |
| `DOWNSTREAM_10MER` | constant 10-mer immediately 3′ of the randomized PAM |
| `N_COUNT` | number of randomized positions (6 or 10) |
| `ADAPTER_F`, `ADAPTER_R` | sequencing adapters |
| `MIN_INSERT_LEN`, `MAX_INSERT_LEN` | insert-length window used in stage 2 |

## Data availability

The sequencing data analyzed in the accompanying paper are deposited in the NCBI Sequence
Read Archive under BioProject
[PRJNA1429032](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1429032), and can be used to
reproduce the published PAM analysis with the configuration shown above.
