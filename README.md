# antigen_cage

A toolkit for combining antigen and cage DNA sequences and inserting them into a GenBank expression backbone.

## Overview

The workflow has two steps:

1. **`combine.py`** — selects antigen and cage sequences from FASTA files, strips overhangs, and joins them into a combined sequence.
2. **`insert.py`** — inserts each combined sequence into a GenBank backbone at a specified site, updates all feature annotations, and translates the inserted CDS.

---

## Directory structure

```
antigen_cage/
├── antigen/          # Antigen FASTA files
├── cage/             # Cage FASTA files
├── backbone/         # GenBank backbone file(s)
├── output/           # Generated combined FASTA files
│   └── inserted/     # Generated GenBank and protein files
├── combine.py
└── insert.py
```

---

## FASTA file format

Both antigen and cage FASTA files follow the same convention:

| Position | Header format | Description |
|----------|--------------|-------------|
| 1st entry | `>A_5overhang` / `>C_5overhang` | 5′ overhang sequence |
| 2nd entry | `>A_3overhang` / `>C_3overhang` | 3′ overhang sequence |
| Remaining | `>A{ID}_{name}` / `>C{ID}_{name}` | Sequences, where ID is a unique integer |

Example:
```
>A_5overhang
GATCAGATATCGCGGCCGCTCTAGAGCCACC

>A_3overhang
GGGGAGCTGGCATATAAACTC

>A1_TH-WI22-M181I
GATCAGATATCGCGGCCGCTCTAGAGCCACC...GGGGAGCTGGCATATAAACTC

>A2_His6-TH-WI22-M181I
...
```

---

## Step 1 — Combine antigen and cage sequences

### Usage

```bash
python combine.py <antigen.fasta> <cage.fasta> A{x}C{y} [A{x}C{y} ...] [--output-dir DIR]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `antigen.fasta` | Antigen FASTA file |
| `cage.fasta` | Cage FASTA file |
| `A{x}C{y}` | One or more combinations specifying antigen ID `x` and cage ID `y` |
| `--output-dir` | Output directory (default: `output`) |

### Example

```bash
python combine.py \
    antigen/20260430_twist.fasta \
    cage/20260430_twist.fasta \
    A1C6 A2C7 A3C1 \
    --output-dir output
```

### Output

`output/combined.fasta` — FASTA file with one record per combination.
Headers follow the format `>{A_ID}{C_ID}_{antigen_name}_{cage_name}`:

```
>A1C6_TH-WI22-M181I_dn5B-hisless
ATGGATTCC...
```

### Overhang handling

| Sequence | 5′ overhang | 3′ overhang |
|----------|------------|------------|
| Antigen | Removed | Removed |
| Cage | Kept | Removed |

The combined sequence is: `antigen_core + cage_without_3overhang`.
An exception is raised if a sequence does not contain the expected overhang.

---

## Step 2 — Insert into GenBank backbone

> Requires BioPython. Run via the Singularity container:
> ```bash
> singularity exec /net/software/containers/universal.sif python insert.py ...
> ```

### Usage

```bash
singularity exec /net/software/containers/universal.sif python insert.py \
    <insertion.fasta> <backbone.gb> <insertion_site> [--output-dir DIR]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `insertion.fasta` | FASTA file with sequences to insert (e.g. output of `combine.py`) |
| `backbone.gb` | GenBank backbone file |
| `insertion_site` | 1-based position; insertion happens **after** this position |
| `--output-dir` | Output directory (default: `output`) |

### Example

```bash
singularity exec /net/software/containers/universal.sif python insert.py \
    output/combined.fasta \
    backbone/CMVR_insertto1386.gb \
    1386 \
    --output-dir output/inserted
```

### Output

| File | Description |
|------|-------------|
| `output/inserted/{ID}.gb` | One GenBank file per inserted sequence, with all downstream features shifted |
| `output/inserted/proteins.faa` | FASTA file with translated protein sequences for all insertions |

### Sequence validation

Each insertion sequence is validated and translated. Warnings are printed to stderr if:
- The sequence does not start with `ATG`
- There are premature stop codons
- The sequence does not end with a stop codon

The output is still written regardless of warnings.

---

## Full workflow example

```bash
# Step 1: combine antigen and cage sequences
python combine.py \
    antigen/20260430_twist.fasta \
    cage/20260430_twist.fasta \
    A1C6 A2C7 A3C1 \
    --output-dir output

# Step 2: insert into backbone at site 1386
singularity exec /net/software/containers/universal.sif python insert.py \
    output/combined.fasta \
    backbone/CMVR_insertto1386.gb \
    1386 \
    --output-dir output/inserted
```
