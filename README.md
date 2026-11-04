# antigen_cage

A toolkit for combining antigen and cage DNA sequences and inserting them into a GenBank expression backbone.

## Overview

The workflow has two steps:

1. **`combine.py`** — prepares the sequence(s) to insert. In the default **two-part mode** it selects antigen and cage sequences from two FASTA files, strips overhangs, and joins them. In **single-part mode** (`--single`) it takes one FASTA file, strips both overhangs off each selected sequence, and writes them out without joining anything.
2. **`insert.py`** — inserts each combined sequence into a GenBank backbone at a specified site and updates all feature annotations. Requires BioPython; must be run via the Singularity container (see below).

---

## Directory structure

```
antigen_cage/
├── antigen/          # Antigen FASTA files
├── cage/             # Cage FASTA files
├── full_insertion/   # Full-length inserts ordered as a single part (--single)
├── backbone/         # GenBank backbone file(s)
├── output/           # Generated combined FASTA files
│   └── inserted/     # Generated GenBank files (one per combination)
├── combine.py        # No dependencies — runs with plain python
└── insert.py         # Requires BioPython — run via Singularity
```

## Requirements

| Script | Runtime |
|--------|---------|
| `combine.py` | Plain `python` (no dependencies) |
| `insert.py` | `singularity exec /net/software/containers/universal.sif python` (BioPython required) |

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

## Step 1 — Combine antigen and cage sequences (two-part mode)

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
| `--output-name` | Output file name (default: `combined.fasta`) |

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

The antigen's 3′ overhang and the cage's 5′ overhang are the same linker sequence,
which is why it is kept exactly once — on the cage side — when the two parts are joined.

---

## Step 1 (alternative) — Single-part insert

Use this when the insert is **one part** rather than an antigen joined to a cage —
typically a full-length antigen–cage fusion ordered as a single fragment (see
`full_insertion/`). There is no junction to preserve, so both overhangs are
trimmed and the bare core is written out.

### Usage

```bash
python combine.py --single <part.fasta> [<ID> ...] [--output-dir DIR] [--output-name NAME]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--single` | Enables single-part mode |
| `part.fasta` | FASTA file holding the part(s): a 5′ overhang entry, a 3′ overhang entry, then the sequences |
| `<ID>` | Optional. Which sequences to take — **omit to take every sequence in the file**. Accepts a full entry name (`38g_floppy_IA8`), a prefixed ID (`A3` → `A3_38g-Yoan`), or a bare number (`7` → `C7_dn5B`) |
| `--output-dir` | Output directory (default: `output`) |
| `--output-name` | Output file name (default: `combined.fasta`) |

The two overhang entries are matched on their `_5overhang` / `_3overhang` suffix, and
**their prefixes need not agree**. A full-length insert spans the antigen's 5′ end and
the cage's 3′ end, so it legitimately carries `A_5overhang` together with `C_3overhang`.

### Example

```bash
# every sequence in the file
python combine.py --single full_insertion/20261104_twist.fasta \
    --output-dir output --output-name full_insertion.fasta

# just two of them
python combine.py --single full_insertion/20261104_twist.fasta \
    38g_floppy_IA8 'SOSIP_SC(8xGS)_EXT-1hep_A3C3' --output-dir output
```

### Output

The same format as two-part mode, so **Step 2 runs unchanged**. Each entry keeps its
own name as the header, which becomes the `.gb` filename in Step 2:

```
>38g_floppy_A3C3
ATGGAC...
```

Quote any name containing shell metacharacters, such as the parentheses in
`SOSIP_SC(8xGS)_EXT-1hep_A3C3`.

### Overhang handling

| Sequence | 5′ overhang | 3′ overhang |
|----------|------------|------------|
| Part | Removed | Removed |

---

## Step 2 — Insert into GenBank backbone

> **Requires BioPython.** BioPython is not available in the default Python environment and must be accessed via the Singularity container:
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
| `insertion_site` | 1-based position; insertion starts **at** this position |
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

### Sequence validation

Warnings are printed to stderr if:
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

For a one-part (full-length) insert, only Step 1 changes:

```bash
# Step 1: trim both overhangs off every sequence in the file
python combine.py --single \
    full_insertion/20261104_twist.fasta \
    --output-dir output \
    --output-name full_insertion.fasta

# Step 2: same as above, pointed at the Step 1 output
singularity exec /net/software/containers/universal.sif python insert.py \
    output/full_insertion.fasta \
    backbone/CMVR_insertto1386.gb \
    1386 \
    --output-dir output/inserted
```
