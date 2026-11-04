import argparse
import re
from pathlib import Path


def parse_fasta(filepath):
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        sequences[current_header] = "".join(current_seq)
    return sequences


OVERHANG_RE = re.compile(r"^(\w+)_([53])overhang$")


def find_overhangs(seqs, fasta_path):
    """Return (5overhang, 3overhang) from the overhang entries of a FASTA.

    The two prefixes need not match: a full-length insert spans the antigen's
    5' end and the cage's 3' end, so it carries A_5overhang and C_3overhang.
    """
    found = {}
    for header, seq in seqs.items():
        m = OVERHANG_RE.match(header)
        if not m:
            continue
        end = m.group(2)
        if end in found:
            raise ValueError(
                f"More than one {end}' overhang entry found in {fasta_path}"
            )
        found[end] = seq

    for end in ("5", "3"):
        if end not in found:
            raise ValueError(
                f"No '<prefix>_{end}overhang' entry found in {fasta_path}"
            )
    return found["5"], found["3"]


def resolve_selection(selection, parts, fasta_path):
    """Resolve one user selection to a (header, seq) pair.

    A selection may be a full entry name ('38g_floppy_A3C3'), a prefixed ID
    ('A3' -> 'A3_38g-Yoan'), or a bare number ('3' -> 'C3_...').
    """
    if selection in parts:
        return selection, parts[selection]

    m = re.match(r"^([A-Za-z]*)(\d+)$", selection)
    if m:
        prefix, number = m.group(1), m.group(2)
        pattern = rf"^{re.escape(prefix)}{number}_" if prefix else rf"^[A-Za-z]*{number}_"
        matches = [h for h in parts if re.match(pattern, h)]
        if len(matches) == 1:
            return matches[0], parts[matches[0]]
        if len(matches) > 1:
            raise ValueError(
                f"Selection '{selection}' is ambiguous in {fasta_path}; "
                f"it matches: {', '.join(matches)}"
            )

    raise ValueError(
        f"Selection '{selection}' not found in {fasta_path}.\n"
        f"  Available: {', '.join(parts)}"
    )


def load_fasta(fasta_path, prefix):
    """Return (5overhang, 3overhang, {id_int: (name, seq)})."""
    seqs = parse_fasta(fasta_path)

    five_oh = seqs.get(f"{prefix}_5overhang")
    three_oh = seqs.get(f"{prefix}_3overhang")

    if five_oh is None:
        raise ValueError(f"No '{prefix}_5overhang' entry found in {fasta_path}")
    if three_oh is None:
        raise ValueError(f"No '{prefix}_3overhang' entry found in {fasta_path}")

    id_map = {}
    for header, seq in seqs.items():
        m = re.match(rf"^{re.escape(prefix)}(\d+)_(.+)$", header)
        if m:
            id_map[int(m.group(1))] = (m.group(2), seq)

    return five_oh, three_oh, id_map


def strip_seq(seq, five_oh, three_oh, label, strip_five, strip_three):
    if strip_five:
        if not seq.upper().startswith(five_oh.upper()):
            raise ValueError(
                f"{label}: 5' overhang not found at start.\n"
                f"  Expected : {five_oh}\n"
                f"  Got      : {seq[:len(five_oh)]}"
            )
        seq = seq[len(five_oh):]

    if strip_three:
        if not seq.upper().endswith(three_oh.upper()):
            raise ValueError(
                f"{label}: 3' overhang not found at end.\n"
                f"  Expected : {three_oh}\n"
                f"  Got      : {seq[-len(three_oh):]}"
            )
        seq = seq[: len(seq) - len(three_oh)]

    return seq


def write_records(records, output_dir, output_name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / output_name

    with open(output_file, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n\n")

    print(f"Wrote {len(records)} sequence(s) to {output_file}")


def combine_single(part_fasta, selections, output_dir, output_name):
    """Single-part mode: trim both overhangs off each selected sequence.

    The result has the same shape as the two-part output, so insert.py runs
    on it unchanged.
    """
    seqs = parse_fasta(part_fasta)
    five_oh, three_oh = find_overhangs(seqs, part_fasta)

    parts = {h: s for h, s in seqs.items() if not OVERHANG_RE.match(h)}
    if not parts:
        raise ValueError(f"No sequence entries found in {part_fasta}")

    if selections:
        chosen = [resolve_selection(sel, parts, part_fasta) for sel in selections]
    else:
        chosen = list(parts.items())

    records = []
    for header, seq in chosen:
        core = strip_seq(
            seq, five_oh, three_oh, header, strip_five=True, strip_three=True
        )
        records.append((header, core))

    write_records(records, output_dir, output_name)


def combine(antigen_fasta, cage_fasta, combinations, output_dir, output_name):
    a_five, a_three, antigen_map = load_fasta(antigen_fasta, "A")
    c_five, c_three, cage_map = load_fasta(cage_fasta, "C")

    records = []
    for combo in combinations:
        m = re.match(r"^A(\d+)C(\d+)$", combo)
        if not m:
            raise ValueError(
                f"Invalid combination '{combo}'. Expected format: A<int>C<int>"
            )
        a_id, c_id = int(m.group(1)), int(m.group(2))

        if a_id not in antigen_map:
            raise ValueError(f"Antigen A{a_id} not found in {antigen_fasta}")
        if c_id not in cage_map:
            raise ValueError(f"Cage C{c_id} not found in {cage_fasta}")

        a_name, a_seq = antigen_map[a_id]
        c_name, c_seq = cage_map[c_id]

        # Antigen: remove both overhangs
        a_core = strip_seq(a_seq, a_five, a_three, f"A{a_id}", strip_five=True, strip_three=True)
        # Cage: remove only the 3' overhang, keep the 5' overhang
        c_no3 = strip_seq(c_seq, c_five, c_three, f"C{c_id}", strip_five=False, strip_three=True)

        header = f"A{a_id}C{c_id}_{a_name}_{c_name}"
        records.append((header, a_core + c_no3))

    write_records(records, output_dir, output_name)


def main():
    parser = argparse.ArgumentParser(
        description="Combine antigen and cage DNA sequences, or prepare a "
        "single-part insert.",
        epilog="Two-part mode:   combine.py antigen.fasta cage.fasta A1C6 A2C7\n"
        "Single-part mode: combine.py --single full_insertion/20261104_twist.fasta",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Single-part mode: one FASTA file, optionally followed by "
        "sequence IDs or names (default: every sequence in the file). "
        "Both overhangs are trimmed and no joining is done.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="ARG",
        help="Two-part mode: antigen.fasta cage.fasta AxCy [AxCy ...]. "
        "Single-part mode: part.fasta [ID ...]",
    )
    parser.add_argument("--output-dir", default="output", help="Output directory")
    parser.add_argument(
        "--output-name",
        default="combined.fasta",
        help="Output FASTA file name (default: combined.fasta)",
    )
    args = parser.parse_args()

    if args.single:
        combine_single(
            args.inputs[0], args.inputs[1:], args.output_dir, args.output_name
        )
    else:
        if len(args.inputs) < 3:
            parser.error(
                "two-part mode needs antigen.fasta, cage.fasta and at least one "
                "AxCy combination (use --single for a one-part insert)"
            )
        combine(
            args.inputs[0],
            args.inputs[1],
            args.inputs[2:],
            args.output_dir,
            args.output_name,
        )


if __name__ == "__main__":
    main()
