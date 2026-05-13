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


def combine(antigen_fasta, cage_fasta, combinations, output_dir):
    a_five, a_three, antigen_map = load_fasta(antigen_fasta, "A")
    c_five, c_three, cage_map = load_fasta(cage_fasta, "C")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "combined.fasta"

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

    with open(output_file, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n\n")

    print(f"Wrote {len(records)} sequence(s) to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Combine antigen and cage DNA sequences."
    )
    parser.add_argument("antigen_fasta", help="Antigen FASTA file")
    parser.add_argument("cage_fasta", help="Cage FASTA file")
    parser.add_argument(
        "combinations",
        nargs="+",
        metavar="AxCy",
        help="One or more combinations, e.g. A1C6 A2C7",
    )
    parser.add_argument("--output-dir", default="output", help="Output directory")
    args = parser.parse_args()

    combine(args.antigen_fasta, args.cage_fasta, args.combinations, args.output_dir)


if __name__ == "__main__":
    main()
