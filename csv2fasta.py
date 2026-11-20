#!/usr/bin/env python3
"""
Convert a downloaded order CSV (Twist, IDT, ...) into a FASTA file that
combine.py can read.

Each CSV row becomes one FASTA record: the 'Name' column is the header and the
'Insert Sequence' column is the sequence. The overhangs are not taken from the
CSV -- they are supplied separately, either as a FASTA file holding
'<prefix>_5overhang' and '<prefix>_3overhang' entries or inline with
--five/--three, and are written at the top of the output file so the result has
the layout the rest of the pipeline expects.

Usage:
    python csv2fasta.py <order.csv> <overhangs.fasta> [--output-dir DIR]
    python csv2fasta.py <order.csv> --five GATCA... --three CCTAG...
"""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from combine import OVERHANG_RE, find_overhangs, parse_fasta  # noqa: E402

DEFAULT_NAME_COLUMN = "Name"
DEFAULT_SEQUENCE_COLUMN = "Insert Sequence"

# Order CSVs are inconsistent about capitalisation ('Insert sequence' vs
# 'Insert Sequence'), so columns are matched case- and space-insensitively.
def _normalize(field):
    return re.sub(r"[\s_]+", "", field).lower()


def find_column(fieldnames, wanted, csv_path):
    target = _normalize(wanted)
    matches = [f for f in fieldnames if _normalize(f) == target]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Column '{wanted}' is ambiguous in {csv_path}; "
            f"it matches: {', '.join(matches)}"
        )
    raise ValueError(
        f"No '{wanted}' column found in {csv_path}.\n"
        f"  Available: {', '.join(fieldnames)}"
    )


def clean_sequence(raw):
    """Drop whitespace (CSV cells may wrap) and upper-case."""
    return re.sub(r"\s+", "", raw).upper()


def read_csv_records(csv_path, name_column, sequence_column):
    """Return [(name, seq)] in file order."""
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"{csv_path} is empty")

        name_field = find_column(reader.fieldnames, name_column, csv_path)
        seq_field = find_column(reader.fieldnames, sequence_column, csv_path)

        records = []
        seen = {}
        for lineno, row in enumerate(reader, start=2):
            name = (row.get(name_field) or "").strip()
            seq = clean_sequence(row.get(seq_field) or "")

            if not name and not seq:
                continue  # trailing blank line
            if not name:
                raise ValueError(f"{csv_path} line {lineno}: empty '{name_field}'")
            if not seq:
                raise ValueError(
                    f"{csv_path} line {lineno}: empty '{seq_field}' for '{name}'"
                )
            if name in seen:
                raise ValueError(
                    f"{csv_path} line {lineno}: duplicate name '{name}' "
                    f"(first seen on line {seen[name]})"
                )

            seen[name] = lineno
            records.append((name, seq))

    if not records:
        raise ValueError(f"No data rows found in {csv_path}")
    return records


def check_record(name, seq, five_oh, three_oh):
    """Warn about anything that would break the downstream steps."""
    bad = sorted(set(re.sub(r"[ACGT]", "", seq)))
    if bad:
        print(
            f"WARNING [{name}]: non-ACGT character(s) in sequence: "
            f"{', '.join(bad)}",
            file=sys.stderr,
        )

    if not seq.startswith(five_oh.upper()):
        print(
            f"WARNING [{name}]: does not start with the 5' overhang\n"
            f"  Expected : {five_oh}\n"
            f"  Got      : {seq[:len(five_oh)]}",
            file=sys.stderr,
        )
    if not seq.endswith(three_oh.upper()):
        print(
            f"WARNING [{name}]: does not end with the 3' overhang\n"
            f"  Expected : {three_oh}\n"
            f"  Got      : {seq[-len(three_oh):]}",
            file=sys.stderr,
        )


def load_overhangs(overhang_fasta, five, three):
    """Return (five_header, five_seq, three_header, three_seq)."""
    if overhang_fasta:
        seqs = parse_fasta(overhang_fasta)
        # find_overhangs validates that exactly one entry of each end is present
        find_overhangs(seqs, overhang_fasta)
        headers = {}
        for header, seq in seqs.items():
            m = OVERHANG_RE.match(header)
            if m:
                headers[m.group(2)] = (header, clean_sequence(seq))
        (five_header, five_seq) = headers["5"]
        (three_header, three_seq) = headers["3"]
        return five_header, five_seq, three_header, three_seq

    return "A_5overhang", clean_sequence(five), "C_3overhang", clean_sequence(three)


def write_fasta(entries, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n{seq}\n\n")
    print(f"Wrote {len(entries)} sequence(s) to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("csv_file", help="Downloaded order CSV")
    parser.add_argument(
        "overhang_fasta",
        nargs="?",
        help="FASTA file with '<prefix>_5overhang' and '<prefix>_3overhang' "
        "entries (alternative to --five/--three)",
    )
    parser.add_argument("--five", help="5' overhang sequence, given inline")
    parser.add_argument("--three", help="3' overhang sequence, given inline")
    parser.add_argument(
        "--name-column",
        default=DEFAULT_NAME_COLUMN,
        help=f"CSV column holding the record name (default: {DEFAULT_NAME_COLUMN})",
    )
    parser.add_argument(
        "--sequence-column",
        default=DEFAULT_SEQUENCE_COLUMN,
        help="CSV column holding the sequence "
        f"(default: {DEFAULT_SEQUENCE_COLUMN})",
    )
    parser.add_argument(
        "--output-dir",
        default="full_insertion",
        help="Output directory (default: full_insertion)",
    )
    parser.add_argument(
        "--output-name",
        help="Output FASTA file name (default: the CSV's name with a .fasta "
        "extension)",
    )
    args = parser.parse_args()

    if args.overhang_fasta and (args.five or args.three):
        parser.error("give either an overhang FASTA file or --five/--three, not both")
    if not args.overhang_fasta and not (args.five and args.three):
        parser.error(
            "overhangs are required: pass an overhang FASTA file, or both "
            "--five and --three"
        )

    try:
        five_header, five_seq, three_header, three_seq = load_overhangs(
            args.overhang_fasta, args.five, args.three
        )
        records = read_csv_records(
            args.csv_file, args.name_column, args.sequence_column
        )
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    for name, seq in records:
        check_record(name, seq, five_seq, three_seq)

    output_name = args.output_name or f"{Path(args.csv_file).stem}.fasta"
    entries = [(five_header, five_seq), (three_header, three_seq)] + records
    write_fasta(entries, Path(args.output_dir) / output_name)


if __name__ == "__main__":
    main()
