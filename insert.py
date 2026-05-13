#!/usr/bin/env python3
"""
Insert sequences from a FASTA file into a GenBank backbone at a given site.

Usage:
    python insert.py <insertion.fasta> <backbone.gb> <insertion_site> [--output-dir DIR]

insertion_site: 1-based position; insertion starts AT this position.
"""

import argparse
import re
import sys
from copy import deepcopy
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


# ---------------------------------------------------------------------------
# Sequence validation
# ---------------------------------------------------------------------------

def validate(dna: str, label: str) -> None:
    dna = dna.upper()
    issues = []

    if not dna.startswith("ATG"):
        issues.append(f"does not start with ATG (starts with '{dna[:3]}')")

    protein = str(Seq(dna).translate())

    if not protein.endswith("*"):
        issues.append("does not end with a stop codon")
    else:
        premature = protein[:-1].count("*")
        if premature:
            issues.append(f"contains {premature} premature stop codon(s)")

    for msg in issues:
        print(f"WARNING [{label}]: insertion sequence {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# GenBank post-processing
# ---------------------------------------------------------------------------

def fix_locus_name(gb_path: Path, full_name: str) -> None:
    """Replace the LOCUS name (truncated by BioPython) with the full label."""
    text = gb_path.read_text()
    text = re.sub(
        r"^(LOCUS\s+)\S+",
        lambda m: m.group(1) + full_name,
        text,
        count=1,
        flags=re.MULTILINE,
    )
    gb_path.write_text(text)


# ---------------------------------------------------------------------------
# Location shifting
# ---------------------------------------------------------------------------

def _shift_interval(start: int, end: int, split: int, ins_len: int):
    """Shift a half-open [start, end) interval given a 0-based split point."""
    if start >= split:
        return start + ins_len, end + ins_len
    if end > split:
        return start, end + ins_len   # feature spans insertion → extend
    return start, end


def shift_location(loc, split: int, ins_len: int):
    """Recursively shift a FeatureLocation or CompoundLocation."""
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [shift_location(p, split, ins_len) for p in loc.parts],
            loc.operator,
        )
    s, e = _shift_interval(int(loc.start), int(loc.end), split, ins_len)
    return FeatureLocation(s, e, strand=loc.strand)


# ---------------------------------------------------------------------------
# Backbone insertion
# ---------------------------------------------------------------------------

def insert_into_backbone(
    backbone: SeqRecord, ins_seq: str, site_1based: int, label: str
) -> SeqRecord:
    """
    Insert ins_seq into backbone after 1-based position site_1based.
    Returns a new SeqRecord with updated features.
    """
    split = site_1based - 1      # 0-based split index: insert BEFORE 1-based position N
    ins_len = len(ins_seq)
    original = str(backbone.seq)

    if split < 0 or split > len(original):
        raise ValueError(
            f"Insertion site {site_1based} is out of range for backbone of length {len(original)}"
        )

    new_seq_str = original[:split] + ins_seq + original[split:]

    # Shift existing features
    new_features = []
    for feat in backbone.features:
        new_loc = shift_location(feat.location, split, ins_len)
        new_feat = SeqFeature(
            location=new_loc,
            type=feat.type,
            qualifiers=deepcopy(feat.qualifiers),
            id=feat.id,
        )
        new_features.append(new_feat)

    # Annotation feature for the inserted CDS
    ins_feature = SeqFeature(
        location=FeatureLocation(split, split + ins_len, strand=1),
        type="CDS",
        qualifiers={
            "label": [label],
            "note": ["inserted sequence"],
        },
    )
    new_features.append(ins_feature)

    new_record = SeqRecord(
        Seq(new_seq_str),
        id=label,
        name=label[:16],           # LOCUS name field, ≤ 16 chars
        description=label,
        features=new_features,
        annotations=deepcopy(backbone.annotations),
    )
    new_record.annotations["molecule_type"] = backbone.annotations.get(
        "molecule_type", "DNA"
    )
    new_record.annotations["topology"] = "circular"
    return new_record


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("insertion_fasta", help="FASTA file with sequences to insert")
    parser.add_argument("backbone_gb", help="GenBank backbone file")
    parser.add_argument(
        "insertion_site",
        type=int,
        help="1-based position; insertion happens AFTER this position",
    )
    parser.add_argument("--output-dir", default="output", help="Output directory (default: output)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    backbone = SeqIO.read(args.backbone_gb, "genbank")
    insertions = list(SeqIO.parse(args.insertion_fasta, "fasta"))

    if not insertions:
        print("ERROR: no sequences found in insertion FASTA.", file=sys.stderr)
        sys.exit(1)

    for rec in insertions:
        label = rec.id
        ins_seq = str(rec.seq).upper()

        validate(ins_seq, label)

        new_record = insert_into_backbone(backbone, ins_seq, args.insertion_site, label)

        gb_out = output_dir / f"{label}.gb"
        SeqIO.write(new_record, gb_out, "genbank")
        fix_locus_name(gb_out, label)
        print(f"Written {gb_out}  ({len(new_record)} bp)")


if __name__ == "__main__":
    main()
