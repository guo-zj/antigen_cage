#!/usr/bin/env python3
"""
Insert sequences from a FASTA file into a GenBank backbone at a given site,
and translate each insertion sequence.

Usage:
    python insert.py <insertion.fasta> <backbone.gb> <insertion_site> [--output-dir DIR]

insertion_site: 1-based position; sequences are inserted AFTER this position.
"""

import argparse
import sys
from copy import deepcopy
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


# ---------------------------------------------------------------------------
# Translation helpers
# ---------------------------------------------------------------------------

def translate_and_validate(dna: str, label: str) -> str:
    """
    Translate dna, check for ATG start and stop codon, warn on premature stops.
    Returns the protein string (with trailing *).
    """
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

    return protein


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
    split = site_1based          # 0-based split index (insert after 1-based N = split at 0-based N)
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
            "translation": [str(Seq(ins_seq.upper()).translate(to_stop=True))],
        },
    )
    new_features.append(ins_feature)

    new_record = SeqRecord(
        Seq(new_seq_str),
        id=backbone.id,
        name=backbone.name[:16],   # GenBank name ≤ 16 chars
        description=f"{backbone.description} | {label}",
        features=new_features,
        annotations=deepcopy(backbone.annotations),
    )
    new_record.annotations["molecule_type"] = backbone.annotations.get(
        "molecule_type", "DNA"
    )
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

    protein_records = []

    for rec in insertions:
        label = rec.id
        ins_seq = str(rec.seq).upper()

        # Validate and translate
        protein = translate_and_validate(ins_seq, label)

        # Build new GenBank record
        new_record = insert_into_backbone(backbone, ins_seq, args.insertion_site, label)

        # Write GenBank file
        gb_out = output_dir / f"{label}.gb"
        SeqIO.write(new_record, gb_out, "genbank")
        print(f"Written {gb_out}  ({len(new_record)} bp)")

        # Collect protein
        protein_records.append(
            SeqRecord(
                Seq(protein),
                id=label,
                description="translated insertion",
            )
        )

    # Write all protein translations to a single FASTA
    prot_out = output_dir / "proteins.faa"
    SeqIO.write(protein_records, prot_out, "fasta")
    print(f"Written protein translations to {prot_out}")


if __name__ == "__main__":
    main()
