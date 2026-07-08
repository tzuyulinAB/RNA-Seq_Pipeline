#!/usr/bin/env python3

import argparse
import csv
import sys


REQUIRED_COLUMNS = ("sample_id", "read1", "read2")


def rows_from_tsv(path):
    with open(path, newline="") as handle:
        lines = (line for line in handle if not line.lstrip().startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Samplesheet is empty")

        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"Samplesheet is missing required column(s): {', '.join(missing)}")

        for row in reader:
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            if any(cleaned.values()):
                yield cleaned


def main():
    parser = argparse.ArgumentParser(description="Validate metatranscriptomics sample TSV.")
    parser.add_argument("--samples", required=True, help="Tab-delimited sample sheet")
    parser.add_argument("--base-dir", help="Deprecated; retained for compatibility")
    args = parser.parse_args()

    sample_ids = set()
    count = 0

    for row in rows_from_tsv(args.samples):
        count += 1
        sample_id = row["sample_id"]
        if not sample_id:
            raise ValueError(f"Row {count} has an empty sample_id")
        if sample_id in sample_ids:
            raise ValueError(f"Duplicate sample_id: {sample_id}")
        sample_ids.add(sample_id)

        for read_col in ("read1", "read2"):
            read_path = row[read_col]
            if not read_path:
                raise ValueError(f"Sample {sample_id} has an empty {read_col}")

    if count == 0:
        raise ValueError("Samplesheet contains no samples")

    print(f"Validated {count} sample(s) from {args.samples}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
