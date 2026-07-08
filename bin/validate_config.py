#!/usr/bin/env python3

import argparse
import csv
import os
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


def resolve_sample_path(samplesheet, value, base_dir):
    if os.path.isabs(value):
        return value
    launch_relative = os.path.abspath(os.path.join(base_dir, value))
    sheet_relative = os.path.abspath(os.path.join(os.path.dirname(samplesheet), value))
    return launch_relative if os.path.exists(launch_relative) else sheet_relative


def main():
    parser = argparse.ArgumentParser(description="Validate metatranscriptomics sample TSV.")
    parser.add_argument("--samples", required=True, help="Tab-delimited sample sheet")
    parser.add_argument("--base-dir", default=os.getcwd(), help="Directory for resolving relative read paths")
    args = parser.parse_args()

    samplesheet = os.path.abspath(args.samples)
    sample_ids = set()
    count = 0

    for row in rows_from_tsv(samplesheet):
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
            resolved = resolve_sample_path(samplesheet, read_path, args.base_dir)
            if not os.path.exists(resolved):
                raise FileNotFoundError(f"Sample {sample_id} {read_col} not found: {read_path}")

    if count == 0:
        raise ValueError("Samplesheet contains no samples")

    print(f"Validated {count} sample(s) from {args.samples}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
