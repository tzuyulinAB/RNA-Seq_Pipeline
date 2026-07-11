#!/usr/bin/env python3

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description="Merge per-sample assigned read tables.")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("tables", nargs="+", help="Per-sample assigned read TSV files")
    args = parser.parse_args()

    with open(args.output, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "ORF_id", "assignedReads"])

        for table in args.tables:
            with open(table, newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                for row in reader:
                    writer.writerow([row["Sample"], row["ORF_id"], row["assignedReads"]])


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
