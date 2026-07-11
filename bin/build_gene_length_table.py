#!/usr/bin/env python3

import argparse
import gzip
import sys


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def clean_gene_id(header):
    return header[1:].split("#", 1)[0].strip()


def emit_record(out, gene_id, length):
    if gene_id:
        out.write(f"{gene_id}\t{length}\n")


def main():
    parser = argparse.ArgumentParser(description="Build a two-column gene length table from FASTA files.")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("fastas", nargs="+", help="Reference FASTA files")
    args = parser.parse_args()

    with open(args.output, "w") as out:
        out.write("gene_id\tlength\n")
        for fasta in args.fastas:
            gene_id = None
            length = 0
            with open_text(fasta) as handle:
                for line in handle:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        emit_record(out, gene_id, length)
                        gene_id = clean_gene_id(line)
                        length = 0
                    else:
                        length += len(line.strip())
                emit_record(out, gene_id, length)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
