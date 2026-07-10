#!/usr/bin/env python3

import argparse
import gzip
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request


DATABASE_URLS = (
    "https://github.com/sortmerna/sortmerna/releases/download/v4.3.6/database.tar.gz",
    "https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz",
)

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def copy_fasta(path, output):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as in_handle, open(output, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)


def output_name_for_archive_member(filename):
    if filename.endswith(".gz"):
        filename = filename[:-3]
    return filename


def extract_all_databases(archive_path, outdir):
    with tempfile.TemporaryDirectory() as tmpdir:
        with tarfile.open(archive_path, "r:gz") as archive:
            archive.extractall(tmpdir)

        fastas = []
        for root, _, files in os.walk(tmpdir):
            for filename in files:
                if filename.endswith(FASTA_EXTENSIONS):
                    fastas.append(os.path.join(root, filename))
        if not fastas:
            raise RuntimeError("Downloaded SortMeRNA archive did not contain FASTA files")

        extracted = []
        for fasta in sorted(fastas):
            out_name = output_name_for_archive_member(os.path.basename(fasta))
            out_path = os.path.join(outdir, out_name)
            copy_fasta(fasta, out_path)
            extracted.append(out_path)
        return extracted


def main():
    parser = argparse.ArgumentParser(description="Extract SortMeRNA database FASTAs.")
    parser.add_argument("--outdir", required=True, help="Directory for database output")
    parser.add_argument("--output", required=True, help="Selected reference FASTA copy for the workflow")
    parser.add_argument("--selected", required=True, help="Selected database FASTA basename to use")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    output = args.output if os.path.isabs(args.output) else os.path.join(args.outdir, args.output)
    selected = os.path.basename(args.selected)

    errors = []
    for url in DATABASE_URLS:
        with tempfile.NamedTemporaryFile(suffix=".tar.gz") as tmp:
            print(f"Downloading SortMeRNA database from {url}")
            try:
                urllib.request.urlretrieve(url, tmp.name)
                extracted = extract_all_databases(tmp.name, args.outdir)
                by_name = {os.path.basename(path): path for path in extracted}
                print("Extracted SortMeRNA database FASTAs:")
                for name in sorted(by_name):
                    print(f"  {name}")

                if selected not in by_name:
                    raise RuntimeError(
                        f"Selected database {selected!r} was not found in the archive. "
                        f"Available: {', '.join(sorted(by_name))}"
                    )

                shutil.copyfile(by_name[selected], output)
                print(f"Selected {selected} for this workflow run as {output}")
                return
            except Exception as exc:
                errors.append(f"{url}: {exc}")

    raise RuntimeError(
        "Could not extract SortMeRNA databases from release downloads. "
        "Provide an existing reference with --sortmerna_ref. Attempts: " + " | ".join(errors)
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
