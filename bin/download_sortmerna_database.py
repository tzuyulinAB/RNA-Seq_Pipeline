#!/usr/bin/env python3

import argparse
import glob
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


def find_local_fastas():
    patterns = [
        "/opt/conda/share/sortmerna*/**/*.fasta",
        "/opt/conda/share/sortmerna*/**/*.fa",
        "/usr/local/share/sortmerna*/**/*.fasta",
        "/usr/local/share/sortmerna*/**/*.fa",
        "/usr/share/sortmerna*/**/*.fasta",
        "/usr/share/sortmerna*/**/*.fa",
    ]
    paths = []
    for pattern in patterns:
        paths.extend(glob.glob(pattern, recursive=True))
    return sorted(path for path in set(paths) if os.path.isfile(path))


def concatenate_fastas(paths, output):
    with open(output, "wb") as out_handle:
        for path in paths:
            with open(path, "rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
                out_handle.write(b"\n")


def extract_fastas(archive_path, output):
    with tempfile.TemporaryDirectory() as tmpdir:
        with tarfile.open(archive_path, "r:gz") as archive:
            archive.extractall(tmpdir)
        fastas = []
        for root, _, files in os.walk(tmpdir):
            for filename in files:
                if filename.endswith((".fa", ".fasta", ".fna")):
                    fastas.append(os.path.join(root, filename))
        if not fastas:
            raise RuntimeError("Downloaded SortMeRNA archive did not contain FASTA files")
        concatenate_fastas(sorted(fastas), output)
        return len(fastas)


def main():
    parser = argparse.ArgumentParser(description="Prepare a concatenated SortMeRNA reference FASTA.")
    parser.add_argument("--outdir", required=True, help="Directory for database output")
    parser.add_argument("--output", required=True, help="Output FASTA path")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    output = args.output if os.path.isabs(args.output) else os.path.join(args.outdir, args.output)

    if os.path.exists(output):
        print(f"SortMeRNA reference already exists: {output}")
        return

    local_fastas = find_local_fastas()
    if local_fastas:
        concatenate_fastas(local_fastas, output)
        print(f"Concatenated {len(local_fastas)} local SortMeRNA FASTA file(s) into {output}")
        return

    errors = []
    for url in DATABASE_URLS:
        with tempfile.NamedTemporaryFile(suffix=".tar.gz") as tmp:
            print(f"Downloading SortMeRNA database from {url}")
            try:
                urllib.request.urlretrieve(url, tmp.name)
                count = extract_fastas(tmp.name, output)
                print(f"Extracted and concatenated {count} FASTA file(s) into {output}")
                return
            except Exception as exc:
                errors.append(f"{url}: {exc}")

    raise RuntimeError(
        "Could not prepare SortMeRNA database from local files or release downloads. "
        "Provide an existing reference with --sortmerna_ref. Attempts: " + " | ".join(errors)
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
