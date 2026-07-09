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

DIRECT_FASTA_URLS = (
    "https://github.com/sortmerna/sortmerna/releases/download/v4.3.6/{name}",
    "https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/{name}",
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


def find_named_fasta(root_dir, output_name):
    matches = []
    for root, _, files in os.walk(root_dir):
        for filename in files:
            if filename == output_name:
                matches.append(os.path.join(root, filename))
    return sorted(matches)


def extract_requested_database(archive_path, output, output_name):
    with tempfile.TemporaryDirectory() as tmpdir:
        with tarfile.open(archive_path, "r:gz") as archive:
            archive.extractall(tmpdir)

        exact_matches = find_named_fasta(tmpdir, output_name)
        if exact_matches:
            shutil.copyfile(exact_matches[0], output)
            return f"copied {output_name}"

        fastas = []
        for root, _, files in os.walk(tmpdir):
            for filename in files:
                if filename.endswith((".fa", ".fasta", ".fna")):
                    fastas.append(os.path.join(root, filename))
        if not fastas:
            raise RuntimeError("Downloaded SortMeRNA archive did not contain FASTA files")

        concatenate_fastas(sorted(fastas), output)
        return f"concatenated {len(fastas)} FASTA file(s)"


def download_direct_fasta(output_name, output):
    errors = []
    for template in DIRECT_FASTA_URLS:
        url = template.format(name=output_name)
        print(f"Trying direct SortMeRNA FASTA download from {url}")
        try:
            with tempfile.NamedTemporaryFile(delete=False) as tmp:
                tmp_path = tmp.name
            try:
                urllib.request.urlretrieve(url, tmp_path)
                if os.path.getsize(tmp_path) < 100:
                    with open(tmp_path, "rb") as handle:
                        snippet = handle.read(100).decode("utf-8", errors="replace")
                    raise RuntimeError(f"downloaded file is unexpectedly small: {snippet!r}")
                shutil.move(tmp_path, output)
                return url
            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)
        except Exception as exc:
            errors.append(f"{url}: {exc}")
    raise RuntimeError("Direct FASTA download failed. Attempts: " + " | ".join(errors))


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

    output_name = os.path.basename(output)
    direct_errors = None
    if output_name.startswith("smr_"):
        try:
            url = download_direct_fasta(output_name, output)
            print(f"Downloaded {output_name} from {url}")
            return
        except Exception as exc:
            direct_errors = str(exc)
            print(f"Direct FASTA download unavailable for {output_name}: {exc}")

    errors = []
    for url in DATABASE_URLS:
        with tempfile.NamedTemporaryFile(suffix=".tar.gz") as tmp:
            print(f"Downloading SortMeRNA database from {url}")
            try:
                urllib.request.urlretrieve(url, tmp.name)
                message = extract_requested_database(tmp.name, output, output_name)
                print(f"Extracted SortMeRNA archive and {message} into {output}")
                return
            except Exception as exc:
                errors.append(f"{url}: {exc}")

    local_fastas = find_local_fastas()
    if local_fastas:
        concatenate_fastas(local_fastas, output)
        print(f"Concatenated {len(local_fastas)} local SortMeRNA FASTA file(s) into {output}")
        return

    if direct_errors:
        errors.insert(0, direct_errors)

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
