#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys
import urllib.request


ADAPTER_BASE_URL = "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters"


def candidate_paths(selected):
    candidates = []
    trimmomatic = shutil.which("trimmomatic")
    if trimmomatic:
        try:
            output = subprocess.check_output([trimmomatic], stderr=subprocess.STDOUT, text=True)
        except Exception:
            output = ""
        for token in output.replace(":", " ").split():
            if token.endswith(".jar") and os.path.exists(token):
                candidates.append(os.path.join(os.path.dirname(token), "adapters", selected))

    common_dirs = [
        "/usr/share/trimmomatic",
        "/usr/share/trimmomatic/adapters",
        "/usr/local/share/trimmomatic",
        "/opt/conda/share/trimmomatic",
        "/opt/conda/share/trimmomatic/adapters",
    ]
    candidates.extend(os.path.join(directory, selected) for directory in common_dirs)
    return candidates


def main():
    parser = argparse.ArgumentParser(description="Prepare a Trimmomatic adapter FASTA.")
    parser.add_argument("--adapter-dir", required=True, help="Directory for adapter FASTA output")
    parser.add_argument("--selected", required=True, help="Adapter FASTA name, e.g. NexteraPE-PE.fa")
    args = parser.parse_args()

    os.makedirs(args.adapter_dir, exist_ok=True)
    output = os.path.join(args.adapter_dir, args.selected)

    if os.path.exists(output):
        print(f"Adapter already exists: {output}")
        return

    for candidate in candidate_paths(args.selected):
        if os.path.exists(candidate):
            shutil.copyfile(candidate, output)
            print(f"Copied adapter from {candidate} to {output}")
            return

    url = f"{ADAPTER_BASE_URL}/{args.selected}"
    print(f"Downloading adapter from {url}")
    try:
        urllib.request.urlretrieve(url, output)
    except Exception as exc:
        raise RuntimeError(
            "Could not find adapter in local Trimmomatic paths and download failed. "
            f"Provide --adapter_fasta or check network access. Details: {exc}"
        )

    print(f"Wrote adapter to {output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
