#!/usr/bin/env python3

import argparse
import shutil


def command_status(command):
    return "ok" if shutil.which(command) else "missing"


def main():
    parser = argparse.ArgumentParser(description="Write a TSV report of command availability.")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--include-bbmap", action="store_true", help="Include bbmap.sh in required command list")
    args = parser.parse_args()

    commands = ["python3", "trimmomatic", "fastqc", "sortmerna"]
    if args.include_bbmap:
        commands.append("bbmap.sh")

    with open(args.output, "w", newline="") as handle:
        handle.write("command\tstatus\tpath\n")
        for command in commands:
            path = shutil.which(command) or ""
            handle.write(f"{command}\t{command_status(command)}\t{path}\n")

    print(f"Wrote dependency report for {len(commands)} command(s) to {args.output}")


if __name__ == "__main__":
    main()
