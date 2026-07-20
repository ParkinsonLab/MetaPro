#!/usr/bin/env python
# output_read_accounting.py
#
# Cumulative read-accounting logger for MetaPro's filtering pipeline.
# Appends one row per (role, label, file) entry to a shared TSV log, plus a
# per-invocation summary block showing input/output totals and the
# unaccounted delta between them (input_total - output_total).
#
# This script places NO marker and does no resume-skip logic of its own --
# by design (per your instruction) it re-runs and re-appends every time the
# stage that calls it executes, including on resumed pipeline runs. That
# means the log can accumulate duplicate rows for a stage across multiple
# runs of the pipeline over the same output directory.
#
# "qf_dedup" and "repop" are expected to balance to (near) zero
# unaccounted, same as every other stage: qf_dedup's accounting entries
# include the .clstr cluster file as an explicit output (the duplicates
# it removes from the deduped FASTQ, not zero the caller loses), and
# repop's accounting entries include the same .clstr file as an explicit
# input (the duplicates it reintroduces). A small positive "unaccounted"
# can still appear for "repop" specifically: if host/vector filtering
# removed a cluster's representative read before repop runs, that
# representative's duplicates are never reintroduced (repop only looks
# up cluster members for representatives present in the post-filtering
# mRNA/other files) -- those reads are still counted as a repop "input"
# here (since they're still sitting in the .clstr file) but never make
# it to a repop "output". That is a real, expected consequence of
# host/vector filtering running before repop, not an accounting bug.
# Everywhere else, a positive "unaccounted" value means some reads could
# not be traced from input file(s) to output file(s).
#
# Usage:
#   python output_read_accounting.py <stage_name> <log_path> <role> <label> <filepath> [<role> <label> <filepath> ...]
#     role: "input" or "output"

import sys
import os
from datetime import datetime as dt


def fastq_count(item):
    """
    Counts reads in a FASTQ file (line_count / 4).
    Returns 0 if the file doesn't exist.
    """
    if not os.path.exists(item):
        return 0
    try:
        with open(item, "rb") as f:
            lines = sum(1 for _ in f)
    except Exception as e:
        print(f"Error counting {item}: {e}")
        return 0

    if lines % 4 != 0:
        print(dt.today(), "Warning: non-standard FASTQ line count in:", item, lines / 4)

    return lines // 4


def clstr_duplicate_count(item):
    """
    Counts "duplicate" reads recorded in a cd-hit-dup .clstr file: reads
    that QF's dedup step collapsed into a representative, so they are NOT
    present in the deduped output FASTQ, but ARE recoverable from this
    cluster file -- which is exactly what read_repopulation.py's
    repopulate_single() reads to reintroduce them.

    This mirrors that function's own cluster_map parsing convention
    line-for-line (skip blank/">"-header lines; a member line starting
    with "0" is the cluster representative, already counted in the
    deduped output FASTQ; any other member line is a duplicate) so the
    two counts cannot drift apart.

    Returns 0 if the file doesn't exist.
    """
    if not os.path.exists(item):
        return 0
    duplicate_count = 0
    try:
        with open(item, "r") as clstr_in:
            for line in clstr_in:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                if line.startswith("0"):
                    # cluster representative -- already in the deduped
                    # output FASTQ, not a duplicate.
                    continue
                duplicate_count += 1
    except Exception as e:
        print(f"Error counting {item}: {e}")
        return 0

    return duplicate_count


def main():
    # 2 fixed args (stage_name, log_path) + N triplets of (role, label, filepath)
    if len(sys.argv) < 6 or (len(sys.argv) - 3) % 3 != 0:
        print(f"Usage: python {sys.argv[0]} <stage_name> <log_path> <role> <label> <filepath> [<role> <label> <filepath> ...]")
        sys.exit(1)

    stage_name = sys.argv[1]
    log_path = sys.argv[2]
    raw_entries = sys.argv[3:]

    entries = []
    for i in range(0, len(raw_entries), 3):
        role = raw_entries[i]
        label = raw_entries[i + 1]
        filepath = raw_entries[i + 2]
        if role not in ("input", "output"):
            print(f"Error: unrecognized role '{role}' for label '{label}'. Must be 'input' or 'output'.")
            sys.exit(1)
        entries.append((role, label, filepath))

    timestamp = str(dt.today())
    log_is_new = not os.path.exists(log_path)

    input_total = 0
    output_total = 0
    rows = []
    for role, label, filepath in entries:
        if filepath.endswith(".clstr"):
            count = clstr_duplicate_count(filepath)
        else:
            count = fastq_count(filepath)
        rows.append((stage_name, timestamp, role, label, str(int(count)), filepath))
        if role == "input":
            input_total += count
        else:
            output_total += count

    unaccounted = input_total - output_total

    with open(log_path, "a") as out:
        if log_is_new:
            out.write("\t".join(["stage", "timestamp", "role", "label", "read_count", "source_file"]) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
        out.write("\t".join([stage_name, timestamp, "summary", "input_total", str(int(input_total)), ""]) + "\n")
        out.write("\t".join([stage_name, timestamp, "summary", "output_total", str(int(output_total)), ""]) + "\n")
        out.write("\t".join([stage_name, timestamp, "summary", "unaccounted", str(int(unaccounted)), ""]) + "\n")


if __name__ == "__main__":
    main()