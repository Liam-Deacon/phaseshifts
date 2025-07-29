#!/usr/bin/env python3

"""Helper script to split the Barbieri/Van Hove phshift2007 package AB files into individual program files."""

import sys
import os
import re


def split_ab(input_path, output_dir):
    # type: (str, str) -> None
    with open(input_path, "r") as f:
        lines = f.readlines()
    prog_re = re.compile(r"^C\s+program\s+(\S+)")
    current_file = None
    out_f = None
    for line in lines:
        m = prog_re.match(line)
        if m:
            if out_f:
                out_f.close()
            fname = m.group(1).lower()
            out_path = os.path.join(output_dir, fname)
            out_f = open(out_path, "w")
        if out_f:
            out_f.write(line)
    if out_f:
        out_f.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: split_ab.py <input_ab_file> <output_dir>")
        sys.exit(1)
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    split_ab(input_path, output_dir)
