#!/usr/bin/env python3

"""
Add INFO fields to VCF Header
"""

import argparse
import os
from typing import List, Tuple


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Add INFO fields to a VCF header.")
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to the input VCF file.",
    )
    parser.add_argument(
        "--info",
        action="append",
        nargs=4,
        metavar=("ID", "Number", "Type", "Description"),
        required=True,
        help=(
            "Define an INFO field in the header. Provide ID, Number, Type, and Description. "
            "This argument can be repeated to add multiple INFO fields."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output VCF file.",
    )
    return parser.parse_args()


def generate_info_lines(info_fields: List[List[str]]) -> List[str]:
    info_lines = []
    for info in info_fields:
        id_, number, type_, description = info
        info_line = f'##INFO=<ID={id_},Number={number},Type={type_},Description="{description}">'
        info_lines.append(info_line)
    return info_lines


def add_info_to_vcf(input_vcf: str, output_vcf: str, info_lines: List[str]) -> None:
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#CHROM"):
                # Insert INFO lines just before the column header line
                for info_line in info_lines:
                    outfile.write(info_line + "\n")
            outfile.write(line)


def main() -> None:
    args = parse_arguments()

    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"Input VCF file '{args.vcf}' does not exist.")

    info_lines_to_be_added = generate_info_lines(args.info)
    add_info_to_vcf(args.vcf, args.output, info_lines_to_be_added)


if __name__ == "__main__":
    main()
