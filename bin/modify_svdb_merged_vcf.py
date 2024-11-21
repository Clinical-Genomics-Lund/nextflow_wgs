#!/usr/bin/env python3

"""
This script ensured that the SV `svdb_origin` and `set` INFO values can be
parsed by scout.

It reduces `svdb_origin=manta1|manta2|tiddit1` -> `svdb_origin=manta|tiddit`
and does the same for the `set` field
"""

import sys
import gzip
from collections import OrderedDict
from typing import Optional

SVDB_ORIGIN_SEPARATOR = "|"
SVDB_ORIGIN_KEY = "svdb_origin"
SVDB_SET_KEY = "set"
SVDB_SET_SEPARATOR = "-"


class Callers:
    MANTA = "manta"
    GATK = "gatk"
    TIDDIT = "tiddit"

    @staticmethod
    def match_svdb_annotation_with_caller(
        unparsed_svdb_origin_annotation: str,
    ) -> Optional[str]:
        """
        Convert caller annotations manta1|gatk1|gatk2 -> manta|gatk

        Returns:
         - caller_name (str) i
         - None if caller not supported
        """
        callers = [Callers.MANTA, Callers.GATK, Callers.TIDDIT]

        for caller in callers:
            if caller in unparsed_svdb_origin_annotation:
                return caller

        return None


def main() -> None:
    if len(sys.argv) != 2:
        print(
            "Usage: modify_svdb_merged_vcf.py file.vcf(.gz) > output.vcf",
            file=sys.stderr,
        )
        sys.exit(1)

    infile = sys.argv[1]

    # Open the VCF file (supports gzip)
    if infile.endswith(".gz"):
        f = gzip.open(infile, "rt")
    else:
        f = open(infile, "r")

    for line in f:
        # Retain header lines
        if line.startswith("#"):
            print(line, end="")
            continue

        fields = line.rstrip("\n").split("\t")
        chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
        others = fields[8:]  # If there are additional columns

        info_entries = info.split(";")
        info_dict = OrderedDict()
        for entry in info_entries:
            if "=" in entry:
                key, value = entry.split("=", 1)
                info_dict[key] = value
            else:
                info_dict[entry] = True  # Flag without a value

        # Modify 'set' and 'svdb_origin' in place
        if SVDB_SET_KEY in info_dict:
            info_dict[SVDB_SET_KEY] = reduce_to_set_of_unique_callers(
                info_dict[SVDB_SET_KEY], separator=SVDB_SET_SEPARATOR
            )

        if SVDB_ORIGIN_KEY in info_dict:
            # Modify 'svdb_origin' value as needed
            # Example: info_dict['svdb_origin'] = 'new_value'
            info_dict[SVDB_ORIGIN_KEY] = reduce_to_set_of_unique_callers(
                info_dict[SVDB_ORIGIN_KEY], separator=SVDB_ORIGIN_SEPARATOR
            )

        # Reconstruct the INFO field while retaining the order
        new_info = []
        for key, value in info_dict.items():
            if value is True:
                new_info.append(f"{key}")
            else:
                new_info.append(f"{key}={value}")
        info_str = ";".join(new_info)

        # Reconstruct and print the modified VCF line
        new_fields = fields[:7] + [info_str] + others
        print("\t".join(new_fields))

    f.close()


def reduce_to_set_of_unique_callers(svdb_annotation: str, separator: str) -> str:
    callers = svdb_annotation.split(separator)
    parsed_callers = set()

    for caller in callers:
        matched_caller = Callers.match_svdb_annotation_with_caller(caller)

        # Add caller as-is if not supported
        if matched_caller is None:
            parsed_callers.add(caller)
            continue

        parsed_callers.add(matched_caller)

    # Sort the diffs to avoid random/order diff between runs.
    sorted_parsed_callers = sorted(parsed_callers)
    return separator.join(sorted_parsed_callers)


if __name__ == "__main__":
    main()
