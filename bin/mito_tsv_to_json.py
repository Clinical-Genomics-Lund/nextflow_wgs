#!/usr/bin/env python3

import argparse
import csv
import json
import sys
from typing import Dict

SCRIPT_DESCRIPTION = (
    "Take TSV output generated by process sentieon_mitochondrial_qc and convert "
    " to JSON for later merge with other sample QC JSON data."
)


AVERAGE_COVERAGE_FIELD_SUFFIX = "_mean_cvg"
PCT_ABOVE_500X_FIELD_SUFFIX = "_%_above_500"

# Parse command-line arguments
parser = argparse.ArgumentParser(description=SCRIPT_DESCRIPTION)
parser.add_argument("mito_coverage_tsv", help="Path to the input TSV file")
args = parser.parse_args()


def extract_value_by_suffix(row: Dict[str, str], key_suffix: str) -> str:
    """
    Iterate over dict and pull out values based on key/field suffix

    The fields are prefixed with a sample-name derived from the BAM-file,
    which can deviate from the group/sample ids specified in the input CSV
    e.g. when the pipeline is started from an old BAM but with new sample ids
    specified in the CSV.

    The TSVs are generated per-sample, meaning there is no risk of the same
    TSV containing multiple samples, which would require a more precise
    dict access.
    """

    keys_with_matching_suffix = [key for key in row.keys() if key.endswith(key_suffix)]

    if len(keys_with_matching_suffix) == 0:
        raise KeyError(f"No field in coverage TSV with suffix {key_suffix}.")
    if len(keys_with_matching_suffix) > 1:
        raise KeyError(
            f"Multiple fields in row with suffix {key_suffix}: {keys_with_matching_suffix}"
        )

    return row[keys_with_matching_suffix[0]]


# Extract relevant data for correct sample:
with open(args.mito_coverage_tsv, "r") as tsv_file:
    reader = csv.DictReader(tsv_file, delimiter="\t")
    row = next(reader)

average_coverage = extract_value_by_suffix(row, AVERAGE_COVERAGE_FIELD_SUFFIX)
pct_above_500x = extract_value_by_suffix(row, PCT_ABOVE_500X_FIELD_SUFFIX)

json_data = {
    "mito_qc": {
        "mean_coverage": float(average_coverage),
        "pct_above_x": {"500": float(pct_above_500x)},
    }
}

# And spit it (std.)out:
json.dump(json_data, sys.stdout, indent=4)
