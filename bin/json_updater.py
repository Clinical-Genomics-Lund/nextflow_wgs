#!/usr/bin/env python3

import argparse
import json
import sys
from typing import Optional

"""
Author:      Alexander Koc <alexander.koc@skane.se> (with lots of help from ChatGPT!)
Description: Add or set values in .JSON files from the command line.
Project URL: https://github.com/alkc/json-updater
"""

__version__ = "1.0.1"


def update_json_file(
    json_file: str, data: dict, parent_key: Optional[str], edit_file_in_place: bool
):
    # Load existing JSON data from the file
    with open(json_file, "r") as file:
        json_data = json.load(file)

    # If parent_key is provided, nest the new key-value pairs under it
    if parent_key:
        # Create the parent key if it doesn't exist
        if parent_key not in json_data:
            json_data[parent_key] = {}
        elif not type(json_data[parent_key]) is dict:
            print(
                f"ERROR: Cannot assign values to parent_key '{parent_key}' as it already exists "
                "in input json, but is not an object/dict!",
                file=sys.stderr,
            )
            sys.exit(1)

        # Update the JSON data with the new key-value pairs under the parent key
        for key, value in data.items():
            json_data[parent_key][key] = value
    else:
        # Update the JSON data with the new key-value pairs directly
        for key, value in data.items():
            json_data[key] = value

    # Write the updated JSON data back to the file or print to stdout
    if edit_file_in_place:
        with open(json_file, "w") as file:
            json.dump(json_data, file, indent=4)
            print(f"Updated JSON file: {json_file}")
    else:
        json.dump(json_data, sys.stdout, indent=4)
        print(file=sys.stdout)  # Print a new line


def main():
    parser = argparse.ArgumentParser(description="Update JSON file with key-value pairs")
    parser.add_argument("json_file", help="path to the JSON file")
    parser.add_argument(
        "-d",
        "--data",
        nargs="+",
        metavar="key=value",
        help="key-value pairs to add/update in the JSON file",
    )
    parser.add_argument(
        "--parent-key", help="parent key under which new key-value pairs will be nested"
    )

    parser.add_argument(
        "-i",
        "--inplace",
        action="store_true",
        help="modify the JSON file in place instead of printing to stdout",
    )

    args = parser.parse_args()

    if not args.data:
        parser.error("Please provide key-value pairs.")

    data = {}
    for item in args.data:
        key, value = item.split("=")
        data[key] = value

    update_json_file(args.json_file, data, args.parent_key, args.inplace)


if __name__ == "__main__":
    main()
