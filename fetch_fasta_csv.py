#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
===============================================================================
Script Name: fetch_proteins_from_csv.py
Description:
    Reads a CSV file containing UniProt IDs (or accession numbers) and fetches
    their canonical FASTA amino acid sequences directly from the UniProt REST API.

Features:
    - Automatically detects variations of UniProt ID column headers:
      ['uniprot', 'UNIPROT', 'uniprot_id', 'uniprot_ids', 'UNIPROT_ID', 'UNIPROT_IDs']
    - Saves individual `.fasta` files for each protein under the output directory.
    - Consolidates all retrieved sequences into a single `multi_fasta.fasta` file.
    - Built-in `--help` section explaining usage and CLI arguments.
    - Tab completion support via argcomplete for input/output path arguments.

Usage Examples:
    python3 fetch_proteins_from_csv.py -i legionella_effectors_368.csv
    python3 fetch_proteins_from_csv.py --input targets.csv --output my_proteins
===============================================================================
"""

import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request

try:
    import argcomplete
except ImportError:
    argcomplete = None


def fetch_fasta_by_uniprot_id(uniprot_id):
    """Fetch the canonical FASTA sequence directly using a UniProt Accession ID."""
    uniprot_id = uniprot_id.strip()
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    fasta_req = urllib.request.Request(
        fasta_url, headers={"User-Agent": "Python-UniProt-CSV-Fetcher"}
    )

    try:
        with urllib.request.urlopen(fasta_req) as fasta_resp:
            fasta_data = fasta_resp.read().decode("utf-8")
        return uniprot_id, fasta_data
    except Exception as e:
        print(f"[ERROR] Exception during FASTA retrieval for {uniprot_id}: {e}")
        return uniprot_id, None


def main():
    parser = argparse.ArgumentParser(
        description="Fetch protein FASTA sequences from UniProtKB using a CSV file containing UniProt IDs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 fetch_proteins_from_csv.py -i legionella_effectors_368.csv
  python3 fetch_proteins_from_csv.py --input targets.csv --output custom_dir
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="csv_file",
        default="targets.csv",
        help="Path to the input CSV file containing UniProt IDs (default: targets.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        default="proteins",
        help="Directory to store output FASTA files (default: proteins)",
    )

    if argcomplete:
        argcomplete.autocomplete(parser)

    args = parser.parse_args()
    csv_file = args.csv_file
    output_dir = args.output_dir

    if not os.path.exists(csv_file):
        print(f"Error: Could not find CSV file '{csv_file}'.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    multi_fasta_path = os.path.join(output_dir, "multi_fasta.fasta")

    # Clear/initialize the multi_fasta.fasta file
    with open(multi_fasta_path, "w", encoding="utf-8") as mf:
        mf.write("")

    print(f"Reading targets from {csv_file}...\n")

    valid_uniprot_headers = {
        "uniprot",
        "UNIPROT",
        "uniprot_id",
        "uniprot_ids",
        "UNIPROT_ID",
        "UNIPROT_IDs",
    }

    successful_count = 0
    failed_count = 0
    with open(csv_file, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        # Identify which UniProt ID column exists in the CSV
        uniprot_col = None
        if reader.fieldnames:
            for field in reader.fieldnames:
                if field.strip() in valid_uniprot_headers:
                    uniprot_col = field
                    break

        if not uniprot_col:
            print(
                f"[ERROR] No valid UniProt ID column found in {csv_file}."
            )
            print(
                f"Accepted column names: {', '.join(sorted(valid_uniprot_headers))}"
            )
            sys.exit(1)

        for row in reader:
            uniprot_id = row[uniprot_col].strip()
            if not uniprot_id:
                continue

            out_filename = f"{uniprot_id}.fasta"
            out_path = os.path.join(output_dir, out_filename)

            print(f"Fetching FASTA for UniProt ID: {uniprot_id} ...")
            accession, fasta_seq = fetch_fasta_by_uniprot_id(uniprot_id)

            if fasta_seq:
                # Ensure fasta content ends with a newline
                formatted_fasta = (
                    fasta_seq if fasta_seq.endswith("\n") else fasta_seq + "\n"
                )

                # Save single FASTA file
                with open(out_path, "w", encoding="utf-8") as out_f:
                    out_f.write(formatted_fasta)

                # Append to multi_fasta.fasta file
                with open(multi_fasta_path, "a", encoding="utf-8") as mf:
                    mf.write(formatted_fasta)

                successful_count += 1
                print(
                    f" -> [SUCCESS] UniProt ID: {accession} Saved -> {out_path}\n"
                )
            else:
                failed_count += 1
                print(
                    f" -> [FAILED] Could not retrieve FASTA sequence for UniProt ID: {uniprot_id}\n"
                )

    print("-" * 50)
    print(
        f"Done! Aggregated {successful_count} sequences into '{multi_fasta_path}' ({failed_count} failed)."
    )


if __name__ == "__main__":
    main()
