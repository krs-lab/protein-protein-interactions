#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Extract and rank ipTM, pTM, and pLDDT scores from ColabFold rank_001 output JSON files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 get_ipTM.py -i output_predictions
  python3 get_ipTM.py --input output_predictions/
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        type=str,
        help="Path to the directory containing ColabFold prediction output JSON files."
    )
    
    args = parser.parse_args()

    # Clean input directory path to handle trailing slashes seamlessly
    input_dir = os.path.normpath(args.input)

    if not os.path.isdir(input_dir):
        print(f"Error: Directory '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Search for rank_001 score JSON files
    search_pattern = os.path.join(input_dir, "*_scores_rank_001_*.json")
    files = glob.glob(search_pattern)

    if not files:
        print(f"No rank_001 JSON score files found in '{input_dir}'.", file=sys.stderr)
        sys.exit(1)

    results = []

    for f in files:
        with open(f, 'r') as fn:
            data = json.load(fn)
            name = os.path.basename(f).split('_scores_rank_001')[0]
            iptm = data.get('iptm', 0.0)
            ptm = data.get('ptm', 0.0)
            
            plddt_list = data.get('plddt', [])
            plddt = sum(plddt_list) / len(plddt_list) if plddt_list else 0.0

            results.append({
                'complex': name,
                'iptm': iptm,
                'ptm': ptm,
                'plddt': plddt
            })

    # Sort results by ipTM in descending order
    results.sort(key=lambda x: x['iptm'], reverse=True)

    # Print terminal output
    print(f'{"Complex":<40} | {"ipTM":<8} | {"pTM":<8} | {"pLDDT":<8}')
    print('-' * 72)
    for res in results:
        print(f"{res['complex']:<40} | {res['iptm']:<8.3f} | {res['ptm']:<8.3f} | {res['plddt']:<8.1f}")

    # Write sorted results to CSV
    output_csv = os.path.join(input_dir, "summary_iptm_scores.csv")
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Complex', 'ipTM', 'pTM', 'pLDDT']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for res in results:
            writer.writerow({
                'Complex': res['complex'],
                'ipTM': f"{res['iptm']:.3f}",
                'pTM': f"{res['ptm']:.3f}",
                'pLDDT': f"{res['plddt']:.1f}"
            })

    print(f"\nSaved sorted CSV report to: {output_csv}")

if __name__ == "__main__":
    main()
