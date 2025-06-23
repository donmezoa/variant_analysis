import requests
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper
from hgvs.exceptions import HGVSDataNotAvailableError
import csv
import argparse
import re
import sys

def fetch_variant_from_rs(rsid):
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid[2:]}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Error fetching data for {rsid}: {response.status_code}")
    
    data = response.json()
    hgvs_list = []

    try:
        placements = data['primary_snapshot_data']['placements_with_allele']
        for placement in placements:
            for allele in placement.get('alleles', []):
                hgvs_values = allele.get('hgvs', [])
                if isinstance(hgvs_values, str):
                    hgvs_list.append(hgvs_values)
                elif isinstance(hgvs_values, list):
                    hgvs_list.extend(hgvs_values)
    except (KeyError, IndexError):
        pass

    return hgvs_list

def filter_hgvs(hgvs_list, level="all"):
    if level == "all":
        return [expr for expr in hgvs_list if "=" not in expr]
    elif level == "transcript":
        return [expr for expr in hgvs_list if expr.startswith("NM_") and "=" not in expr]
    elif level == "protein":
        return [expr for expr in hgvs_list if expr.startswith("NP_") and "=" not in expr]
    elif level == "genomic":
        return [expr for expr in hgvs_list if expr.startswith("NC_") and "=" not in expr]
    else:
        raise ValueError(f"Unknown level: {level}")

def hgvs_to_igv(hgvs_expr):
    match = re.match(r"NC_0{0,5}0*(\d+)\.\d+:g\.(\d+)", hgvs_expr)
    if match:
        chrom = f"chr{int(match.group(1))}"
        pos = int(match.group(2))
        start = max(1, pos - 5)
        end = pos + 5
        return f"{chrom}:{start}-{end}"
    return ""

def classify_variant_type(hgvs_expr):
    if 'delins' in hgvs_expr:
        return "Complex Substitution"
    elif 'del' in hgvs_expr:
        return "Deletion"
    elif 'ins' in hgvs_expr:
        return "Insertion"
    elif 'dup' in hgvs_expr:
        return "Duplication"
    elif ':' in hgvs_expr and '>' in hgvs_expr:
        return "SNV"
    elif 'fs' in hgvs_expr or 'frameshift' in hgvs_expr:
        return "Frameshift"
    elif ':' in hgvs_expr and '.' in hgvs_expr:
        return "Substitution"
    else:
        return "Other"


def process_file(input_file, output_file, level="all", include_igv=False):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        headers = ["rsID", "HGVS", "Variant_Type"]
        if include_igv:
            headers.append("IGV_Region")
        writer.writerow(headers)

        for line in infile:
            rsid = line.strip()
            if not rsid:
                continue
            try:
                hgvs_exprs = fetch_variant_from_rs(rsid)
                filtered_exprs = filter_hgvs(hgvs_exprs, level=level)

                for expr in sorted(set(filtered_exprs)):
                    variant_type = classify_variant_type(expr)
                    row = [rsid, expr, variant_type]
                    if include_igv:
                        igv_region = hgvs_to_igv(expr) if expr.startswith("NC_") else ""
                        row.append(igv_region)
                    writer.writerow(row)
            except Exception as e:
                row = [rsid, f"ERROR: {str(e)}", "", "", "", "", ""]
                if include_igv:
                    row.insert(3, "")  # empty IGV field
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Convert rsIDs to HGVS with classification, and IGV region")
    parser.add_argument("input_file", help="File with rsIDs (one per line)")
    parser.add_argument("output_file", help="CSV file to write results")
    parser.add_argument("--level", choices=["all", "transcript", "protein", "genomic"], default="all",
                        help="Type of HGVS expression to include")
    parser.add_argument("--igv", action="store_true", help="Include IGV region (only valid for 'all' or 'genomic')")

    args = parser.parse_args()

    if args.igv and args.level not in {"all", "genomic"}:
        print("Error: --igv option is only valid when --level is 'all' or 'genomic'.", file=sys.stderr)
        sys.exit(1)

    process_file(args.input_file, args.output_file, level=args.level, include_igv=args.igv)

if __name__ == "__main__":
    main()

