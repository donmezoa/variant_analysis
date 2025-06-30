import requests
import time
import argparse


# Constants
GENE_LIST_FILE = 'genes.txt'
OUTPUT_FILE = 'snps_by_gene.txt'
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "omer.donmez@cchmc.org"
API_KEY = "b5818d5ea7bb82bc40bed45a681211208e08"
RATE_LIMIT_DELAY = 0.11



def read_genes(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("Gene")]

def fetch_gene_id(gene_name):
    params = {
        "db": "gene",
        "term": f"{gene_name}[Gene Name] AND 9606[Taxonomy ID]",
        "retmode": "json",
        "email": EMAIL,
        "api_key": API_KEY
    }
    r = requests.get(NCBI_BASE_URL + "esearch.fcgi", params=params)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None

def fetch_gene_coordinates(gene_id):
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "email": EMAIL,
        "api_key": API_KEY
    }
    r = requests.get(NCBI_BASE_URL + "esummary.fcgi", params=params)
    r.raise_for_status()
    doc = r.json()
    gene_data = doc["result"].get(gene_id)
    if not gene_data:
        return None

    genome = gene_data["genomicinfo"][0]
    chr_num = genome["chraccver"]
    chr_start = genome["chrstart"]
    chr_end = genome["chrstop"]

    # Convert RefSeq to chromosome number (e.g., NC_000017.11 → 17)
    chrom_map = {
        "NC_000001": "1", "NC_000002": "2", "NC_000003": "3", "NC_000004": "4",
        "NC_000005": "5", "NC_000006": "6", "NC_000007": "7", "NC_000008": "8",
        "NC_000009": "9", "NC_000010": "10", "NC_000011": "11", "NC_000012": "12",
        "NC_000013": "13", "NC_000014": "14", "NC_000015": "15", "NC_000016": "16",
        "NC_000017": "17", "NC_000018": "18", "NC_000019": "19", "NC_000020": "20",
        "NC_000021": "21", "NC_000022": "22", "NC_000023": "X", "NC_000024": "Y"
    }

    refseq_prefix = chr_num.split('.')[0]
    chr_number = chrom_map.get(refseq_prefix, refseq_prefix)
    start = min(chr_start, chr_end)
    end = max(chr_start, chr_end)
    return chr_number, start, end

def fetch_snps_in_range(chrom, start, end):
    term = f"{chrom}[CHR] AND {start}:{end}[CHRPOS]"
    base_params = {
        "db": "snp",
        "term": term,
        "retmax": 5000,
        "retmode": "json",
        "email": EMAIL,
        "api_key": API_KEY
    }
    
    all_snps = []
    retstart = 0
    total = None

    while True:
        params = base_params.copy()
        params["retstart"] = retstart

        r = requests.get(NCBI_BASE_URL + "esearch.fcgi", params=params)
        r.raise_for_status()
        data = r.json()["esearchresult"]

        if total is None:
            total = int(data["count"])
            print(f"  Total SNPs reported: {total}")

        batch = data.get("idlist", [])
        if not batch:
            break

        all_snps.extend(batch)
        retstart += len(batch)

        if len(all_snps) >= total:
            break

        time.sleep(RATE_LIMIT_DELAY)

    return all_snps

def fetch_snp_metadata_batch(rsid_list):
    summaries = []
    chunk_size = 300
    for i in range(0, len(rsid_list), chunk_size):
        chunk = rsid_list[i:i+chunk_size]
        ids = ",".join(chunk)
        data = {
            "db": "snp",
            "id": ids,
            "retmode": "json",
            "email": EMAIL,
            "api_key": API_KEY
        }
        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }
        r = requests.post(NCBI_BASE_URL + "esummary.fcgi", data=data, headers=headers)
        r.raise_for_status()
        results = r.json().get("result", {})

        for rsid in chunk:
            if rsid in results:
                record = results[rsid]
                snp_info = {
                    "rsID": f"rs{rsid}",
                    "chr": record.get("chr", "NA"),
                    "position": record.get("chrpos", "NA"),
                    "variant_class": record.get("snp_class", "NA"),
                    "fxn_class": ",".join(record.get("fxn_class", [])) if isinstance(record.get("fxn_class"), list) else record.get("fxn_class", "NA"),
                }
                summaries.append(snp_info)
        time.sleep(RATE_LIMIT_DELAY)
    return summaries





def parse_args():
    parser = argparse.ArgumentParser(description="Fetch SNPs for genes using dbSNP and NCBI APIs.")
    parser.add_argument("--genes", type=str, default="genes.txt", help="Path to gene list file.")
    parser.add_argument("--output", type=str, default="snps_by_gene.txt", help="Path to output file.")
    parser.add_argument("--limit", type=int, default=None, help="Optional maximum number of SNPs per gene.")
    parser.add_argument("--rsid-only", action="store_true", help="If set, only output rsIDs (one per line).")
    return parser.parse_args()




def main():
    args = parse_args()
    genes = read_genes(args.genes)

    with open(args.output, 'w') as out:
        if not args.rsid_only:
            out.write("Gene\trsID\tChr\tPosition\tType\tFunction\n")

        for gene in genes:
            print(f"Processing gene: {gene}")
            try:
                gene_id = fetch_gene_id(gene)
                if not gene_id:
                    print(f"  Gene ID not found.")
                    continue
                coord = fetch_gene_coordinates(gene_id)
                if not coord:
                    print(f"  Coordinates not found.")
                    continue
                chrom, start, end = coord
                snps = fetch_snps_in_range(chrom, start, end)
                print(f"  Found {len(snps)} SNPs in chr{chrom}:{start}-{end}")

                if args.limit:
                    snps = snps[:args.limit]
                    print(f"  Limited to {args.limit} SNPs.")

                if args.rsid_only:
                    for rsid in snps:
                        out.write(f"rs{rsid}\n")
                    continue

                snp_info = fetch_snp_metadata_batch(snps)
                for info in snp_info:
                    out.write(f"{gene}\t{info['rsID']}\t{info['chr']}\t{info['position']}\t{info['variant_class']}\t{info['fxn_class']}\n")

            except Exception as e:
                print(f"  Error processing {gene}: {e}")
                time.sleep(2)

if __name__ == "__main__":
    main()