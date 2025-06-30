These are programs to help with variant analysis. 

I found all the genes that MGPS uses for CinCSeq Comprehensive Cancer Panel. 
https://www.cincinnatichildrens.org/clinical-labs/our-labs/pathology/molecular-genomic-pathology-services/-/media/C99B1C995FC8428DA7080059BBF3AD70.ashx
I named the file mgps_genes.txt



I made a program to find all the SNPs in all the genes. For the example and the example uploaded to GitHub, I limited it to 2 SNPs per gene for size and speed considerations. Also, I saved the output in a text file.

```{}
python gene2snp.py --genes mgps_genes.txt --output snps_by_gene.txt --limit 2 > output.txt
```

I then only wanted to save the rsid into a single text file, to use for further analysis.

```{}
python gene2snp.py --genes mgps_genes.txt --output rsID.txt --limit 2 --rsid-only
```

Then I used the rsID.txt file as input for the rsid_to_hgvs.py script to convert the rsID to HGVS nomenclature.

```{}
python rsid_to_hgvs.py rsID.txt hgvs_output.txt --level genomic --igv
```

Help Command
```{}
python gene2snp.py --help
```
Help Output
```{}
usage: gene2snp.py [-h] [--genes GENES] [--output OUTPUT] [--limit LIMIT] [--rsid-only]

Fetch SNPs for genes using dbSNP and NCBI APIs.

options:
  -h, --help       show this help message and exit
  --genes GENES    Path to gene list file.
  --output OUTPUT  Path to output file.
  --limit LIMIT    Optional maximum number of SNPs per gene.
  --rsid-only      If set, only output rsIDs (one per line).
```



Help Command 
```{}
python rsid_to_hgvs.py --help
```
Help Output
```{}
usage: rsid_to_hgvs.py [-h] [--level {all,transcript,protein,genomic}] [--igv] input_file output_file

Convert rsIDs to HGVS with classification, and IGV region

positional arguments:
  input_file            File with rsIDs (one per line)
  output_file           CSV file to write results

options:
  -h, --help            show this help message and exit
  --level {all,transcript,protein,genomic}
                        Type of HGVS expression to include
  --igv                 Include IGV region (only valid for 'all' or 'genomic')
```




