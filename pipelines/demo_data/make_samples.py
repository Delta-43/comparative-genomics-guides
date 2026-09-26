#!/usr/bin/env python3
"""Turn an ENA run table (the *.ena.tsv files next to this script) into
   (1) download.sh  - resumable curl downloads + md5 verification
   (2) samples.tsv  - ready for pipelines/rnaseq, chip_atac or hic

Usage:  python make_samples.py E-MTAB-13253_ATAC-seq.ena.tsv --pipeline chip_atac --outdir /data/ciuba2025/atac
Refresh a table: curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP150225&result=read_run&fields=run_accession,sample_alias,library_layout,read_count,base_count,fastq_md5,fastq_ftp&format=tsv"

One row is written per sequencing RUN (sample = <line>_<run accession>) so no data is silently merged.
Several runs of the same library are re-sequencing (technical replicates): sum their counts in R
(RNA-seq), merge BAMs (the chip_atac `group` column pools them for peak calling), or merge .hic
files with Juicer's mega.sh (Hi-C).
"""
import argparse, csv, os, re, sys

def load_species_map(path):
    """TSV: regex<TAB>species<TAB>genome. The first regex matching the sample alias wins."""
    rows = []
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            pat, sp, g = line.rstrip("\n").split("\t")[:3]
            rows.append((pat, sp, g))
    return rows

def species_of(line, species_map):
    for pat, sp, g in species_map:
        if re.match(pat, line, re.I):
            return sp, g
    sys.exit(f"no species-map entry matches sample '{line}': add a line to your --species-map file")

ap = argparse.ArgumentParser()
ap.add_argument("ena_tsv")
ap.add_argument("--pipeline", required=True, choices=["rnaseq", "chip_atac", "hic"])
ap.add_argument("--outdir", required=True, help="where FASTQs will be downloaded")
ap.add_argument("--species-map", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "species_map.tsv"),
                help="TSV of regex, species, genome used to label samples (default: species_map.tsv next to this script)")
ap.add_argument("--min-reads", type=float, default=5e6, help="skip shallow runs, e.g. small QC libraries (default 5M read pairs)")
a = ap.parse_args()
SPECIES = load_species_map(a.species_map)

rows = list(csv.DictReader(open(a.ena_tsv), delimiter="\t"))
os.makedirs(a.outdir, exist_ok=True)
dl = open(os.path.join(a.outdir, "download.sh"), "w")
dl.write("#!/usr/bin/env bash\nset -euo pipefail\ncd \"$(dirname \"$0\")\"\n")
md5 = open(os.path.join(a.outdir, "md5sums.txt"), "w")
out = open(os.path.join(a.outdir, "samples.tsv"), "w")
hdr = {"rnaseq": "sample\tfastq_1\tfastq_2\tline\tspecies",
       "chip_atac": "sample\tfastq_1\tfastq_2\tgroup\tcontrol",
       "hic": "sample\tfastq_1\tfastq_2\tgenome"}[a.pipeline]
out.write(hdr + "\n")
kept = skipped = 0
for r in sorted(rows, key=lambda r: (r["sample_alias"], r["run_accession"])):
    if r["library_layout"] != "PAIRED" or int(r["read_count"] or 0) < a.min_reads:
        skipped += 1
        continue
    line = re.sub(r"[^A-Za-z0-9]+", "", r["sample_alias"].split(":", 1)[-1].split("_")[0])
    sp, genome = species_of(line, SPECIES)
    rep = "_".join(r["sample_alias"].split(":", 1)[-1].split("_")[1:]) or "1"
    name = f"{sp}_{line}_{re.sub(r'[^A-Za-z0-9]+', '', rep)}_{r['run_accession']}"
    urls, sums = r["fastq_ftp"].split(";"), r["fastq_md5"].split(";")
    fqs = []
    for i, (u, m) in enumerate(zip(urls, sums), start=1):
        fq = f"{name}_R{i}.fastq.gz"
        dl.write(f"curl -sS -C - -o {fq} https://{u}\n")
        md5.write(f"{m}  {fq}\n")
        fqs.append(os.path.join(os.path.abspath(a.outdir), fq))
    extra = {"rnaseq": f"{line}\t{sp}", "chip_atac": f"{sp}\t", "hic": genome}[a.pipeline]
    out.write(f"{name}\t{fqs[0]}\t{fqs[1]}\t{extra}\n")
    kept += 1
dl.write("md5sum -c md5sums.txt\n")
print(f"{kept} runs written to {a.outdir}/samples.tsv ({skipped} skipped: single-end or < {a.min_reads:.0f} reads)")
print(f"download with: bash {a.outdir}/download.sh")
