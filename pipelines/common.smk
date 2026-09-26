# Shared helpers for every pipeline in this folder.
import csv, os

def read_samples(path):
    """samples.tsv -> {sample: {column: value}}. Lines starting with '#' are ignored.
    Relative paths in columns named fastq_* / bam* are resolved against the
    directory holding samples.tsv, so a sample sheet can live next to its data."""
    samples = {}
    base = os.path.dirname(os.path.abspath(path))
    with open(path) as fh:
        rows = csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")
        for r in rows:
            r = {k.strip(): (v or "").strip() for k, v in r.items()}
            for k, v in r.items():
                if v and k.startswith(("fastq", "bam")) and not os.path.isabs(v):
                    r[k] = os.path.join(base, v)
            if r["sample"] in samples:
                raise ValueError(f"duplicate sample name in {path}: {r['sample']}")
            samples[r["sample"]] = r
    if not samples:
        raise ValueError(f"no samples found in {path}")
    return samples

def threads(name, default=4):
    return int(config.get("threads", {}).get(name, config.get("threads", {}).get("default", default)))

def mem_mb(name, default=8000):
    return int(config.get("mem_mb", {}).get(name, default))
