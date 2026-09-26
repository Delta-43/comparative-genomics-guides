---
title: "Hi-C: 3D genome organisation"
subtitle: "Contact maps, compartments, domains, loops and how they differ"
---

::: {.callout-tip}
## What this guide covers
- What a Hi-C contact map measures, and how to plan depth for the features of interest.
- Processing reads into a `.hic` file with Juicer, and judging library quality.
- Reading and normalising contact matrices in R; plotting maps and distance-decay curves.
- Calling and interpreting compartments, TADs and loops, and validating loops with aggregate analysis.
- Comparing 3D structure between conditions or species, and the limits of each comparison.

**Time:** about 4 hours. **You need:** the set-up from [Getting started](00-getting-started.md);
a large workstation or HPC for processing (Hi-C libraries are often billions of reads).
:::

## 1. What Hi-C measures

Hi-C crosslinks chromatin, cuts it with a restriction enzyme, ligates fragments that were close
together in the nucleus, and sequences the junctions
([Lieberman-Aiden *et al.* 2009](https://doi.org/10.1126/science.1181369)). Each read pair is one
**contact** between two loci. Binning the genome (e.g. 5 kb) turns the contacts into a matrix.

Features you can see at increasing resolution:

| Feature | Scale | What it is |
|---|---|---|
| A/B compartments | 1–100 Mb (analysed at 100–500 kb bins) | Active (A) and inactive (B) chromatin that preferentially contact their own kind |
| TADs / contact domains | 100 kb–2 Mb (analysed at 10–50 kb) | Regions that contact themselves more than their neighbours ([Dixon *et al.* 2012](https://doi.org/10.1038/nature11082)) |
| Loops | Anchors 5–10 kb, spans ≤ a few Mb | Point-like enrichments, mostly between convergent CTCF sites ([Rao *et al.* 2014](https://doi.org/10.1016/j.cell.2014.11.021)) |

Loop extrusion by cohesin, stopped at CTCF sites facing each other, explains most loops and
domains ([Fudenberg *et al.* 2016](https://doi.org/10.1016/j.celrep.2016.04.085)).

**It can tell you:** which loci are in contact more often than expected, in a population average.
**It can't tell you:** whether a contact happens in all cells or a few, whether it's functional,
or the direction of causality between structure and expression.

## 2. Experimental design

| Factor | Recommendation |
|---|---|
| Replicates | ≥ 2 biological replicates. Merge them for feature calling, keep them separate for testing |
| Depth | Compartments: ~50–100 M valid pairs. TADs: a few hundred million. Loops at 5 kb: ~1–2 billion ([Rao *et al.* 2014](https://doi.org/10.1016/j.cell.2014.11.021)) |
| Enzyme | 4-cutters (MboI/DpnII) or multi-enzyme kits give finer resolution than 6-cutters (HindIII) |
| Protocol | *In situ* Hi-C (ligation in intact nuclei) has much less random-ligation noise than dilution Hi-C |
| Check first | Shallow-sequence each library (~10–20 M pairs) and check QC before sequencing deeply |

## 3. From FASTQ to a `.hic` file

```bash
cd pipelines/hic
$EDITOR config.yaml samples.tsv     # genome FASTA, chromosome list, enzyme
snakemake --use-conda --cores 16 -n && snakemake --use-conda --cores 16
```

The pipeline runs [Juicer](https://github.com/aidenlab/juicer) ([Durand *et al.* 2016](https://doi.org/10.1016/j.cels.2016.07.002))
in CPU mode:

1. BWA aligns each read end separately. Chimeric reads that span the ligation junction are
   split and kept.
2. Each pair is assigned to restriction fragments. Pairs that can't be real ligation products
   (e.g. both ends in one fragment) are removed.
3. PCR duplicates are removed.
4. A `.hic` file is built from pairs where **both ends have MAPQ ≥ 30** (`inter_30.hic`), with
   multiple resolutions and balancing normalisations stored inside.

It then dumps raw contacts per chromosome at 5 kb for R, calls TopDom domains at 25 kb, and calls
loops with HiCCUPS or Mustache.

::: {.callout-note collapse="true"}
## Alternative toolchain: pairtools + cooler
The [Open2C](https://open2c.github.io/) ecosystem does the same job in Python: `bwa mem -SP5M` →
`pairtools parse/sort/dedup` → `cooler cload` → `.mcool` ([Abdennur & Mirny 2020](https://doi.org/10.1093/bioinformatics/btz540)).
Analyses then use `cooltools` (expected, compartments, insulation, pile-ups). **Pros:**
scriptable, fast, rich analysis library. **Cons:** a different file format (`.mcool` vs `.hic`),
though `hic2cool` converts between them. Many analyses (and Juicebox visualisation) assume `.hic`.
:::

### Library QC

Juicer writes statistics to `results/hic/<sample>.inter_30.txt`; the `inter.txt` in each sample's
`aligned/` folder covers all MAPQ values. What to look at:

| Statistic | Good sign | Bad sign |
|---|---|---|
| Unique alignable pairs | High fraction of sequenced pairs | Low: adapter problems, contamination |
| Duplicates | Low in the shallow run and rising slowly with depth | High early: low complexity, so deeper sequencing is wasted |
| Intra- vs inter-chromosomal | Most contacts intra-chromosomal | High inter-chromosomal fraction: random ligation noise |
| Long-range intra (> 20 kb) | A large share of intra-chromosomal contacts | Mostly < 20 kb: undigested or unligated DNA |
| Map resolution | Smallest bin where ≥ 80 % of bins have ≥ 1,000 contacts (Rao 2014; Juicer `calculate_map_resolution.sh`) | Coarser than the features you want to study |

## 4. Reading contact matrices in R

`pipelines/hic` dumps **raw** counts (`observed NONE`) as sparse triplets: bin1 start, bin2
start, count, upper triangle only. Reading one chromosome into a sparse matrix:

```r
library(Matrix); library(data.table)
read_dump <- function(file, chrom_length, binsize) {
  d <- fread(file, col.names = c("x", "y", "count"))
  n <- ceiling(chrom_length / binsize)
  m <- sparseMatrix(i = d$x / binsize + 1, j = d$y / binsize + 1, x = d$count, dims = c(n, n))
  m + t(m) - Diagonal(n, diag(m))                  # symmetrise without doubling the diagonal
}
m <- read_dump("results/dump_5kb/ctrl_rep1/ctrl_rep1_inter_30_chr21_NO_5KB_dump.txt",
               chrom_length = 46709983, binsize = 5000)
```

For quick exploration without dumping, [strawr](https://github.com/aidenlab/straw) reads
`.hic` files directly: `strawr::straw("KR", "sample.hic", "21", "21", "BP", 25000)`.

### Distance decay

Contact frequency falls with genomic distance. The slope of that curve, P(s), is a basic
property of the sample:

```r
ps <- function(m, binsize, max_bins = 2000) {
  d <- summary(m); d <- d[d$j >= d$i & (d$j - d$i) <= max_bins, ]
  agg <- aggregate(x ~ I(j - i), data = d, FUN = sum)
  names(agg) <- c("bins", "contacts")
  agg$distance <- agg$bins * binsize
  agg$p <- agg$contacts / (nrow(m) - agg$bins)      # mean contacts per bin pair at that distance
  agg[agg$bins > 0, ]
}
p <- ps(m, 5000)
plot(p$distance, p$p, log = "xy", type = "l", xlab = "genomic distance (bp)", ylab = "contact probability")
```

Replicates should give overlapping curves. A different slope (e.g. between cell-cycle stages or
after cohesin loss) is itself a biological result, and must be accounted for before comparing
individual contacts.

### Normalisation (matrix balancing)

Raw counts carry biases that aren't about 3D structure: restriction-site density, GC content,
mappability. **Matrix balancing** assumes each locus should, overall, be equally visible. It
iteratively rescales rows and columns until they all have the same sum (ICE,
[Imakaev *et al.* 2012](https://doi.org/10.1038/nmeth.2148); Knight–Ruiz "KR" and Juicer's
"SCALE" solve the same problem differently):

```r
ice <- function(m, iterations = 200, min_coverage = 0.01) {
  m <- as(m, "CsparseMatrix"); n <- nrow(m)
  cov <- rowSums(m)
  keep <- cov > quantile(cov[cov > 0], min_coverage)          # drop empty / very sparse bins
  mask <- Diagonal(x = as.numeric(keep))
  mm <- mask %*% m %*% mask                                   # dropped bins set to 0, dimensions kept
  b <- rep(1, n)
  for (k in seq_len(iterations)) {
    s <- rowSums(mm); s[keep] <- s[keep] / mean(s[keep]); s[!keep] <- 1
    f <- sqrt(s)                                              # split the correction between row and column
    mm <- Diagonal(x = 1 / f) %*% mm %*% Diagonal(x = 1 / f)
    b <- b * f
    if (max(abs(s[keep] - 1)) < 1e-3) break
  }
  if (k == iterations) warning("ICE did not converge; increase `iterations`")
  bias <- b; bias[!keep] <- NA
  list(matrix = mm, bias = bias, kept = keep, iterations = k)  # same bin indices as the input
}
bal <- ice(m)
summary(rowSums(bal$matrix)[bal$kept])      # all ~equal after balancing
```

In practice, use the normalisations Juicer stores in the `.hic` file (KR or SCALE) unless you
need your own. Always normalise each sample and chromosome separately.

## 5. Visualising contact maps

[plotgardener](https://phanstiellab.github.io/plotgardener/) reads `.hic` files and lays out
maps, gene tracks and loops on one page:

```r
library(plotgardener)
pageCreate(width = 6, height = 6, default.units = "inches")
hic <- plotHicSquare(data = "results/hic/ctrl_rep1.inter_30.hic", chrom = "chr21",
                     chromstart = 30e6, chromend = 33e6, resolution = 10000, norm = "KR",
                     x = 0.5, y = 0.5, width = 5, height = 5, default.units = "inches")
annoHeatmapLegend(plot = hic, x = 5.6, y = 0.5, width = 0.15, height = 1.5, default.units = "inches")
pageGuideHide()
```

Show the same region at the same colour scale for every condition. Differences that only
appear because one map is scaled to a brighter maximum are a common artefact in figures.

## 6. Compartments, domains and loops

**Compartments.** Juicer computes the first eigenvector of the observed/expected correlation
matrix per chromosome:

```bash
java -jar juicer_tools.jar eigenvector KR sample.inter_30.hic chr21 BP 250000 chr21_eigen.txt
```

The sign of an eigenvector is arbitrary. Orient it so that positive values correlate with gene
density or GC content (= A compartment) before comparing samples.

**Domains (TopDom).** The pipeline calls TopDom ([Shin *et al.* 2016](https://doi.org/10.1093/nar/gkv1505))
at 25 kb with `window.size = 4`, a 100 kb window on each side of each bin. Larger windows give
fewer, larger domains. Domain calls depend strongly on resolution and caller, so compare
**boundaries** between samples with a tolerance (± 1–2 bins) rather than exact domains.

A simple **insulation score** (lower = stronger boundary) is easy to compute yourself:

```r
insulation <- function(m, w = 20) {                 # w bins on each side (20 x 5 kb = 100 kb)
  n <- nrow(m); s <- rep(NA_real_, n)
  for (i in (w + 1):(n - w)) s[i] <- mean(m[(i - w):(i - 1), (i + 1):(i + w)])
  log2(s / mean(s, na.rm = TRUE))
}
ins <- insulation(bal$matrix)
```

**Loops.** HiCCUPS ([Rao *et al.* 2014](https://doi.org/10.1016/j.cell.2014.11.021)) finds
pixels enriched over several local backgrounds. The pipeline runs its CPU version
(`hiccups --cpu`), which juicer_tools labels experimental and which only searches within 8 Mb of
the diagonal. [Mustache](https://github.com/ay-lab/mustache) ([Roayaei Ardakany *et al.* 2020](https://doi.org/10.1186/s13059-020-02167-0))
is a fast, scale-space alternative that often finds more loops at lower depth. Different callers
agree on strong loops and disagree on weak ones.

### Validating loops: aggregate peak analysis (APA)

Averaging the contact map around many loop anchor pairs shows whether the set, *as a whole*, is
enriched. Always compare against a matched random control:

```r
apa <- function(m, loops_bins, flank = 10) {        # loops_bins: data.frame(i, j) of bin indices
  n <- nrow(m); acc <- matrix(0, 2 * flank + 1, 2 * flank + 1); k <- 0
  for (r in seq_len(nrow(loops_bins))) {
    i <- loops_bins$i[r]; j <- loops_bins$j[r]
    if (i - flank < 1 || j + flank > n || j - i < 2 * flank) next
    acc <- acc + as.matrix(m[(i - flank):(i + flank), (j - flank):(j + flank)]); k <- k + 1
  }
  acc / k
}
# loops on this chromosome from a BEDPE (chr1 x1 x2 chr2 y1 y2 ...), converted to 5 kb bin indices
bedpe <- read.table("results/loops/ctrl_rep1.hiccups.bedpe", comment.char = "#")
bedpe <- bedpe[bedpe$V1 %in% c("21", "chr21") & bedpe$V1 == bedpe$V4, ]
loops_bins <- data.frame(i = floor((bedpe$V2 + bedpe$V3) / 2 / 5000) + 1,
                         j = floor((bedpe$V5 + bedpe$V6) / 2 / 5000) + 1)
set.seed(1)
random <- data.frame(i = sample(20:(nrow(m) - 400), nrow(loops_bins)))
random$j <- random$i + (loops_bins$j - loops_bins$i)          # same distances as the real loops
real_apa <- apa(bal$matrix, loops_bins); rand_apa <- apa(bal$matrix, random)
c(real = real_apa[11, 11] / mean(real_apa[18:21, 1:4]),       # centre vs lower-left corner
  random = rand_apa[11, 11] / mean(rand_apa[18:21, 1:4]))
image(real_apa, main = "APA")
```

Real loops give a bright centre (ratio well above 1). The random control, at matched
distances, should be near 1. A bright random centre means a distance or normalisation artefact.

## 7. Comparing conditions

| Situation | Approach | Pros | Cons |
|---|---|---|---|
| ≥ 2–3 replicates per condition | Count-based tests on binned contacts or loop pixels: [diffHic](https://bioconductor.org/packages/diffHic/), [multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare/) (edgeR/GLM) | Real p-values, handles replicates | Needs depth and replicates; many tests at high resolution |
| Loops called separately per condition | Overlap loop sets, then test APA / pixel signal at the union of loops | Intuitive | Absence of a call ≠ absence of a loop (calling depends on depth) |
| One sample per group | Descriptive only: require a difference larger than between any replicate pair, with a random-loop control for the false-positive rate | Works with n = 1 | No error model; state results as candidates |

Before calling a loop differential, **match depth** (downsample the deeper sample, or compare
normalised signal), check P(s) curves are similar, and inspect the region in maps from every
replicate.

## 8. Connecting structure to function

- **Genes and loops/domains:** are differentially expressed genes (see the
  [RNA-seq guide](01-rna-seq.md)) enriched inside differential domains or at differential loop
  anchors? Test against all expressed genes, not the whole genome.
- **CTCF:** overlap loop anchors with CTCF ChIP-seq peaks and scan them for the CTCF motif.
  Anchors with motifs in **convergent** orientation (pointing towards each other) are the
  canonical extrusion-stopped loops, and should give the strongest APA signal.
- **Enhancers and promoters:** loops linking a differential enhancer
  ([ATAC-seq guide](02-atac-seq.md)) to the promoter of a differential gene are the most
  interpretable candidates.

## 9. What you can and can't conclude

- Hi-C averages millions of cells. A 20 % change in a loop can be a small change in every cell,
  or a large change in a few.
- Contact ≠ regulation. Many loops, when removed, change expression very little.
- Differential feature *calls* are fragile. Compare signal at a shared set of features.
- Depth, enzyme, protocol and normalisation all change what gets called, so harmonise them
  before comparing samples.

## 10. Extension: comparing species

Hi-C maps from different species live in different coordinates, and genome rearrangements
break simple coordinate mapping.

1. **Map each species on its own genome** (one Juicer run per species).
2. **Lift loop anchors** to the other genome separately. Use a strict match for short anchors,
   then require that (a) both anchors map, (b) to the same chromosome, (c) with a span within a
   tolerance (e.g. ± 25 %) of the original. Loops passing this are **cross-mappable**.
3. A cross-mappable loop is **conserved** only if the other species has its own loop call there,
   or significant APA / pixel signal at the lifted coordinates. Coordinate mapping alone says
   nothing about structure.
4. Compare TAD boundaries the same way: lift 500 bp around each boundary, require most of it to
   map, and allow a tolerance of a bin or two.
5. Remember known karyotype differences: chimpanzee chromosomes 2A and 2B correspond to human
   chromosome 2, for example.

```bash
# split BEDPE into anchors, lift each, then re-pair by loop ID in R
awk 'BEGIN{OFS="\t"} NR>1{print $1,$2,$3,"loop"NR"_L"; print $4,$5,$6,"loop"NR"_R"}' loops.bedpe > anchors.bed
liftOver -minMatch=0.9 anchors.bed hg38ToPanTro6.over.chain.gz anchors_panTro6.bed anchors_unmapped.bed
```

```r
read_anchors <- function(f) {
  a <- read.table(f, col.names = c("chr", "start", "end", "id"))
  a$loop <- sub("_[LR]$", "", a$id); a$side <- sub(".*_", "", a$id)
  l <- merge(a[a$side == "L", ], a[a$side == "R", ], by = "loop", suffixes = c("_L", "_R"))
  l$span <- abs(l$start_R - l$start_L); l
}
orig   <- read_anchors("anchors.bed")               # source genome
lifted <- read_anchors("anchors_panTro6.bed")       # target genome (loops with both anchors lifted)
lifted$orig_span <- orig$span[match(lifted$loop, orig$loop)]
cross_mappable <- subset(lifted, chr_L == chr_R & abs(span - orig_span) <= 0.25 * orig_span)
nrow(cross_mappable) / nrow(orig)                   # fraction of loops that are cross-mappable
```

::: {.callout-note}
## Example of results
Ciuba *et al.* (2025) generated in situ Hi-C from human primary fetal astrocytes and from human
and chimpanzee iPSC-derived astrocytes (ArrayExpress E-MTAB-13259). They called TopDom domains at
25 kb. TAD architecture was largely preserved between species, and species-specific boundaries
mostly reflected stronger insulation at positions with some boundary potential in the other
species, rather than new boundaries. Another public comparative dataset: human and chimpanzee
iPSC Hi-C from [Eres *et al.* 2019](https://doi.org/10.1371/journal.pgen.1008278).
:::

## 11. Reproducibility checklist and exercises

- [ ] Genome FASTA version, restriction enzyme and Juicer commit recorded; `snakemake --report` saved.
- [ ] Library statistics (`results/hic/<sample>.inter_30.txt`: valid pairs, cis/trans, long-range share) reviewed per sample.
- [ ] Resolution, normalisation (e.g. KR) and caller parameters (TopDom window, HiCCUPS resolutions) saved with every domain and loop set.
- [ ] `set.seed()` called before random background or permutation tests; `sessionInfo()` saved.
- [ ] Samples compared at matched depth, or the depth difference noted.

**Exercises**

1. Your library's intra-chromosomal contacts are 90 % shorter than 20 kb. What went wrong?
2. Why can't you compare two samples' maps using KR-normalised values alone, if one has twice
   the depth and a different P(s) slope?
3. A human loop lifts over to chimpanzee but no loop is called there. List three explanations.

::: {.callout-note collapse="true"}
## Answers
1. Mostly unligated or undigested DNA, and self-ligation products (dangling ends). Check the
   digestion efficiency and fill-in/ligation steps.
2. KR equalises total visibility per locus, not depth or distance dependence. Compare
   observed/expected, or depth-matched values after checking P(s).
3. The loop is genuinely absent. It's present but below calling sensitivity (depth, caller).
   Or the lifted anchors land in an unmappable or rearranged region. Test the signal with APA or
   pixel counts before deciding.
:::

**Further reading:** [Juicer wiki](https://github.com/aidenlab/juicer/wiki),
[Open2C tutorials](https://open2c.github.io/), [plotgardener vignettes](https://phanstiellab.github.io/plotgardener/).
