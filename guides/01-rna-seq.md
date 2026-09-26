---
title: "RNA-seq: from raw reads to differentially expressed genes"
subtitle: "Bulk RNA-seq with STAR, featureCounts and DESeq2"
---

::: {.callout-tip}
## What this guide covers
- Designing a bulk RNA-seq experiment with enough biological replication to answer the question.
- Processing paired-end FASTQ files into a gene count matrix, with checks at each step.
- Exploring sample structure (PCA, clustering) to spot outliers and batch effects.
- Testing for differential expression with DESeq2: design, contrasts and fold-change shrinkage.
- Visualising and functionally interpreting the results, and what they can and can't show.
- Extending the analysis to a comparison across species.

**Time:** about 3 hours. **You need:** a terminal, R ≥ 4.4 with Bioconductor, and the set-up
from [Getting started](00-getting-started.md).
:::

## 1. What RNA-seq measures

Bulk RNA-seq sequences a random sample of the (usually poly-A-selected) RNA in a population of
cells. After alignment and counting, each gene has a **count**: the number of sequenced fragments
assigned to it. Counts are proportional to expression *and* to gene length, sequencing depth and
library composition, so they are only comparable between samples after normalisation, and
between genes only with length correction (TPM).

**It can tell you:** which genes change in average expression between conditions, by how much,
and with what confidence.
**It can't tell you:** which cells changed (a shift in cell-type proportions looks like
differential expression; see the [scRNA-seq guide](05-scrna-seq.md)), whether protein levels
changed, or *why* a gene changed.

## 2. Experimental design

Decisions made here matter more than any later analysis choice.

| Factor | Recommendation | Why |
|---|---|---|
| Biological replicates | ≥ 3 per group; 4–6 for subtle effects | Variance between replicates is what the statistics estimate. Adding replicates buys much more power than adding depth ([Liu *et al.* 2014](https://doi.org/10.1093/bioinformatics/btt688); [Schurch *et al.* 2016](https://doi.org/10.1261/rna.053959.115)) |
| Technical replicates | Rarely needed | Re-sequencing a library adds little; sum its counts with the original |
| Depth | 20–30 M read pairs per sample for gene-level DE | Enough for most expressed genes; go deeper for isoforms or low-expressed genes |
| Read layout | Paired-end, stranded | Stranded libraries resolve overlapping antisense genes |
| Batches | Balance conditions across batches (RNA extraction day, library prep, lane) | A batch that coincides with a condition can't be separated from it |
| Covariates | Record sex, age, passage, RIN… | Anything you record you can model; anything you don't becomes noise |

::: {.callout-warning}
## Confounding can't be fixed afterwards
If every treated sample was prepared on Monday and every control on Tuesday, no model can tell
the treatment effect from the day effect. Randomise or balance *before* the experiment.
:::

## 3. Quality control of raw reads

The pipeline runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on every
FASTQ before and after trimming, and gathers everything into one
[MultiQC](https://multiqc.info) report (`results/multiqc/multiqc_report.html`). How to read the
main FastQC modules for RNA-seq:

| Module | Normal for RNA-seq | Worry if |
|---|---|---|
| Per-base quality | Median Phred > 28 across most of the read, dropping slightly at the 3′ end | Large drops mid-read, or a whole run is low |
| Per-base sequence content | Uneven first ~10 bases (random-hexamer priming bias) | Unevenness along the whole read |
| Per-sequence GC | Roughly one broad peak | A second peak (contamination, e.g. bacteria or rRNA) |
| Duplication | **High** (30–60 %) is expected: highly expressed genes give many identical reads | Nearly everything duplicated (very low input) |
| Overrepresented sequences / adapters | A few highly expressed transcripts; adapter at the 3′ end of short inserts | rRNA dominating, or adapter in most reads |

FastQC's red "fail" flags are calibrated for random genomic DNA. For RNA-seq, "sequence
duplication" and "per-base sequence content" routinely fail on perfectly good data.

## 4. From FASTQ to a count matrix

```bash
cd pipelines/rnaseq
$EDITOR config.yaml samples.tsv       # genome FASTA + GTF, your samples
snakemake --use-conda --cores 16 -n   # dry run
snakemake --use-conda --cores 16
```

Reference files: use a genome FASTA and a GTF **from the same source and release**, e.g. Ensembl
`Homo_sapiens.GRCh38.dna.primary_assembly.fa` + `Homo_sapiens.GRCh38.<release>.gtf`, or GENCODE's
`GRCh38.primary_assembly.genome.fa` + `gencode.v<N>.primary_assembly.annotation.gtf`.

What each step does, and the main alternatives:

**4.1 Trimming (Trim Galore).** Removes adapter read-through and low-quality 3′ ends. STAR
soft-clips unaligned read ends anyway, so for alignment-based counting trimming changes little.
It matters more for short inserts and for pseudo-aligners.

**4.2 Alignment (STAR).** A splice-aware aligner: RNA-seq reads span exon–exon junctions, so a
DNA aligner (bowtie2, BWA) would fail on them. The genome index is built once, with the GTF so
STAR knows the annotated junctions.

**4.3 Strandedness check.** featureCounts must be told the library's strand orientation (`-s`),
and a wrong setting silently loses or mis-assigns most reads. The pipeline infers it from STAR's
`ReadsPerGene.out.tab`, which lists counts for both orientations:

```bash
# column 3 = counts if read 1 is on the gene's strand (featureCounts -s 1),
# column 4 = counts if read 2 is (featureCounts -s 2). The larger one wins.
awk 'NR>4 {s1+=$3; s2+=$4} END {print "s1:", s1, " s2:", s2, " ratio:", s1/s2}' \
    results/star/<sample>.<genome>.ReadsPerGene.out.tab
cat results/qc/strandedness.txt       # the pipeline's call: 0, 1 or 2
```

Most current stranded kits (dUTP-based: Illumina Stranded mRNA, NEBNext Ultra II Directional,
KAPA mRNA HyperPrep) are **reverse-stranded, `-s 2`**. A ratio near 1 means unstranded (`-s 0`).

**4.4 Counting (featureCounts).** Counts each read pair once for the gene whose exons it
overlaps. By default, pairs overlapping two genes and multi-mapping pairs are left uncounted.
That is conservative but loses reads from paralogous families and overlapping genes.

::: {.callout-note collapse="true"}
## Alternative: alignment-free quantification (Salmon, kallisto)
Tools like [Salmon](https://combine-lab.github.io/salmon/) quantify transcripts directly
against the transcriptome, typically 10–20× faster than STAR. They model multi-mapping reads
across isoforms statistically rather than discarding them. Import the result into DESeq2 with
[`tximport`](https://bioconductor.org/packages/tximport/). HBC's
[Intro-to-bulk-RNAseq](https://hbctraining.github.io/Intro-to-bulk-RNAseq/) course uses this
route. **Pros:** speed, better isoform handling, gene-length bias correction. **Cons:** no
genome BAM for browsing or QC, it misses reads from unannotated genes, and it's harder to use with
custom or masked genomes such as the cross-species reference in the appendix. Both routes give
very similar gene-level DE results for well-annotated organisms.
:::

**Checkpoints after processing** (all in the MultiQC report):

| Metric | Typical good value | If not |
|---|---|---|
| STAR uniquely mapped | > 70–80 % (human/mouse) | Contamination, wrong genome, or degraded RNA |
| STAR multi-mapped | < 10–20 % | Very high: rRNA contamination |
| featureCounts assigned | > 60 % | Wrong strandedness, wrong GTF, or many intronic reads (nuclear RNA, degradation) |
| Strandedness ratio | > 4 or < 0.25 for stranded kits | Near 1: the library is unstranded |
| Genes detected | Similar across samples | An outlier sample may have failed |

## 5. Loading counts into R

The code from here on runs as written on the **airway** dataset (Himes *et al.* 2014): airway
smooth-muscle cells from four donors, treated or not with dexamethasone. It ships with
Bioconductor, so you can practise before using your own data.

```r
# BiocManager::install(c("DESeq2", "airway", "apeglm", "pheatmap", "clusterProfiler", "org.Hs.eg.db"))
library(DESeq2)
library(airway)
data("airway")
counts <- assay(airway)                                   # genes x samples, raw integer counts
coldata <- as.data.frame(colData(airway))[, c("cell", "dex")]
coldata$dex <- relevel(coldata$dex, ref = "untrt")        # the reference level = the denominator
```

With **your own pipeline output**, read the featureCounts table instead:

```r
fc <- read.delim("results/featureCounts/GRCh38_gene_counts.txt", comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, -(1:6)])                          # columns 1-6: Geneid, Chr, Start, End, Strand, Length
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\..*Aligned.sortedByCoord.out.bam$", "", basename(colnames(counts)))
coldata <- read.delim("sample_info.tsv", row.names = 1)    # one row per sample: condition, batch, ...
stopifnot(all(colnames(counts) == rownames(coldata)))      # same samples, same order
```

::: {.callout-important}
## Keep raw counts raw
DESeq2 needs **unnormalised integer counts**. Never give it TPM, FPKM or already-normalised values.
Normalise only copies used for plotting.
:::

## 6. Sample-level QC: do the replicates agree?

```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ cell + dex)
smallestGroupSize <- 4                                     # samples in the smallest group
dds <- dds[rowSums(counts(dds) >= 10) >= smallestGroupSize, ]   # drop near-empty genes

vsd <- vst(dds, blind = TRUE)            # variance-stabilised log2-like scale, for plots only
plotPCA(vsd, intgroup = c("dex", "cell"))

library(pheatmap)
d <- as.matrix(dist(t(assay(vsd))))
pheatmap(d, annotation_col = coldata, clustering_distance_rows = "euclidean")
```

How to read these plots:

- **PCA:** samples should separate by the biology you expect, and replicates should sit
  together. In airway, PC1 separates treated from untreated and PC2 separates donors, so donor
  belongs in the design.
- **An outlier** far from its group: check its QC metrics before removing it. Only remove a
  sample for a documented technical reason, never because it weakens the result.
- **Batch:** if samples group by processing date rather than condition, add batch to the
  design (next section). Don't "correct" counts before DESeq2.

## 7. Differential expression with DESeq2

### 7.1 The design formula

The design says which variables explain variation in expression. Put the variable of interest
**last**; everything before it is controlled for.

| Situation | Design |
|---|---|
| Two groups | `~ condition` |
| Two groups + known batch / paired donors | `~ batch + condition` (airway: `~ cell + dex`) |
| Does the treatment effect differ between genotypes? | `~ genotype + treatment + genotype:treatment` |
| Several levels / time course | Wald contrasts between levels, or the likelihood-ratio test (below) |

### 7.2 Fit and extract results

```r
dds <- DESeq(dds)                                    # size factors, dispersions, GLM fit
resultsNames(dds)                                    # names of the fitted coefficients
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
summary(res)

# Shrunken log2 fold changes: better for ranking and plotting (p-values don't change)
res_shr <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
sig <- subset(res_shr, padj < 0.05 & abs(log2FoldChange) > 1)
sig <- sig[order(sig$padj), ]
head(sig)
```

What the steps do:

1. **Size factors** (median-of-ratios) correct for depth and composition. A few very highly
   expressed genes can't distort them the way total-count scaling can.
2. **Dispersion** is the gene-wise variance beyond Poisson noise. With few replicates each gene's
   estimate is noisy, so DESeq2 shrinks it towards a trend fitted across all genes.
3. **Wald test** per gene for the chosen contrast, then Benjamini–Hochberg adjustment. `padj`
   is the false discovery rate (FDR), and `alpha` should match the FDR cut-off you'll use.
4. **Automatic filters** (`NA` padj): genes with an extreme single-sample outlier (Cook's
   distance), or too few counts to ever be significant (independent filtering).

::: {.callout-note collapse="true"}
## More than two groups, time courses: the likelihood-ratio test
To ask "does this gene change *at all* across any of several conditions or time points", compare
a full model against a reduced one:

```r
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ cell)   # full design = ~ cell + dex
res_lrt <- results(dds_lrt)
```
The LRT p-value tests every level at once. The fold change it reports is only for one
coefficient, so use Wald contrasts to describe *which* levels differ.
:::

### 7.3 Pros and cons of the main choices

| Choice | Pros | Cons / when not |
|---|---|---|
| DESeq2 (negative binomial GLM) | Robust with 3–10 replicates; handles complex designs; widely used | Assumes most genes don't change. Very large global shifts (e.g. transcriptional amplification) need spike-ins |
| edgeR / limma-voom | Equally well validated; limma-voom is fast and flexible with many samples and random effects | Different defaults, so results differ slightly in edge cases |
| LFC shrinkage (apeglm) | Stops noisy low-count genes topping the ranked list | Shrunken values are estimates for ranking; report the test's padj |
| Fold-change threshold after testing | Simple | Filtering on \|LFC\| after the test isn't the same as testing against it; use `results(..., lfcThreshold = 1)` for that |
| Pre-filtering low counts | Faster, smaller multiple-testing burden | Only an efficiency step: independent filtering already handles power |

## 8. Visualising results

```r
plotMA(res_shr, ylim = c(-5, 5))                         # LFC vs mean expression

library(ggplot2)
df <- as.data.frame(res_shr)
df$sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
ggplot(df, aes(log2FoldChange, -log10(pvalue), colour = sig)) +          # volcano plot
  geom_point(size = 0.6) + scale_colour_manual(values = c("grey70", "firebrick")) +
  theme_bw()

top <- head(rownames(sig), 30)                                            # heatmap of top genes
mat <- assay(vst(dds, blind = FALSE))[top, ]
pheatmap(mat - rowMeans(mat), annotation_col = coldata)

plotCounts(dds, gene = rownames(sig)[1], intgroup = "dex")               # one gene, all samples
```

| Plot | Question it answers | Red flag |
|---|---|---|
| MA plot | Are changes spread across expression levels? | A cloud skewed up or down: normalisation problem |
| Volcano | How many genes change, how strongly, how confidently? | Everything significant with tiny fold changes: huge n or a confound |
| Heatmap of top genes | Do all replicates agree? | One sample driving the pattern |
| p-value histogram (`hist(res$pvalue)`) | Is the test well calibrated? | Not flat with a peak near 0: hidden batch or a mis-specified design |

## 9. Functional interpretation

**Over-representation analysis (ORA)** asks whether a gene set (e.g. a GO term) appears more
often among your significant genes than chance predicts. The *universe* must be all genes that
could have been detected, not the whole genome.

```r
library(clusterProfiler); library(org.Hs.eg.db)
universe <- rownames(res)[!is.na(res$padj)]                   # tested genes only
up <- rownames(subset(res_shr, padj < 0.05 & log2FoldChange > 1))
ego <- enrichGO(up, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                universe = universe, pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(ego, showCategory = 15)
```

**Gene set enrichment analysis (GSEA)** uses *all* genes ranked by a statistic, so it has no
significance threshold and can detect coordinated small shifts:

```r
ranks <- res$stat[!is.na(res$stat)]; names(ranks) <- rownames(res)[!is.na(res$stat)]
ranks <- sort(ranks, decreasing = TRUE)
gse <- gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", seed = TRUE)   # a warning about tied ranks is harmless
```

Caveats: enrichment tells you which *annotated* processes are over-represented, and annotation
is biased towards well-studied genes. Treat GO results as hypotheses, and prefer a few specific,
coherent terms over long lists of generic ones.

## 10. What you can and can't conclude

- **Significant ≠ large.** With many replicates, tiny changes become significant. Look at fold
  changes and expression levels, not just padj.
- **Not significant ≠ unchanged.** Absence of evidence at your sample size isn't evidence of
  equal expression. To claim "no change", test for equivalence
  (`results(dds, lfcThreshold = 0.5, altHypothesis = "lessAbs")`).
- **Composition effects.** In tissue, a gene "going up" may mean the cells expressing it became
  more abundant. Check marker genes of the main cell types, or deconvolve.
- **n = 1 per group** gives no estimate of biological variance. Results are descriptive at best.
- **Validate the story.** Confirm key genes by an independent method (qPCR, protein, a second
  dataset) and check they're expressed at a level that could matter.

## 11. Extension: comparing species

Cross-species RNA-seq needs two extra steps. Genes must be put in a shared ID space, and a
change must be assigned to the right lineage.

**Shared gene space: two options.**

1. **Map each species to its own genome, then keep 1:1 orthologs.** Count each species on its own
   genome and GTF, then subset every matrix to genes with one-to-one orthologs (Ensembl
   Compara via `biomaRt`) and rename rows to one species' IDs. This is simple and standard, but
   loses genes without clean 1:1 orthologs, and annotation quality differs between species.
2. **Map every species to one masked reference.** For closely related species (e.g. great apes),
   build a *consensus genome*: the human assembly with every base that differs from the other
   species masked to N. Reads from all species then align with no reference bias and are counted
   with one GTF ([appendix](appendix-comparative-genomics.md);
   [`pipelines/consensus_genome`](../pipelines/README.md)). This needs a whole-genome alignment
   and only works for species close enough to align (roughly within the primates).

**Which lineage changed? Use an outgroup.** A human–chimpanzee difference alone doesn't say
which species changed. With a third, more distant species (e.g. macaque), a gene is a candidate
**human-lineage** change only if human differs from *both* others in the same direction:

```r
# dds_hc: human vs chimp; dds_hm: human vs macaque (separate fits, each with ~ species)
r_hc <- results(dds_hc, contrast = c("species", "human", "chimp"))
r_hm <- results(dds_hm, contrast = c("species", "human", "macaque"))
g <- intersect(rownames(r_hc), rownames(r_hm))
both <- data.frame(lfc_hc = r_hc[g, "log2FoldChange"], padj_hc = r_hc[g, "padj"],
                   lfc_hm = r_hm[g, "log2FoldChange"], padj_hm = r_hm[g, "padj"], row.names = g)
human_up   <- subset(both, padj_hc < 0.01 & padj_hm < 0.01 & lfc_hc > 0 & lfc_hm > 0)
human_down <- subset(both, padj_hc < 0.01 & padj_hm < 0.01 & lfc_hc < 0 & lfc_hm < 0)
plot(both$lfc_hc, both$lfc_hm, pch = ".", xlab = "log2FC human vs chimp", ylab = "log2FC human vs macaque")
```

Extra cautions for cross-species work:

- **Cell models must be comparable.** With iPSC-derived cells, differentiation efficiency can
  differ between lines and species. Score each sample on independent maturity or identity marker
  sets before trusting a species difference, and drop genes that are themselves maturity markers.
- **Individuals, not cells, are replicates.** Several cultures from one individual don't count
  as biological replicates of that species.

::: {.callout-note}
## Example of results
Ciuba *et al.* (2025, *Cell Stem Cell*) compared iPSC-derived astrocytes from human, chimpanzee
and macaque this way (consensus-genome alignment, pairwise DESeq2, congruence across both
outgroup comparisons). They report 677 up- and 486 down-regulated congruent genes, which were
then filtered for expression in primary human astrocytes. Their raw data are in ArrayExpress
(E-MTAB-13252); [Getting started](00-getting-started.md#example-data) shows how to download them.
:::

## 12. Optional: alternative splicing

Gene-level counts miss changes in *which* isoforms are used. [rMATS](https://github.com/Xinglab/rmats-turbo)
detects differential splicing events (skipped exons, retained introns, alternative 5′/3′ sites,
mutually exclusive exons) directly from the STAR BAMs this pipeline produces:

```bash
rmats.py --b1 group1_bams.txt --b2 group2_bams.txt --gtf genes.gtf -t paired \
         --readLength 100 --libType fr-firststrand --nthread 8 --od rmats_out --tmp rmats_tmp
```

`--libType fr-firststrand` corresponds to featureCounts `-s 2` (dUTP kits), and `fr-secondstrand`
to `-s 1`. Alternatives: [DEXSeq](https://bioconductor.org/packages/DEXSeq/) (exon usage) or
isoform-level quantification with Salmon + [DRIMSeq](https://bioconductor.org/packages/DRIMSeq/).

## 13. Reproducibility checklist and exercises

- [ ] Genome FASTA and GTF versions recorded; `snakemake --report` saved.
- [ ] Strandedness taken from `results/qc/strandedness.txt`, not assumed.
- [ ] Sample sheet and metadata in version control; sample order checked with `stopifnot()`.
- [ ] `sessionInfo()` saved with the results.
- [ ] Every excluded sample and its technical reason written down.

**Exercises** (airway data)

1. Fit `~ dex` without `cell`. How many genes pass padj < 0.05, compared with `~ cell + dex`? Why?
2. Make the p-value histogram for both models. Which looks better calibrated?
3. Rank genes by `res$log2FoldChange` and then by `res_shr$log2FoldChange`. Which genes fall out
   of the top 20 after shrinkage, and what do they have in common?

::: {.callout-note collapse="true"}
## Answers
1. The donor-aware design finds more genes. Donor (`cell`) differences are large, and
   modelling them removes that variation from the residual noise, so the treatment effect is
   estimated more precisely (a paired design).
2. Both are flat with a spike near 0. The donor-aware model's spike is taller: more real signal
   is detected.
3. Low-count genes with huge but unreliable fold changes drop out after shrinkage. That is what
   shrinkage is for.
:::

**Further reading:** the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
the [RNA-seq workflow](https://bioconductor.org/packages/rnaseqGene/), and HBC's
[Intro-to-DGE](https://hbctraining.github.io/Intro-to-DGE/) course.
