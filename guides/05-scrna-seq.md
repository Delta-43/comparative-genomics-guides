---
title: "scRNA-seq: cell types, states and their differences"
subtitle: "From 10x FASTQ to annotated clusters, pseudobulk DE and abundance testing"
---

::: {.callout-tip}
## Learning objectives
- Design a single-cell experiment: cells vs samples, multiplexing, batches.
- Generate count matrices with Cell Ranger, including antibody-capture (CITE-seq / hashtag) data.
- Run QC, doublet removal, normalisation, dimensionality reduction and clustering in R, and
  judge clusterings with numbers, not just UMAP pictures.
- Annotate cell types; test expression differences with pseudobulk and abundance differences
  with compositional methods.
- Extend the design and analysis across species.

**Time:** about 4 hours. **You need:** R ≥ 4.4 with Bioconductor, the set-up from
[Getting started](00-getting-started.md), and ~16 GB RAM for the example data.
:::

## 1. What scRNA-seq measures

Droplet methods (10x Genomics Chromium) capture single cells with a barcoded bead. Every
transcript gets a cell barcode and a unique molecular identifier (UMI), so counts are
**molecules per gene per cell**. Only a fraction (~10–30 %) of each cell's mRNA is captured,
which makes the data sparse and zeros common.

**It can tell you:** which cell types and states are present, their relative abundance, and how
expression differs within a cell type. **It can't tell you:** absolute cell numbers in the
tissue (dissociation favours some cell types), spatial arrangement, or reliable information about
lowly expressed genes in single cells.

## 2. Experimental design

| Factor | Recommendation | Why |
|---|---|---|
| Biological replicates | ≥ 3 samples (donors/animals) per condition | Samples, not cells, are the unit of replication for any comparison between conditions ([Squair *et al.* 2021](https://doi.org/10.1038/s41467-021-25960-2)) |
| Cells per sample | 3,000–10,000 | Enough to see populations at ~1 % |
| Depth | ≥ 20,000 reads per cell for 3′ gene expression | Check sequencing saturation in the Cell Ranger report |
| Multiplexing | Pool samples in one run with cell hashing ([Stoeckius *et al.* 2018](https://doi.org/10.1186/s13059-018-1603-1)) or genotype-based demultiplexing | Removes run-to-run batch effects between samples, identifies cross-sample doublets, and cuts cost |
| Batches | If samples must be run separately, spread conditions across runs | A run that coincides with a condition can't be corrected away |

## 3. From FASTQ to a count matrix (Cell Ranger)

```bash
cd pipelines/scrna
$EDITOR config.yaml runs.tsv templates/library.csv    # reference, runs, FASTQ folders
snakemake --cores 32 -n && snakemake --cores 32
```

Use 10x's prebuilt reference for human or mouse (`refdata-gex-GRCh38-2024-A`), or build one
with `cellranger mkref` (the pipeline does this if no prebuilt reference is set). The core
command:

```bash
cellranger count --id=<run> --libraries=library.csv --feature-ref=feature_ref.csv \
    --transcriptome=<reference> --create-bam=true --localcores=16 --localmem=64
```

- `library.csv` lists each FASTQ folder and its `library_type` (`Gene Expression`,
  `Antibody Capture`). A library sequenced on several runs is simply listed once per run folder.
- `feature_ref.csv` lists antibody barcodes (hashtags and/or surface-protein antibodies). The
  template holds BioLegend TotalSeq-B human hashtags 1–4 as listed by 10x.
- Separate libraries from different samples are counted separately and, if needed, combined in R
  (or with `cellranger aggr`).

**Checkpoints** (`outs/web_summary.html`): estimated cells close to what you loaded;
reads mapped confidently to the transcriptome > 60–70 %; fraction of reads in cells > 70 %; a
barcode-rank "knee" plot with a clear drop between cells and empty droplets. 10x's report flags
values outside its expected ranges.

## 4. Loading data and quality control

The code below runs as written on the 10x **PBMC 3k** dataset (from Bioconductor's
`TENxPBMCData`), a standard teaching dataset. Load your own Cell Ranger output with
`DropletUtils::read10xCounts()` instead (shown in the comments).

```r
library(TENxPBMCData); library(scater); library(scran); library(scuttle)
sce <- TENxPBMCData("pbmc3k")
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL_ID, rowData(sce)$Symbol_TENx)
counts(sce) <- as(counts(sce), "dgCMatrix")                     # in-memory sparse matrix
# Your data:
# sce <- DropletUtils::read10xCounts("results/<run>/outs/filtered_feature_bc_matrix")
# rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
# sce <- splitAltExps(sce, rowData(sce)$Type)                   # antibody capture -> altExp(sce)

is_mito <- grepl("^MT-", rownames(sce))                          # human; mouse: "^mt-"
qc <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))
filt <- perCellQCFilters(qc, sub.fields = "subsets_Mito_percent")  # outliers: 3 MADs from the median
colSums(as.matrix(filt))
colData(sce) <- cbind(colData(sce), qc)
sce$discard <- filt$discard
plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "discard") + scale_x_log10()
sce <- sce[, !sce$discard]
```

**Why adaptive (MAD-based) thresholds?** Fixed cut-offs (e.g. "> 500 UMIs, < 10 % mito") are
simple, but what counts as "low" differs by tissue, chemistry and depth. Outlier thresholds adapt
to each dataset. Do check the plots: if one real cell type (e.g. metabolically active
cardiomyocytes, or small lymphocytes) is removed wholesale, relax the filter for that dataset.
Always run QC **per sample**, since samples differ in depth.

### Doublets

Two cells in one droplet give a hybrid profile that can look like a novel intermediate cell
type. [scDblFinder](https://doi.org/10.12688/f1000research.73600.2) simulates doublets and scores
real cells against them:

```r
library(scDblFinder)
set.seed(100)
sce <- scDblFinder(sce)                    # with several samples: scDblFinder(sce, samples = "sample")
table(sce$scDblFinder.class)
sce <- sce[, sce$scDblFinder.class == "singlet"]
```

With cell hashing, cross-sample doublets are also identified directly from the hashtags (§9).

## 5. Normalisation, feature selection, dimensionality reduction

```r
set.seed(101)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)   # pooling/deconvolution size factors
sce <- logNormCounts(sce)

dec <- modelGeneVar(sce)                             # mean-variance trend
hvg <- getTopHVGs(dec, prop = 0.1)                   # top 10% most variable genes
set.seed(102)
sce <- fixedPCA(sce, subset.row = hvg, rank = 30)
sce <- runUMAP(sce, dimred = "PCA")
```

- **Deconvolution size factors** ([Lun *et al.* 2016](https://doi.org/10.1186/s13059-016-0947-7))
  pool cells to get around sparse counts, then deconvolve per-cell factors. They are robust to
  composition differences between cell types.
- **Highly variable genes** carry the biological signal. Restricting PCA to them removes noise.
- **PCA first, UMAP only for display.** Cluster on PCs. UMAP distances between separated
  clusters don't mean anything quantitative, and the picture changes with `n_neighbors` and the
  seed.

## 6. Clustering, and judging a clustering

```r
library(bluster)
set.seed(103)
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = SNNGraphParam(k = 10, type = "rank", cluster.fun = "walktrap"))
table(colLabels(sce))
plotReducedDim(sce, "UMAP", colour_by = "label")

sil <- approxSilhouette(reducedDim(sce, "PCA"), colLabels(sce))   # > 0: closer to own cluster
pur <- neighborPurity(reducedDim(sce, "PCA"), colLabels(sce))     # fraction of neighbours in own cluster
data.frame(cluster = levels(colLabels(sce)),
           silhouette = tapply(sil$width, colLabels(sce), mean),
           purity = tapply(pur$purity, colLabels(sce), mean))
```

There is no single "correct" number of clusters. `k` (neighbours) and the algorithm (Walktrap,
Louvain, Leiden) set the resolution. Choose the resolution that answers your question. A cluster
with low silhouette and purity is probably a continuum or a mix; check its markers before
calling it a cell type.

## 7. Marker genes and cell-type annotation

```r
markers <- scoreMarkers(sce, colLabels(sce))
top <- lapply(markers, function(m) head(rownames(m)[order(m$mean.AUC, decreasing = TRUE)], 10))
top[["1"]]
plotExpression(sce, features = c("CD3E", "MS4A1", "LYZ", "NKG7", "PPBP"), x = "label")

library(SingleR); library(celldex)
ref <- celldex::BlueprintEncodeData()
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
table(cluster = colLabels(sce), SingleR = pred$labels)
plotScoreHeatmap(pred)
```

| Approach | Pros | Cons |
|---|---|---|
| Marker genes (manual) | Uses expert knowledge; works for any tissue or species | Subjective; slow for many clusters |
| Reference-based (SingleR, Azimuth, CellTypist) | Fast, reproducible | Only as good as the reference. Novel or disease states get forced into the closest known label |
| Protein (CITE-seq) markers | Surface markers match flow-cytometry definitions | Limited panel; antibody background needs correcting (DSB) |

Use both: a reference for a first pass, then confirm with canonical markers.

## 8. Comparing conditions

### 8.1 Differential expression: pseudobulk

Treating each cell as a replicate inflates significance enormously: cells from one donor aren't
independent. Instead, **sum counts per cell type per sample** and use bulk methods
([Squair *et al.* 2021](https://doi.org/10.1038/s41467-021-25960-2); [muscat](https://doi.org/10.1038/s41467-020-19894-4)):

```r
# Needs sample and condition labels per cell, e.g. from hashtags or separate runs:
# sce$sample <- ...; sce$condition <- ...
library(DESeq2)
pb <- aggregateAcrossCells(sce, ids = colData(sce)[, c("label", "sample")], use.assay.type = "counts")
ct <- "1"                                                           # one cell type at a time
pb1 <- pb[, pb$label == ct & pb$ncells >= 10]                       # drop sample x cluster with few cells
meta <- as.data.frame(colData(pb1))[, c("sample", "condition", "ncells")]   # keep only what you need:
# the single-cell `sizeFactor` column would otherwise be picked up by DESeq2 and break it
dds <- DESeqDataSetFromMatrix(counts(pb1), colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"))
```

From here, everything in the [RNA-seq guide](01-rna-seq.md#7-differential-expression-with-deseq2)
applies: design formulas, shrinkage, plots, enrichment.

### 8.2 Differential abundance

Cell-type proportions sum to one, so they aren't independent. Test them with methods built for
proportions across replicate samples, such as `speckle::propeller`
([Phipson *et al.* 2022](https://doi.org/10.1093/bioinformatics/btac582)), or with neighbourhood
methods that don't need discrete clusters (miloR):

```r
library(speckle)
propeller(clusters = colLabels(sce), sample = sce$sample, group = sce$condition)
```

A change in proportion can come from biology **or** from dissociation, sorting or capture
differences between samples. Be careful with fragile cell types.

::: {.callout-note collapse="true"}
## Practise on a real multi-sample dataset
The PBMC 3k data have one sample, so they can't show between-condition statistics. For practice,
use `muscData::Kang18_8vs8()` (Bioconductor): PBMCs from 8 lupus patients, each unstimulated
and IFN-β-stimulated, multiplexed and demultiplexed by genotype. It is the standard example for
pseudobulk DE and abundance testing (used in the muscat vignette).
:::

## 9. Multiplexed samples: demultiplexing hashtags

With cell hashing, each sample's cells carry a different hashtag antibody, counted in the
antibody-capture library (`altExp(sce)` after `splitAltExps`):

```r
library(DropletUtils)
hto <- counts(altExp(sce))                          # hashtags x cells
dm <- hashedDrops(hto)                              # best hashtag, log-fold change vs second best
table(dm$Confident, dm$Doublet)
sce <- sce[, which(dm$Confident)]                    # confident singlets only
key <- read.csv("templates/hashtag_to_sample.csv")  # id,name
sce$sample <- key$name[match(rownames(hto)[dm$Best[which(dm$Confident)]], key$id)]
```

A droplet with two strong hashtags is a cross-sample doublet. One with no clear hashtag is
ambiguous, usually a low-quality cell or an empty droplet.

## 10. Several batches: integration

If samples were run in separate batches, cells may cluster by batch rather than biology. Check
first by colouring the UMAP by batch. If needed, correct the **embedding** (never the counts
used for DE):

```r
library(batchelor)
set.seed(104)
merged <- fastMNN(sce, batch = sce$batch, subset.row = hvg)     # MNN correction (Haghverdi et al. 2018)
reducedDim(sce, "MNN") <- reducedDim(merged, "corrected")
sce <- runUMAP(sce, dimred = "MNN", name = "UMAP_MNN")
# Alternative: harmony::RunHarmony(sce, group.by.vars = "batch")  -> reducedDim(sce, "HARMONY")
```

Over-correction can erase real differences, for example between conditions processed in
different batches. After integration, check that known cell types still separate, and that
condition-specific populations weren't merged away. Pseudobulk DE (§8.1) should use raw counts
with batch in the design, not corrected values.

::: {.callout-note collapse="true"}
## Alternative: Seurat
[Seurat](https://satijalab.org/seurat/) implements the same steps with different defaults:
`NormalizeData`/`SCTransform` → `FindVariableFeatures` → `RunPCA` → `FindNeighbors`/`FindClusters`
(Louvain/Leiden) → `RunUMAP` → `FindAllMarkers`, with CCA/RPCA or Harmony for integration.
HBC's [Intro-to-scRNAseq](https://hbctraining.github.io/Intro-to-scRNAseq/) course teaches this
route. The statistical advice here (sample-level replication, pseudobulk, compositional
testing) applies equally. `as.Seurat()` / `as.SingleCellExperiment()` convert between the two.
:::

## 11. What you can and can't conclude

- **A cluster is not a cell type** until markers, references, and ideally an orthogonal method
  agree.
- **Cells aren't replicates.** Any condition comparison needs samples as replicates:
  pseudobulk DE, propeller/miloR for abundance.
- **Proportions are relative.** An "increase" in one type can be a decrease in others.
- **UMAP geometry is illustrative**, not quantitative.
- **Absence of a gene's expression** in a cell often just means it wasn't captured
  (sparsity). Reason at the cluster or pseudobulk level.

## 12. Extension: comparing species

1. **Shared reference.** Either map each species to its own genome and keep one-to-one orthologs,
   or map all species to one masked **consensus genome** (reads from every species align without
   reference bias; [appendix](appendix-comparative-genomics.md)). The consensus route needs
   closely related species and a whole-genome alignment.
2. **Pool species in one run.** Hashtagging cells from several species (or genotype-based
   demultiplexing) puts every species through the same capture and sequencing. Species
   differences can then no longer be run effects.
3. **Annotate consistently.** Confirm that marker genes behave as markers in *each* species
   before transferring labels. Some markers change expression pattern between species.
4. **Compare at the sample level.** Pseudobulk per cell type per individual, then use the
   outgroup logic from the [RNA-seq guide](01-rna-seq.md#11-extension-comparing-species) to assign
   changes to a lineage.
5. **Watch QC across species.** Mitochondrial gene names, annotation completeness and
   mappability differ, so check QC distributions separately per species.

::: {.callout-note}
## Example datasets
A large adult primate dataset: Ma *et al.* 2022 ([Science](https://doi.org/10.1126/science.abo7257)),
single-nucleus RNA-seq of the dorsolateral prefrontal cortex in human, chimpanzee, macaque and
marmoset. Ciuba *et al.* (2025) used pseudobulk comparisons of astrocytes from published human and
macaque fetal cortex data (processed with Cell Ranger against a consensus genome) to check
genes found in iPSC-derived astrocytes.
:::

## 13. Exercises

1. Re-cluster PBMC 3k with `k = 5` and `k = 50`. How does the number of clusters change, and what
   happens to the mean silhouette width?
2. SingleR labels one cluster as a mixture of two cell types. What would you check?
3. You have 2 donors per species and 5,000 cells each. How many replicates does a pseudobulk
   comparison between species have?

::: {.callout-note collapse="true"}
## Answers
1. Smaller `k` gives more, smaller clusters; larger `k` merges them. Mean silhouette usually
   drops when a real population is split. Pick the resolution that matches the question and
   shows coherent markers.
2. The markers of both labels in that cluster (sub-cluster it), its doublet scores, and whether
   it is a continuum (e.g. differentiation) rather than two types.
3. Two per species. That's the minimum for any variance estimate, and gives little power.
   Report such results as preliminary.
:::

**Further reading:** [Orchestrating Single-Cell Analysis with Bioconductor (OSCA)](https://bioconductor.org/books/release/OSCA/)
([Amezquita *et al.* 2020](https://doi.org/10.1038/s41592-019-0654-x)), best-practice reviews
([Luecken & Theis 2019](https://doi.org/10.15252/msb.20188746)), HBC's
[Intro-to-scRNAseq](https://hbctraining.github.io/Intro-to-scRNAseq/) and
[Pseudobulk-for-scRNAseq](https://github.com/hbctraining/Pseudobulk-for-scRNAseq).
