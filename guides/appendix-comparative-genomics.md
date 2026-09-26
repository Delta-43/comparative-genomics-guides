---
title: "Appendix: comparing genomics data across species"
subtitle: "Coordinates, orthologs, consensus genomes and lineage assignment"
---

::: {.callout-tip}
## What this appendix covers
Every guide's cross-species extension rests on the same three tools. This page explains each
one once:

1. **liftOver**: translating coordinates between genome assemblies.
2. **Orthologs**: putting genes from different species into one ID space.
3. **Consensus (masked) genomes**: aligning reads from several species to one reference
   without favouring any of them.

It ends with the **outgroup logic** that assigns a difference to a lineage.
:::

## 1. Why cross-species data need special handling

Reads from species B aligned to species A's genome map worse wherever the two genomes differ.
Mismatches lower alignment scores, indels break alignments, and species-specific sequence has
nowhere to go. This **reference bias** makes species B look lower in exactly the regions that
diverged, which are the regions an evolutionary study cares about. Every strategy below is a
way to avoid it.

| Strategy | Use for | Pros | Cons |
|---|---|---|---|
| Each species on its own genome, then compare **coordinates** via liftOver | Peaks, loops, domains (ATAC, ChIP, Hi-C) | Each species gets its best alignment | Only regions that lift over can be compared |
| Each species on its own genome, then compare **genes** via 1:1 orthologs | Gene expression (bulk, single cell) | Standard; works across distant species | Loses genes without clean 1:1 orthologs; annotation quality differs between species |
| All species on one **masked consensus genome** | Expression, binding, accessibility in closely related species | No reference bias; one annotation; directly comparable coordinates | Needs a multiple alignment; masked sites lose reads; only for species close enough to align well |

## 2. liftOver and chain files

A **chain file** describes blocks of aligned sequence between two assemblies. UCSC builds them by
aligning whole genomes (`lastz` → `axtChain` → `chainNet` → `netChainSubset`) and publishes them
for common pairs:

```bash
# https://hgdownload.soe.ucsc.edu/goldenPath/<from>/liftOver/<from>To<To>.over.chain.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToPanTro6.over.chain.gz
liftOver -minMatch=0.5 peaks_hg38.bed hg38ToPanTro6.over.chain.gz peaks_panTro6.bed unmapped.bed
```

- **`-minMatch`** is the fraction of an interval's bases that must map (default 0.95). Short
  intervals (peaks) tolerate looser values (0.5). For paired features (loop anchors) be stricter,
  and then check both ends land on the same chromosome at a similar distance.
- **Lift in both directions.** A region that maps A → B and back B → A to the same place is a
  reciprocal match, which is much more trustworthy than a one-way lift.
- **Unmapped ≠ absent.** Intervals that fail to lift usually overlap insertions, deletions,
  rearrangements or assembly gaps. Report them as their own class.
- In R, `rtracklayer::liftOver(gr, import.chain("hg38ToPanTro6.over.chain"))` does the same on
  `GRanges` (unzip the chain first). It returns a `GRangesList`, because one interval can split
  across chain blocks.

## 3. One-to-one orthologs

For gene-level comparisons across species aligned to their own genomes:

```r
library(biomaRt)
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")   # if Ensembl is slow: mirror = "useast"
orth <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                             "ptroglodytes_homolog_ensembl_gene", "ptroglodytes_homolog_orthology_type"),
              mart = mart)
one2one <- subset(orth, ptroglodytes_homolog_orthology_type == "ortholog_one2one")
# Then: keep rows of the chimp count matrix in one2one$ptroglodytes_homolog_ensembl_gene and
# rename them to the matching human ensembl_gene_id before combining the matrices.
```

Record the Ensembl release you queried (`listEnsemblArchives()`; `useEnsembl(..., version = 113)`
pins one), because orthology calls change between releases. Ensembl's orthology comes from gene
trees (Ensembl Compara). One-to-one orthologs are the safe
set. One-to-many and many-to-many relationships (recent duplications) need a separate treatment,
or should be excluded. Gene lengths and annotation completeness differ between species, so compare
within-gene fold changes (DESeq2), not TPMs across species.

## 4. Consensus (masked) genomes

For closely related species (e.g. human, chimpanzee, macaque), a single reference can be built
that none of them is favoured by:

1. Take one assembly as the coordinate system (e.g. hg38).
2. Build a multiple alignment with the other species: UCSC pairwise alignments (chain → net →
   axt → MAF), combined with `multiz`.
3. **Mask to N** every reference base where any aligned species differs (mismatch), and the
   bases around every insertion or deletion. Unaligned regions stay as in the reference.

The result has the reference's length and coordinates, so its gene annotation (GTF) still
applies, but no species' reads are penalised at divergent sites. The
[`consensus_genome`](../pipelines/README.md) pipeline implements this. Its masking script
(`mask_divergent.py`) takes the flank size around indels (`--indel-flank`, default 6 bp) and
optionally around mismatches (`--mismatch-flank`, default 0), and writes a BED of every masked
interval.

```bash
cd pipelines/consensus_genome
$EDITOR config.yaml               # target: hg38; queries: [panTro6, rheMac10]
snakemake --use-conda --cores 8
cat results/hg38_panTro6_rheMac10.masking_summary.txt     # how much of the genome was masked
```

**Checks and trade-offs**

- **How much is masked?** Roughly the divergence to the most distant species plus indel flanks.
  Human and rhesus macaque genomes are ~93 % identical
  ([Rhesus Macaque Genome Consortium 2007](https://doi.org/10.1126/science.1139247)), and human and
  chimpanzee ~99 %. Expect on the order of 5–10 % of aligned hg38 bases to be masked for
  human/chimp/macaque, then check the pipeline's `masking_summary.txt`. More distant species mask
  much more, so mappability falls and the approach stops paying off.
- **Reads spanning masked bases** align with lower scores or not at all. Every species is
  affected equally, which is the point, but total mapping rates drop a little.
- **Aligners treat N differently:** bowtie2 applies a penalty to N positions (`--np`,
  default 1), and STAR and BWA have their own rules. Keep the aligner and settings identical for
  all species.
- **Species-specific sequence** (insertions relative to the reference) has no place in the
  consensus. Features inside it are invisible by design. Use the own-genome + liftOver route to
  study them.

## 5. Which lineage changed? Outgroup logic

A difference between two species doesn't say which one changed. With a third, more distantly
related species (an **outgroup**), a change is assigned to the lineage of species A only when A
differs from **both** B and the outgroup in the same direction, while B and the outgroup agree:

```r
# lfc_AB, lfc_AO, lfc_BO: log2 fold changes (A vs B, A vs outgroup, B vs outgroup); padj_*: FDR
a_lineage_up <- padj_AB < 0.01 & padj_AO < 0.01 & lfc_AB > 0 & lfc_AO > 0
# stricter: B and the outgroup don't differ
a_lineage_up_strict <- a_lineage_up & (is.na(padj_BO) | padj_BO > 0.1)
```

Caveats: with one individual of the outgroup, individual variation and species differences
can't be separated. More outgroup individuals (or species) make the assignment robust. The same
logic applies to peaks (ATAC/ChIP), loops (Hi-C) and pseudobulk cell types (scRNA-seq).

## Example study

The guides use Ciuba *et al.* 2025, *Cell Stem Cell* 32:426–444
([doi:10.1016/j.stem.2024.12.011](https://doi.org/10.1016/j.stem.2024.12.011)) as a worked
example. It aligned RNA-seq from human, chimpanzee and macaque iPSC-derived astrocytes to an
hg38-based consensus genome masked against panTro6 and rheMac10. For ATAC-seq and ChIP-seq it
mostly aligned each species to its own genome and compared peaks via liftOver, using the consensus
genome only where an analysis needed one shared coordinate system. It assigned expression and
enhancer changes to the human lineage with the outgroup logic above. Raw data: ArrayExpress
E-MTAB-13252 (RNA-seq), E-MTAB-13253 (ATAC-seq), E-MTAB-13254/-13255 (H3K27ac/H3K4me3
ChIP-seq), E-MTAB-13259 (Hi-C).
