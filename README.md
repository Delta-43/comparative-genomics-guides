<h1 align="center">Functional Genomics Guides</h1>

<p align="center"><i>Reference guides and ready-to-run Snakemake pipelines that document one way of taking RNA-seq, ATAC-seq, ChIP-seq, Hi-C and single-cell RNA-seq data from raw reads to biological inference, with the reasoning behind each choice.</i></p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/Content-CC%20BY%204.0-lightgrey.svg" alt="Content licence: CC BY 4.0"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Code-MIT-blue.svg" alt="Code licence: MIT"></a>
  <a href="https://github.com/Delta-43/comparative-genomics-guides/commits/main"><img src="https://img.shields.io/github/last-commit/Delta-43/comparative-genomics-guides" alt="Last commit"></a>
  <a href="https://github.com/Delta-43/comparative-genomics-guides/graphs/contributors"><img src="https://img.shields.io/github/contributors-anon/Delta-43/comparative-genomics-guides" alt="Contributors"></a>
</p>

<p align="center">
  <a href="https://delta-43.github.io/comparative-genomics-guides/"><img src="https://img.shields.io/badge/Read%20the%20guides-online-2ea44f" alt="Read the guides online"></a>
  <a href="https://github.com/Delta-43/comparative-genomics-guides/actions/workflows/publish.yml"><img src="https://img.shields.io/github/actions/workflow/status/Delta-43/comparative-genomics-guides/publish.yml?branch=main&label=site" alt="Site build"></a>
</p>

<p align="center">
  <a href="#-what-it-is">What it is</a> •
  <a href="#-who-its-for">Who it's for</a> •
  <a href="#-author">Author</a> •
  <a href="#-the-guides">The guides</a> •
  <a href="#-project-status">Status</a> •
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-tech-stack">Tech Stack</a> •
  <a href="#-project-structure">Structure</a> •
  <a href="#-contributing">Contributing</a> •
  <a href="#-license">License</a>
</p>

> [!TIP]
> You don't need to install anything to start reading. The guides are easiest to read on the
> website: **https://delta-43.github.io/comparative-genomics-guides/**

---

## 🧬 What it is

Each assay gets two things. The first is a **guide**: a page that walks through experimental
design, raw-read QC, processing, quality checkpoints, statistics, plots and interpretation, with
R code you can run. The second is a **Snakemake pipeline** in [`pipelines/`](pipelines/) that turns your
FASTQ files into the processed files each guide's analysis starts from.

- **One documented path per assay** — each guide follows one workflow built from widely used tools (for example STAR + featureCounts + DESeq2 for RNA-seq). At each key decision point, a pros/cons box names the main alternative (Salmon, cooler, Seurat, ...) and says when to prefer it.
- **Quality checkpoints with numbers** — named metrics (strandedness, FRiP, TSS enrichment, fragment sizes, cis/trans ratio, mitochondrial fraction) with thresholds drawn from ENCODE and tool documentation. They tell you whether your data is good enough to continue.
- **Inference, not just commands** — every guide ends with what your results can and can't support, common pitfalls, and troubleshooting.
- **Tested R code** — every R block in the guides has been executed on public example data (airway, PBMC 3k and others) under R 4.5 / Bioconductor 3.22.
- **Comparative extension** — each guide first covers the single-species analysis, then ends with a cross-species module (liftOver, consensus genomes, outgroup logic). An appendix covers the shared comparative-genomics groundwork.
- **Worked example on open data** — a published study of primate astrocyte evolution (Ciuba et al. 2025, *Cell Stem Cell*, [doi:10.1016/j.stem.2024.12.011](https://doi.org/10.1016/j.stem.2024.12.011)) and its ArrayExpress data illustrate what each analysis produces. It is an example, not a target to reproduce.

## 🎓 Who it's for

Researchers who have, or are planning, a sequencing experiment and want a documented, worked
example of one way to analyse it, to compare with their own approach and adapt. You don't have to
be a bioinformatician.

**A reference, not a standard.** The guides record [the author's](#-author) workflow. Other valid
approaches exist, and the right settings depend on your data, so check each tool's own
documentation before relying on a choice made here.

**What you need to know first:**

- **A terminal** — moving between folders, running a command, editing a text file.
- **Some R** — loading a package, reading a table, calling a function. The guides explain the rest.

If either is new to you, the free, beginner-friendly [Harvard Chan Bioinformatics Core
training materials](https://hbctraining.github.io/main/) cover both. They inspired the structure
of these guides, which link to them wherever they go deeper.

## 👤 Author

These guides were created by **Debadeep Chaudhury**. They are based on five years of experience as a
bioinformatician during a PhD, and on hands-on research experience, written up so that others can
follow, compare and adapt the same workflows.

<p>
  <a href="https://orcid.org/0000-0002-9089-732X"><img src="https://img.shields.io/badge/ORCID-0000--0002--9089--732X-A6CE39?logo=orcid&logoColor=white" alt="ORCID 0000-0002-9089-732X"></a>
  <a href="https://scholar.google.com/citations?user=zDcS-78AAAAJ"><img src="https://img.shields.io/badge/Google%20Scholar-profile-4285F4?logo=googlescholar&logoColor=white" alt="Google Scholar profile"></a>
  <a href="https://www.linkedin.com/in/dchaudhury"><img src="https://img.shields.io/badge/LinkedIn-dchaudhury-0A66C2" alt="LinkedIn profile"></a>
</p>

- **ORCID** — [0000-0002-9089-732X](https://orcid.org/0000-0002-9089-732X)
- **Google Scholar** — [publications and citations](https://scholar.google.com/citations?user=zDcS-78AAAAJ)
- **LinkedIn** — [linkedin.com/in/dchaudhury](https://www.linkedin.com/in/dchaudhury)

## 📚 The guides

Start with **Getting started**, then pick the assay you're working with.

| Guide | The question it answers | Pipeline |
|---|---|---|
| [Getting started](guides/00-getting-started.md) | How do I set up software, reference genomes and example data? | — |
| [RNA-seq](guides/01-rna-seq.md) | Which genes change expression? | [`rnaseq`](pipelines/rnaseq/) |
| [ATAC-seq](guides/02-atac-seq.md) | Which regulatory elements open or close? | [`chip_atac`](pipelines/chip_atac/) (`assay: atac`) |
| [ChIP-seq](guides/03-chip-seq.md) | Where does a protein bind or a histone mark sit, and does it change? | [`chip_atac`](pipelines/chip_atac/) (`assay: chip`) |
| [Hi-C](guides/04-hic.md) | How is the genome folded, and does folding change? | [`hic`](pipelines/hic/) |
| [scRNA-seq](guides/05-scrna-seq.md) | Which cell types exist, and how do they differ? | [`scrna`](pipelines/scrna/) |
| [Comparative genomics](guides/appendix-comparative-genomics.md) | How do I compare species fairly? | [`consensus_genome`](pipelines/consensus_genome/) |

## 📊 Project Status

| Component | Status | Detail |
|---|---|---|
| Guides (R code) | ![Tested](https://img.shields.io/badge/-Tested-brightgreen) | Every R block executed on public example data (R 4.5, Bioconductor 3.22). |
| `rnaseq`, `chip_atac` | ![Tested](https://img.shields.io/badge/-Tested-brightgreen) | Run end to end on simulated data. They recovered the simulated differential genes and 40 of 40 simulated peaks. |
| `hic`, `scrna`, `consensus_genome` | ![Dry run](https://img.shields.io/badge/-Dry--run%20tested-yellow) | Every mode dry-runs and every environment builds. They need full-size inputs (Juicer, Cell Ranger, UCSC genomes) that haven't been run yet. |

No pipeline has yet been run on a full-size real dataset. If you run one, please report what
happened in the [issue tracker](https://github.com/Delta-43/comparative-genomics-guides/issues).
Details are in [`pipelines/README.md`](pipelines/README.md).

## 🚀 Quick Start

**1. Install conda and Snakemake once.** [Miniforge](https://github.com/conda-forge/miniforge)
gives you `conda` and `mamba`. Then:

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake=9
conda activate snakemake
git clone https://github.com/Delta-43/comparative-genomics-guides.git
cd comparative-genomics-guides
```

**2. Set up R for the analysis sections.** This one environment has every R and Bioconductor
package the guides use:

```bash
mamba env create -n guides-r -f pipelines/envs/r_analysis.yaml
```

**3. Process your own data.** Edit the config and sample sheet, dry-run, then run. Each rule
installs its own tools on first use (`--use-conda`), so you never install them by hand.

```bash
cd pipelines/rnaseq                    # or chip_atac, hic, scrna, consensus_genome
$EDITOR config.yaml samples.tsv        # reference files, options, your samples
snakemake --use-conda --cores 16 -n    # dry run: lists every step. Always do this first.
snakemake --use-conda --cores 16       # run
```

Processing needs a workstation or a cluster (e.g. 16 cores and 64 GB RAM for RNA-seq with STAR
on human). The R analysis runs on a laptop. See [Getting started](guides/00-getting-started.md) for
hardware, reference genomes and example data. See [`pipelines/README.md`](pipelines/README.md) for
cluster profiles, Cell Ranger and Juicer set-up.

## 🧰 Tech Stack

<p align="center">
  <a href="https://snakemake.readthedocs.io"><img src="https://img.shields.io/badge/Snakemake-9-039475" alt="Snakemake 9"></a>
  <a href="https://www.r-project.org"><img src="https://img.shields.io/badge/R-4.5-276DC3?logo=r&logoColor=white" alt="R 4.5"></a>
  <a href="https://bioconductor.org"><img src="https://img.shields.io/badge/Bioconductor-3.22-1a81c2" alt="Bioconductor 3.22"></a>
  <a href="https://docs.conda.io"><img src="https://img.shields.io/badge/conda-bioconda-44A833?logo=anaconda&logoColor=white" alt="conda / bioconda"></a>
  <a href="https://quarto.org"><img src="https://img.shields.io/badge/Quarto-website-75AADB?logo=quarto&logoColor=white" alt="Quarto"></a>
</p>

The main tools, with versions pinned in [`pipelines/envs/`](pipelines/envs/): FastQC, Trim Galore,
STAR, featureCounts, Bowtie2, MACS2, deepTools, Juicer, Mustache, Cell Ranger and the UCSC
chain/net tools. The R analysis uses DESeq2, DiffBind, ChIPseeker, clusterProfiler, plotgardener,
scran/scater and muscat.

## 🗂️ Project Structure

```text
comparative-genomics-guides/
├── guides/                 # the guides (Markdown): 00 getting started, 01–05 one per assay, comparative appendix
├── pipelines/              # Snakemake pipelines, one folder per assay
│   ├── rnaseq/             #   each has a Snakefile, config.yaml and an example sample sheet
│   ├── chip_atac/
│   ├── hic/
│   ├── scrna/
│   ├── consensus_genome/
│   ├── envs/               #   pinned conda environments used by the rules (+ r_analysis.yaml for the guides)
│   ├── demo_data/          #   ENA run tables for the example study + make_samples.py
│   └── README.md           #   install, run, cluster use and testing status
├── index.qmd, about.qmd    # website landing and About pages
├── _quarto.yml, assets/    # website configuration and styling
├── tools/build_site.sh     # builds the website into _site/
└── .github/workflows/      # publishes the website to GitHub Pages on every push to main
```

## 🤝 Contributing

Corrections, questions and reports from real runs are all welcome. Please open an
[issue](https://github.com/Delta-43/comparative-genomics-guides/issues). Reports on the parts
marked above as not yet run on full-size data are especially useful. For changes, open a pull
request against `main`. If you edit an R block in a guide, please run it first. To preview the
website locally, build it with `tools/build_site.sh` (it needs [Quarto](https://quarto.org)).

## 📄 License

The guides and website content are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You may reuse and adapt them for any
purpose, including teaching, as long as you credit the source. The code is licensed under
[MIT](LICENSE). See [`LICENSE`](LICENSE).

---

<p align="center">
  <a href="https://github.com/Delta-43/comparative-genomics-guides/graphs/contributors">
    <img src="https://contrib.rocks/image?repo=Delta-43/comparative-genomics-guides" alt="Contributors">
  </a>
</p>

<p align="center"><sub>Content licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>; code under <a href="LICENSE">MIT</a>.</sub></p>
