# Cell Ranger input templates

- `library.csv`: one row per library (FASTQ folder + `sample` prefix + `library_type`). Paths must
  be absolute. Re-sequenced runs of the same library: add one row per FASTQ folder, same `sample`.
- `feature_ref.csv`: antibody barcodes. These rows are BioLegend TotalSeq-B human hashtags 1–4,
  copied from 10x Genomics' own list (`lib/rust/cr_lib/src/stages/data/TotalSeq_Hashtag.csv` in
  github.com/10XGenomics/cellranger). Replace them with your panel's barcodes. The `pattern`
  depends on the chemistry (TotalSeq-B: `5PNNNNNNNNNN(BC)`), so check your kit's documentation.
- `hashtag_to_sample.csv`: your own lookup from hashtag `id` to sample, used in R for demultiplexing.
