#!/usr/bin/env Rscript
# Juicer `dump observed NONE ... BP <res>` sparse triplets -> TopDom dense matrix -> TopDom domains.
# Usage: topdom_from_dump.R <dump.txt> <chrom> <chrom_length> <bin_size> <window_size> <out_prefix>
# Writes <out_prefix>.topdom.bed (domain/gap/boundary rows) and <out_prefix>.topdom.bedpe (domains only).
suppressPackageStartupMessages({ library(data.table); library(TopDom) })
a <- commandArgs(trailingOnly = TRUE)
dump <- a[1]; chrom <- a[2]; clen <- as.numeric(a[3]); bs <- as.numeric(a[4]); ws <- as.integer(a[5]); out <- a[6]

n <- ceiling(clen / bs)
d <- fread(dump, col.names = c("x", "y", "v"))
d[, `:=`(i = x %/% bs + 1L, j = y %/% bs + 1L)]
m <- matrix(0, n, n)
m[cbind(d$i, d$j)] <- d$v
m[cbind(d$j, d$i)] <- d$v                      # dump is upper-triangular; symmetrise
m[is.na(m)] <- 0
starts <- (seq_len(n) - 1) * bs
dense <- data.table(chr = chrom, from.coord = starts, to.coord = pmin(starts + bs, clen), m)
tmp <- paste0(out, ".dense.tmp")
fwrite(dense, tmp, sep = "\t", col.names = FALSE)

td <- TopDom(tmp, window.size = ws)
unlink(tmp)
bed <- td$bed
fwrite(bed, paste0(out, ".topdom.bed"), sep = "\t", col.names = FALSE)
dom <- bed[bed$name == "domain", ]
fwrite(data.table(dom$chrom, dom$chromStart, dom$chromEnd, dom$chrom, dom$chromStart, dom$chromEnd),
       paste0(out, ".topdom.bedpe"), sep = "\t", col.names = FALSE)
