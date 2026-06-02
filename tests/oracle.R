# Regenerate the squeezeVar oracle TSV from a variances file.
# Usage: oracle.R <vars.tsv> <df-or-dffile> <out.tsv>
suppressMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
varfile <- args[1]; dfarg <- args[2]; out <- args[3]
tab <- read.table(varfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- tab[[1]]; v <- tab[[2]]
if (file.exists(dfarg)) {
  dft <- read.table(dfarg, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df <- dft[[2]]
} else {
  df <- as.numeric(dfarg)
}
o <- squeezeVar(v, df, legacy = TRUE)
dp <- if (is.infinite(o$df.prior)) "Inf" else format(o$df.prior, digits = 17)
vp <- format(o$var.prior, digits = 17)
con <- file(out, "w")
writeLines(sprintf("# var.prior=%s df.prior=%s", vp, dp), con)
writeLines("gene\tvar.post", con)
writeLines(paste(genes, format(o$var.post, digits = 17), sep = "\t"), con)
close(con)
