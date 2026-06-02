# Generate the committed golden variance files (scalar-df + per-gene-df cases).
suppressMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
set.seed(2004)
G <- 400
sigma2 <- 1.3 * 8 / rchisq(G, df = 8)
df <- 6
v <- sigma2 * rchisq(G, df = df) / df
genes <- sprintf("g%04d", seq_len(G))
write.table(data.frame(gene = genes, variance = v), file.path(dir, "vars.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# per-gene df case
set.seed(11)
dfv <- sample(3:9, G, replace = TRUE)
vv <- sigma2 * rchisq(G, df = dfv) / dfv
write.table(data.frame(gene = genes, variance = vv), file.path(dir, "vars_pg.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = genes, df = dfv), file.path(dir, "df_pg.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
