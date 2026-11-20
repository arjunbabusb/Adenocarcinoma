# 04_robustness.R -- how fragile is the "transitory cells sit in state X / between states" signal?
# usage: Rscript 04_robustness.R <tag>        e.g.  Rscript 04_robustness.R primary
suppressPackageStartupMessages(library(cluster))
tag <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "primary" }
set.seed(7)
sim  <- read.csv(paste0("racipe_states_", tag, ".csv"), row.names = 1, check.names = FALSE)
cell <- read.csv("myeloid_node_table_EBUS06.csv", row.names = 1, check.names = FALSE)
grp  <- factor(cell$group, levels = c("Specialist", "Generalist", "Transitory")); cell$group <- NULL
det  <- colMeans(cell > 0) * 100
base_nodes <- intersect(colnames(sim), names(det)[det >= 10])

run <- function(nodes, k, sub = 4000) {
  S <- as.matrix(sim[, nodes]); S <- S[complete.cases(S), ]
  S <- S[sample(nrow(S), min(sub, nrow(S))), ]                 # subsample states so hclust stays fast
  cl <- cutree(hclust(dist(S), "ward.D2"), k)
  cent <- t(sapply(1:k, function(i) colMeans(S[cl == i, , drop = FALSE])))
  C <- scale(as.matrix(cell[, nodes])); C[is.na(C)] <- 0
  R <- suppressWarnings(cor(t(C), t(cent)))
  srt <- t(apply(R, 1, sort, decreasing = TRUE))
  # orient the clusters so numbering is stable: state 1 = the cluster with the highest mean NFKB + STAT1 (inflammatory-like),
  # the others follow in decreasing order. Without this, hclust numbering changes from run to run.
  ord <- order(-rowMeans(cent[, intersect(c("NFKB", "STAT1"), nodes), drop = FALSE]))
  cent <- cent[ord, , drop = FALSE]; R <- R[, ord, drop = FALSE]
  list(cent = cent, best = max.col(R, ties.method = "first"), margin = srt[, 1] - srt[, 2])
}
report <- function(label, o) {
  tab <- table(grp, factor(o$best, levels = seq_len(nrow(o$cent))))
  p_fisher <- if (ncol(tab) == 2) fisher.test(tab[c("Transitory", "Specialist"), ])$p.value else chisq.test(tab[c("Transitory", "Specialist"), ], simulate.p.value = TRUE)$p.value
  mw <- wilcox.test(o$margin[grp == "Transitory"], o$margin[grp != "Transitory"])$p.value
  cat(sprintf("%-34s transitory share in each state: %s | vs specialist p=%.3f | median margin transitory %.2f vs others %.2f (MW p=%.3f)\n",
              label, paste(round(100 * prop.table(tab["Transitory", ])), collapse = "/"), p_fisher,
              median(o$margin[grp == "Transitory"]), median(o$margin[grp != "Transitory"]), mw))
}
cat("== 1. State centroids (z), k = 2: read these to name the states ==\n")
o2 <- run(base_nodes, 2); rownames(o2$cent) <- paste0("state", 1:2, c(" (NFKB/STAT1-high)", " (other)")); print(round(o2$cent, 2))
cat(sprintf("correlation between the two centroids: %.3f (about -1 means one axis, not two independent states)\n", cor(o2$cent[1, ], o2$cent[2, ])))
cat("\n== 2. Different number of state clusters ==\n")
for (k in 2:4) report(paste0("k = ", k, ", all nodes"), run(base_nodes, k))
cat("\n== 3. Leave one node out (k = 2). Does the transitory pattern survive? ==\n")
for (n in base_nodes) report(paste0("without ", n), run(setdiff(base_nodes, n), 2))
cat("\n== 4. Permutation null: shuffle group labels 1000x, how often is a 25-cell group this skewed? ==\n")
o <- run(base_nodes, 2); st <- o$best; obs <- mean(st[grp == "Transitory"] == 1)
perm <- replicate(1000, { g <- sample(grp); mean(st[g == "Transitory"] == 1) })
cat(sprintf("observed share of transitory in state 1 (NFKB/STAT1-high) = %.2f; permutation two-sided p = %.3f\n", obs, mean(abs(perm - mean(st == 1)) >= abs(obs - mean(st == 1)))))
cat("\nNote: each run re-clusters a random subsample of states, so numbers wobble slightly from run to run; look at the pattern, not decimals.\n")
