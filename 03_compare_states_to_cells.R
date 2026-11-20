# 03_compare_states_to_cells.R
# Compare sRACIPE steady-state clusters with the Specialist / Generalist / Transitory myeloid cells (EBUS_06).
# usage: Rscript 03_compare_states_to_cells.R <tag> <k>      e.g.  Rscript 03_compare_states_to_cells.R primary 0
#   k = 0 -> choose the number of state clusters by average silhouette (2..8); k >= 2 -> use that k
suppressPackageStartupMessages(library(cluster))
a   <- commandArgs(trailingOnly = TRUE)
tag <- if (length(a) >= 1) a[1] else "primary"
kfx <- if (length(a) >= 2) as.integer(a[2]) else 0L
set.seed(1)

sim  <- read.csv(paste0("racipe_states_", tag, ".csv"), row.names = 1, check.names = FALSE)   # states x nodes (z-scored log2)
cell <- read.csv("myeloid_node_table_EBUS06.csv", row.names = 1, check.names = FALSE)
grp  <- factor(cell$group, levels = c("Specialist", "Generalist", "Transitory")); cell$group <- NULL

# 1. keep only nodes that are in the circuit AND reasonably detected in the cells (>= 10% of cells with value > 0)
det_pct <- colMeans(cell > 0) * 100
nodes   <- intersect(colnames(sim), names(det_pct)[det_pct >= 10])
cat("nodes used (detected in >= 10% of cells):", paste(nodes, collapse = ", "), "\n")
cat("nodes dropped (too sparse or not in the data):", paste(setdiff(colnames(sim), nodes), collapse = ", "), "\n")
stopifnot(length(nodes) >= 5)

# 2. cluster the simulated steady states on those nodes
S <- as.matrix(sim[, nodes]); S <- S[complete.cases(S), ]
hc <- hclust(dist(S), method = "ward.D2")
if (kfx >= 2) { k <- kfx } else {
  ks <- 2:8; d <- dist(S)
  sil <- sapply(ks, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))
  k <- ks[which.max(sil)]; cat("silhouette by k:", paste(ks, round(sil, 3), sep = "=", collapse = "  "), "-> k =", k, "\n")
}
cl <- cutree(hc, k)
cent <- t(sapply(1:k, function(i) colMeans(S[cl == i, , drop = FALSE])))
rownames(cent) <- paste0("state", 1:k)
cat("state-cluster sizes:", paste(table(cl), collapse = " / "), "\n")
write.csv(round(cent, 2), paste0("racipe_cluster_centroids_", tag, ".csv"))

# 3. score every cell against every simulated-state centroid (Pearson over the shared nodes, cells z-scored per node)
C <- scale(as.matrix(cell[, nodes])); C[is.na(C)] <- 0
R <- cor(t(C), t(cent))                        # cells x states
best <- colnames(R)[max.col(R, ties.method = "first")]
srt  <- t(apply(R, 1, sort, decreasing = TRUE))
margin <- srt[, 1] - srt[, 2]

tab <- table(group = grp, best_state = factor(best, levels = rownames(cent)))
cat("\ncells per group x best-matching RACIPE state:\n"); print(tab)
cat("\nrow %:\n"); print(round(100 * prop.table(tab, 1), 0))
mean_r <- aggregate(as.data.frame(R), list(group = grp), mean)
cat("\nmean correlation of each group's cells with each state centroid:\n"); print(format(mean_r, digits = 2))
cat("\nmedian margin between best and second-best state (small = intermediate):\n")
print(round(tapply(margin, grp, median), 3))
cat("Kruskal-Wallis p, margin by group:", signif(kruskal.test(margin ~ grp)$p.value, 3), "\n")

write.csv(data.frame(cell = rownames(cell), group = grp, best_state = best, margin = margin, round(R, 3)),
          paste0("cell_state_assignment_", tag, ".csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(tab), paste0("group_by_state_", tag, ".csv"))
pdf(paste0("fig_states_vs_groups_", tag, ".pdf"), width = 9, height = 5)
op <- par(mfrow = c(1, 2))
image(t(cent[nrow(cent):1, ]), axes = FALSE, col = hcl.colors(50, "RdBu", rev = TRUE), main = "RACIPE state centroids (z)")
axis(1, at = seq(0, 1, length.out = length(nodes)), labels = nodes, las = 2, cex.axis = 0.7)
axis(2, at = seq(0, 1, length.out = k), labels = rev(rownames(cent)), las = 1, cex.axis = 0.7)
barplot(t(100 * prop.table(tab, 1)), beside = FALSE, legend.text = TRUE, ylab = "% of group's cells", main = "Best-matching state",
        args.legend = list(x = "topright", cex = 0.6))
dev.off()
