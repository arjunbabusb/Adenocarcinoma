# 01_run_racipe.R -- simulate the macrophage circuit with sRACIPE (Bioconductor)
# usage:  Rscript 01_run_racipe.R <topology.csv> <tag> <numModels> <nCores>
#   e.g.  Rscript 01_run_racipe.R topology_primary_boolean.csv primary 2000 4
suppressPackageStartupMessages({library(sRACIPE); library(SummarizedExperiment)})
set.seed(42)
a <- commandArgs(trailingOnly = TRUE)
topo_file <- if (length(a) >= 1) a[1] else "topology_primary_boolean.csv"
tag       <- if (length(a) >= 2) a[2] else "primary"
n_models  <- if (length(a) >= 3) as.integer(a[3]) else 2000L
n_cores   <- if (length(a) >= 4) as.integer(a[4]) else 1L

topo <- read.csv(topo_file, stringsAsFactors = FALSE)   # Source, Target, Interaction (1 = activation, 2 = inhibition)
stopifnot(all(c("Source","Target","Interaction") %in% names(topo)))

t0 <- Sys.time()
rSet <- sracipeSimulate(circuit = topo, numModels = n_models,
                        nIC = 5,               # random starting points per model -> finds several stable states per model
                        limitcycles = FALSE,   # skip the (slow) limit-cycle search; we only want stable steady states
                        plots = FALSE, nCores = n_cores)
rSet <- sracipeNormalize(rSet)                  # log2 + z-score each gene across all states
expr <- assay(rSet, 1)                          # genes x (model x stable-state) matrix
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " genes:", nrow(expr), " state-columns:", ncol(expr), "\n")
saveRDS(rSet, paste0("rSet_", tag, ".rds"))
write.csv(t(expr), paste0("racipe_states_", tag, ".csv"))   # rows = simulated states, columns = genes (z-scored log2)
