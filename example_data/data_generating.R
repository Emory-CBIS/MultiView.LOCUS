# Generate example simulation datasets in a simple script flow:
#   1) simu1 -> Data_x_simu1.rds, Data_y_simu1.rds
#   2) simu2 -> Data_x_simu2.rds, Data_y_simu2.rds
#   3) simu3(Real) -> Data_x.rds,      Data_y.rds
#
# Run from repository root:
#   source("example_data/generate_example_data.R")

set.seed(111)

source_dir <- file.path("example_data", "S")
output_dir <- "example_data"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ICcorr is precomputed in this repo.
ICcorr <- readRDS(file.path("data", "ICcorr.rds"))

# -------------------------
# simu1
# -------------------------
common_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "simulated_common_SourceMar7.RData"), envir = common_env)

specific_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "simulatedspecificSourceMar6.RData"), envir = specific_env)

sample1 <- sample(1:237, 100)

mixing_cor_commonx <- t(ICcorr$M[sample1, c(1, 3, 7)]) / apply(ICcorr$M[sample1, c(1, 3, 7)], 2, sd)
mixing_cor_commony <- mixing_cor_commonx + matrix(rnorm(300, 0, sqrt(3)), nrow = 3)
mixing_cor_commonx <- mixing_cor_commonx * apply(ICcorr$M[sample1, c(1, 3, 7)], 2, sd)
mixing_cor_commony <- mixing_cor_commony * apply(ICcorr$M[sample1, c(1, 3, 7)], 2, sd)

mixing_x <- t(rbind(mixing_cor_commonx, t(ICcorr$M[sample1, 8:9])))
mixing_y <- t(rbind(mixing_cor_commony, t(ICcorr$M[sample1, 10:11])))

if (!requireNamespace("MultiView.LOCUS", quietly = TRUE)) {
  stop("Package 'MultiView.LOCUS' is required for Ltrans() in simu1 generation.")
}

S_s <- rbind(
  MultiView.LOCUS::Ltrans(specific_env$S_s1),
  MultiView.LOCUS::Ltrans(specific_env$S_s2),
  MultiView.LOCUS::Ltrans(specific_env$S_s3),
  MultiView.LOCUS::Ltrans(specific_env$S_s4)
)

noise <- 0.1
Data_x <- mixing_x %*% rbind(common_env$S_x, S_s[1:2, ]) + matrix(rnorm(122500, sd = noise), nrow = 100)
Data_y <- mixing_y %*% rbind(common_env$S_y, S_s[3:4, ]) + matrix(rnorm(122500, sd = noise), nrow = 100)

saveRDS(Data_x, file = file.path(output_dir, "Data_x_simu1.rds"))
saveRDS(Data_y, file = file.path(output_dir, "Data_y_simu1.rds"))

# -------------------------
# simu2
# -------------------------
set.seed(111)
noise_level <- 0.1
correlation_level <- 0.4

simu2_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "simulated_sources_type2.Rdata"), envir = simu2_env)

sample1 <- sample(1:237, 100)

mixing_cor_commonx <- t(ICcorr$M[sample1, c(2, 3)]) / apply(ICcorr$M[sample1, c(2, 3)], 2, sd)
mixing_cor_commony <- mixing_cor_commonx + matrix(rnorm(200, 0, sqrt(1 / correlation_level^2 - 1)), nrow = 2)
mixing_cor_commonx <- mixing_cor_commonx * apply(ICcorr$M[sample1, c(2, 3)], 2, sd)
mixing_cor_commony <- mixing_cor_commony * apply(ICcorr$M[sample1, c(2, 3)], 2, sd)

mixing_x <- t(rbind(mixing_cor_commonx, t(ICcorr$M[sample1, 8])))
mixing_y <- t(rbind(mixing_cor_commony, t(ICcorr$M[sample1, 9])))

Data_x <- mixing_x %*% simu2_env$S_x + matrix(rnorm(122500, sd = noise_level), nrow = 100)
Data_y <- mixing_y %*% simu2_env$S_y + matrix(rnorm(122500, sd = noise_level), nrow = 100)

saveRDS(Data_x, file = file.path(output_dir, "Data_x_simu2.rds"))
saveRDS(Data_y, file = file.path(output_dir, "Data_y_simu2.rds"))

# -------------------------
# simu3 (mid settings only)
# -------------------------
set.seed(111)
noise_level <- 0.05
correlation_level <- 0.4

S_x <- readRDS(file.path(source_dir, "S_xreal.rds"))
S_y <- readRDS(file.path(source_dir, "S_yreal.rds"))

sample1 <- sample(1:237, 100)

mixing_cor_commonx <- t(ICcorr$M[sample1, c(5, 6)]) / apply(ICcorr$M[sample1, c(5, 6)], 2, sd)
mixing_cor_commony <- mixing_cor_commonx + matrix(rnorm(200, 0, sqrt(1 / correlation_level^2 - 1)), nrow = 2)

mixing_x <- t(rbind(mixing_cor_commonx, t(ICcorr$M[sample1, 8:9])))
mixing_y <- t(rbind(mixing_cor_commony, t(ICcorr$M[sample1, 10:11])))

Data_x <- mixing_x %*% S_x + matrix(rnorm(3471600, sd = noise_level), nrow = 100)
Data_y <- mixing_y %*% S_y + matrix(rnorm(3471600, sd = noise_level), nrow = 100)

# Match existing file names used in example_data/
saveRDS(Data_x, file = file.path(output_dir, "Data_x.rds"))
saveRDS(Data_y, file = file.path(output_dir, "Data_y.rds"))
