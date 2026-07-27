# Exploratory benchmark workflow. Requires external qgg data and is not
# intended to run during package checks.
#
# Set SBLR_EXAMPLE_DATA_DIR before running this workflow, or edit the fallback
# path below to point to a local directory for the example data.

library(qgg)

data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR", unset = "path/to/example/data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

files <- c("bed", "bim", "fam", "pheno", "covar")

for (f in files) {
 url <- sprintf(
  "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.%s",
  f
 )
 download.file(url, destfile = file.path(data_dir, paste0("human.", f)), mode = "wb")
}

Glist <- gprep(
 study = "Example",
 bedfiles = file.path(data_dir, "human.bed"),
 bimfiles = file.path(data_dir, "human.bim"),
 famfiles = file.path(data_dir, "human.fam")
)

rsids <- gfilter(Glist = Glist, excludeMAF = 0.05, excludeMISS = 0.05,
                 excludeCGAT = TRUE, excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
                 excludeMHC = FALSE)

ldfiles <- file.path(data_dir, "human.ld")
Glist <- gprep(Glist, task = "sparseld", msize = 1000, rsids = rsids, ldfiles = ldfiles,
               overwrite = TRUE)
saveRDS(Glist, file = file.path(data_dir, "Glist_sparseLD_1k.RDS"),
        compress = FALSE)

Glist <- readRDS(file = file.path(data_dir, "Glist_sparseLD_1k.RDS"))

B <- getLD(Glist, chr=1)
ld <- getSparseLD(Glist=Glist, chr=1)

library(RhpcBLASctl)

RhpcBLASctl::blas_set_num_threads(4)
RhpcBLASctl::omp_set_num_threads(2)

system.time(out <- sparseLD_stream_CSR(
 bed_file = Glist$bedfiles[chr],
 n = Glist$n,
 cls = seq_along(Glist$rsids[[chr]]),
 af = Glist$af[[chr]],
 out_prefix = file.path(
  data_dir,
  sprintf("ld_chr%d", chr)
 ),
 pos_bp = Glist$pos[[chr]],
 max_distance_bp = 0,            # disables bp-distance filtering
 max_distance_variants = 5000,
 r2_threshold = 0.01,
 block_size = 512,
 nthreads = 2
))

bench_mkl <- expand.grid(
 blas_threads = c(1, 2, 4, 6, 8, 12),
 omp_threads  = c(1, 2),
 block_size   = c(512, 1024),
 stringsAsFactors = FALSE
)

bench_mkl$elapsed <- NA_real_
bench_mkl$user <- NA_real_
bench_mkl$system <- NA_real_
bench_mkl$nnz <- NA_real_

for (k in seq_len(nrow(bench_mkl))) {
 RhpcBLASctl::blas_set_num_threads(bench_mkl$blas_threads[k])
 RhpcBLASctl::omp_set_num_threads(bench_mkl$omp_threads[k])

 prefix <- file.path(
  data_dir,
  sprintf(
   "ld_chr%d_blas%d_omp%d_b%d",
   chr,
   bench_mkl$blas_threads[k],
   bench_mkl$omp_threads[k],
   bench_mkl$block_size[k]
  )
 )

 tm <- system.time({
  out_k <- sparseLD_stream_CSR(
   bed_file = Glist$bedfiles[chr],
   n = Glist$n,
   cls = seq_along(Glist$rsids[[chr]]),
   af = Glist$af[[chr]],
   out_prefix = prefix,
   pos_bp = Glist$pos[[chr]],
   max_distance_bp = 0,            # disables bp-distance filtering
   max_distance_variants = 5000,
   r2_threshold = 0.001,
   block_size = bench_mkl$block_size[k],
   nthreads = bench_mkl$omp_threads[k]
  )
 })

 bench_mkl$user[k] <- tm[["user.self"]]
 bench_mkl$system[k] <- tm[["sys.self"]]
 bench_mkl$elapsed[k] <- tm[["elapsed"]]
 bench_mkl$nnz[k] <- out_k$nnz
}

bench_mkl[order(bench_mkl$elapsed), ]


csr <- sparseLD_read_CSR(
 prefix = file.path(
  data_dir,
  sprintf("ld_chr%d", chr)
 ),
 one_based = TRUE
)

str(csr)
head(csr$row_ptr)
tail(csr$row_ptr)
csr$nnz

chr <- 1
csr <- readLD_to_CSR(
 Glist$ldfiles[chr],
 as.integer(length(Glist$rsidsLD[[chr]])),
 as.integer(Glist$msize),
 as.numeric(0.01),
 FALSE
)
object.size(csr)
object.size(B)



# =============================================================================
# Simulation of Multi-Trait Bayesian Linear Regression (MT-BLR) Performance
# =============================================================================

# Load necessary libraries
library(qgg)
library(mr.mash.alpha)
library(ggplot2)
library(dplyr)
library(sblr)

# Create
sim <- expand.grid(
 pve = c(0.05, 0.1),
 p = c(100, 1000),
 corr = c(0.5, 1),
 n = list(c(500, 500, 500),
          c(1000, 1000, 1000),
          c(2000, 2000, 2000),
          c(500, 1000, 2000),
          c(500, 500, 1000),
          c(750, 750, 1000))
)
sim$Scenario <- paste0("sim",seq_len(nrow(sim)))
sim <- sim[, c("Scenario","pve", "p", "corr", "n")]
rownames(sim) <- sim$Scenario
print(sim)

scenario <- "sim30"
pve <- sim[scenario,"pve"]
p <- sim[scenario,"p"]
B_cor <- sim[scenario,"corr"]
n_train <- unlist(sim[scenario,"n"])


# -----------------------------------------------------------------------------
# Simulation Setup
# -----------------------------------------------------------------------------

nrep <- 10                        # Number of replicates
#p <- 1000                          # Number of predictors
#pve <- 0.05                       # Proportion of variance explained
p_causal <- 10                    # Number of causal predictors
#n_train <- c(1000, 1000, 1000)  # Observations per trait
nt <- length(n_train)            # Number of traits
n_test <- 1000                   # Number of test samples
n <- sum(n_train) + n_test       # Total sample size

# -----------------------------------------------------------------------------
# Analysis Parameters
# -----------------------------------------------------------------------------

pi <- 0.01                         # Prior inclusion probability
h2 <- 0.05                        # Heritability
nit <- 3000                      # Total MCMC iterations
nburn <- 500                     # Burn-in period

# -----------------------------------------------------------------------------
# Define Train/Test Sets
# -----------------------------------------------------------------------------

train <- split(seq_len(sum(n_train)), rep(seq_len(nt), times = n_train))
test <- (sum(n_train) + 1):n

# -----------------------------------------------------------------------------
# Output Containers
# -----------------------------------------------------------------------------

pa <- NULL          # Accuracy metrics (R², Corr) across replicates and models
pa_impr <- NULL     # Relative improvement of MT-BLR vs ST-BLR


# -----------------------------------------------------------------------------
# Simulation for testing timings
# -----------------------------------------------------------------------------

# Simulate data
dat <- simulate_mr_mash_data(
 n = n, p = p, p_causal = p_causal, r = nt,
 r_causal = list(1:nt), intercepts = rep(0, nt),
 pve = pve, B_cor = B_cor, B_scale = 1,
 w = 1, X_cor = 0.6, X_scale = 1, V_cor = 0
)

sets <- list(
 causal = dat$causal_variables,
 noncausal = setdiff(1:p, dat$causal_variables)
)

X <- dat$X
Y <- scale(dat$Y)
Yadj <- scale(Y, scale = FALSE)

# Generate hidden key
p <- ncol(X)
key <- makeKey(p)

# Summary Statistics
stat <- list(XX = list(), XXr = list(), XXre = list(),
             Xy = list(), Xyr = list(), Xyre = list(),
             yy = numeric(nt),
             n = numeric(nt), Q=list(), w=list())
for (i in 1:nt) {

 stat$XX[[i]] <- crossprod(X[train[[i]], ])
 stat$Xy[[i]] <- crossprod(X[train[[i]], ], Yadj[train[[i]], i])
 stat$yy[i]   <- sum(Yadj[train[[i]], i]^2)
 stat$n[i]    <- length(train[[i]])

 # Rotate stats
 rss <- rotateStat(XX=stat$XX[[i]], Xy=stat$Xy[[i]], stat$n[[i]])
 stat$XXr[[i]] <- rss$XX
 stat$Xyr[[i]] <- rss$Xy

 # Encode stats
 ess <- encodeStat(stat$XX[[i]], stat$Xy[[i]], key)

 # Rotate encoded stats
 ress <- rotateStat(XX=ess$XX, Xy=ess$Xy, stat$n[[i]])

 stat$XXre[[i]] <- ress$XX
 stat$Xyre[[i]] <- ress$Xy


}



system.time(MT1 <- sblr(yy = stat$yy, Xy = stat$Xy, XX = stat$XX, n = stat$n,
                        sets = NULL, model = NULL, algorithm="default",
                        h2 = h2, pi = pi, updateB =TRUE, updatePi = TRUE,
                        method = "bayesc", nit = nit, nburn = nburn, verbose = FALSE))

# -----------------------------------------------------------------------------
# Simulation Loop
# -----------------------------------------------------------------------------

for (rep in 1:nrep) {
 cat("Running replicate", rep, "\n")

 # Simulate data
 dat <- simulate_mr_mash_data(
  n = n, p = p, p_causal = p_causal, r = nt,
  r_causal = list(1:nt), intercepts = rep(0, nt),
  pve = pve, B_cor = B_cor, B_scale = 1,
  w = 1, X_cor = 0.6, X_scale = 1, V_cor = 0
 )

 sets <- list(
  causal = dat$causal_variables,
  noncausal = setdiff(1:p, dat$causal_variables)
 )

 X <- dat$X
 Y <- scale(dat$Y)
 Yadj <- scale(Y, scale = FALSE)

 # Generate hidden key
 p <- ncol(X)
 key <- makeKey(p)

 # Summary Statistics
 stat <- list(XX = list(), XXr = list(), XXre = list(),
              Xy = list(), Xyr = list(), Xyre = list(),
              yy = numeric(nt),
              n = numeric(nt), Q=list(), w=list())
 for (i in 1:nt) {

  stat$XX[[i]] <- crossprod(X[train[[i]], ])
  stat$Xy[[i]] <- crossprod(X[train[[i]], ], Yadj[train[[i]], i])
  stat$yy[i]   <- sum(Yadj[train[[i]], i]^2)
  stat$n[i]    <- length(train[[i]])

  # Rotate stats
  rss <- rotateStat(XX=stat$XX[[i]], Xy=stat$Xy[[i]], stat$n[[i]])
  stat$XXr[[i]] <- rss$XX
  stat$Xyr[[i]] <- rss$Xy

  # Encode stats
  ess <- encodeStat(stat$XX[[i]], stat$Xy[[i]], key)

  # Rotate encoded stats
  ress <- rotateStat(XX=ess$XX, Xy=ess$Xy, stat$n[[i]])

  stat$XXre[[i]] <- ress$XX
  stat$Xyre[[i]] <- ress$Xy


 }

 # ---------------------------------------------------------------------------
 # Fit Single-Trait BLR Models
 # ---------------------------------------------------------------------------

 yobs <- Y[test, ]
 ypred <- NULL

 for (i in 1:nt) {
  fitST <- qgg:::blr(yy = stat$yy[i], Xy = stat$Xy[[i]],
                     XX = stat$XX[[i]], n = stat$n[i],
                     h2 = h2, pi = pi,
                     method = "bayesc", nit = nit, nburn = nburn)
  ypred <- cbind(ypred, crossprod(t(X[test, ]), fitST$bm))
 }

 rownames(ypred) <- rownames(yobs) <- seq_len(nrow(yobs))

 st_pa_rep <- do.call(rbind, lapply(1:nt, function(x) {
  pa_rep <- qgg:::acc(ypred = ypred[, x], yobs = yobs[, x])
  data.frame(Model = "ST", Trait = paste0("T", x), Replicate = rep, pa_rep)
 }))
 st_r2 <- st_pa_rep$R2
 pa <- rbind(pa, st_pa_rep)

 # ---------------------------------------------------------------------------
 # Fit Multi-Trait BLR Models
 # ---------------------------------------------------------------------------

 mt_models <- list(
  MT1 = sblr(yy = stat$yy, Xy = stat$Xy, XX = stat$XX, n = stat$n,
              sets = NULL, model = NULL, algorithm="default",
              h2 = h2, pi = pi, updateB =TRUE, updatePi = TRUE,
              method = "bayesc", nit = nit, nburn = nburn, verbose = FALSE)

 )

 for (model_name in names(mt_models)) {

  bm <- mt_models[[model_name]]$bm
  ypred <- crossprod(t(X[test, ]), bm)
  rownames(ypred) <- rownames(yobs)

  mt_pa_rep <- do.call(rbind, lapply(1:nt, function(x) {
   pa_rep <- qgg:::acc(ypred = ypred[, x], yobs = yobs[, x])
   data.frame(Model = model_name, Trait = paste0("T", x), Replicate = rep, pa_rep)
  }))
  pa <- rbind(pa, mt_pa_rep)

  # Relative Improvement
  mt_r2 <- mt_pa_rep$R2
  impr <- 100 * (mt_r2 - st_r2) / st_r2
  pa_impr <- rbind(pa_impr, data.frame(
   Model = model_name, Trait = paste0("T", 1:nt), Replicate = rep,
   mt_r2 = mt_r2, st_r2 = st_r2, improvement = impr
  ))
 }
}
rownames(pa) <- NULL
rownames(pa_impr) <- NULL





# -----------------------------------------------------------------------------
# Summarize and Plot Results
# -----------------------------------------------------------------------------



# Summarize raw prediction accuracy (R²) across replicates
r2_summary <- pa %>%
 group_by(Model, Trait) %>%
 summarise(
  mean_r2 = mean(R2),                                # Average R² across replicates
  sem = sd(R2) / sqrt(n()),                          # Standard error of the mean
  .groups = "drop"
 )

# Plot R² estimates as points ± SEM
ggplot(r2_summary, aes(x = Trait, y = mean_r2, color = Model)) +
 geom_point(position = position_dodge(width = 0.5), size = 3) +  # Mean R² points
 geom_errorbar(aes(ymin = mean_r2 - sem, ymax = mean_r2 + sem), # SEM error bars
               position = position_dodge(width = 0.5),
               width = 0.2) +
 labs(
  title = "Mean Prediction Accuracy (R²) by Model",     # Title
  y = "Mean R² ± SEM",                                  # Y-axis
  x = "Trait"                                            # X-axis
 ) +
 theme_minimal() +                                       # Clean look
 theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14), # Centered, bold
  axis.title = element_text(size = 12),                           # Label font
  axis.text = element_text(size = 10),                            # Tick font
  legend.position = "top"                                         # Legend on top
 )

# Define a common dodge
pd <- position_dodge(width = 0.6)

# Summarize raw prediction accuracy (R²) across replicates
r2_summary <- pa %>%
  group_by(Model, Trait) %>%
  summarise(
    mean_r2 = mean(R2),                                # Average R² across replicates
    sem = sd(R2) / sqrt(n()),                          # Standard error of the mean
    .groups = "drop"
  )

# Rename Model levels
r2_summary$Model <- factor(
  r2_summary$Model,
  levels = c("MT1","MT2","ST"),
  labels = c("MT SS", "MT rotated SS", "ST SS")
)
r2_summary$Trait <- factor(
 r2_summary$Trait,
 levels = c("T1","T2","T3"),
 labels = c("D1 (n=500)", "D2 (n=1000)", "D3 (n=2000)")
)

# Plot R² estimates as points ± SEM
ggplot(r2_summary, aes(x = Trait, y = mean_r2, color = Model)) +
  geom_point(position = pd, size = 3) +  # Mean R² points
  geom_errorbar(
    aes(ymin = mean_r2 - sem, ymax = mean_r2 + sem),
    position = pd,
    width = 0.25
  ) +
  labs(
    title = "Mean Prediction Accuracy (R²) by Model",     # Title
    y = "Mean R² ± SEM",                                  # Y-axis
    x = "Trait"                                           # X-axis
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), # Centered, bold
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

