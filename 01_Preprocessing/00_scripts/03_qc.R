# =============================================================================
# 03_qc.R  --  DiLeucine transcriptomics (Fall 2026), Stage 03
#
# One job: decide whether to believe this data, and find out what is actually
# driving the variation in it, BEFORE any hypothesis is tested.
#
# Nothing here is a result. Everything here is either an audit (did the data
# arrive intact?) or a description (where does the variance live?). The output
# is a decision -- proceed, or go back -- plus the numbers that justify the
# model choices made at Stage 04.
#
# In : 01_Preprocessing/02_data/01_expression.rds
#      01_Preprocessing/02_data/02_annotation.csv
# Out: 01_Preprocessing/01_reports/03_qc_figures.pdf
#      01_Preprocessing/01_reports/03_qc_summary.txt
#      01_Preprocessing/01_reports/03_pca_scores.csv
#      01_Preprocessing/01_reports/03_pca_loadings.csv
# =============================================================================

library(here)
set.seed(1)

DATA   <- here("01_Preprocessing", "02_data")
REPORT <- here("01_Preprocessing", "01_reports")
dir.create(REPORT, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 500   # permutations for the subject-variance null
N_LOAD <- 20    # genes listed per loading tail


# ---- 0. logging + load -------------------------------------------------------
.log <- character()
say <- function(fmt, ...) { m <- sprintf(fmt, ...); .log <<- c(.log, m); cat(m, "\n", sep=""); invisible(m) }
check <- function(label, ok) {
  ok <- isTRUE(ok); say("[%s] %s", if (ok) "PASS" else "FAIL", label)
  if (!ok) stop("failed check: ", label, call. = FALSE); invisible(TRUE)
}

d <- readRDS(file.path(DATA, "01_expression.rds"))
expr <- d$expr; targets <- d$targets

ann <- read.csv(file.path(DATA, "02_annotation.csv"), stringsAsFactors = FALSE)
ann <- ann[match(rownames(expr), ann$ID), ]
check("annotation is aligned to the matrix", identical(ann$ID, rownames(expr)))

# Label = symbol where we have one, transcript cluster ID where we do not.
# Never silently drop an unnamed cluster: an unnamed gene is still a gene.
gene_label <- ifelse(is.na(ann$symbol) | ann$symbol == "", ann$ID, ann$symbol)
names(gene_label) <- ann$ID

check("matrix and sheet still aligned", identical(colnames(expr), as.character(targets$sample_id)))
say("QC run: %s", format(Sys.time()))
say("%d transcript clusters x %d arrays", nrow(expr), ncol(expr))
say("")

pdf(file.path(REPORT, "03_qc_figures.pdf"), width = 11, height = 8)
on.exit(dev.off(), add = TRUE)


# =============================================================================
# 1. ARRAY-LEVEL DISTRIBUTIONS  --  an audit of the normalisation
#
# SST-RMA quantile-normalises at PROBE level, then applies the signal-space
# transformation and summarises to gene level. That ordering matters for how
# you read these plots: array MEANS come out tightly matched, but medians and
# IQRs need not, because summarisation reshapes the distribution AFTER
# normalisation. So this is an audit of centring, not a guarantee of identical
# shape -- and not a reliable way to find a bad array either, since whatever
# would have flagged one was largely absorbed by the normalisation.
# =============================================================================
med <- apply(expr, 2, median)
iqr <- apply(expr, 2, IQR)

par(mfrow = c(1,1), mar = c(8,4,3,1))
boxplot(expr, las = 2, cex.axis = .4, outline = FALSE, col = "grey90",
        main = "log2 signal by array (must be near-identical after quantile normalisation)")
plot(med, iqr, pch = 19, cex = .7, xlab = "array median", ylab = "array IQR",
     main = "centre vs spread, one point per array")
text(med, iqr, labels = ifelse(abs(scale(med)) > 2.5 | abs(scale(iqr)) > 2.5,
                               targets$sample_id, ""), pos = 3, cex = .5)

say("--- 1. array distributions ---")
mu <- colMeans(expr)
say("array means  : %.3f to %.3f (spread %.3f)  <- what normalisation matched",
    min(mu), max(mu), diff(range(mu)))
say("array medians: %.3f to %.3f (spread %.3f)", min(med), max(med), diff(range(med)))
say("array IQRs   : %.3f to %.3f (spread %.3f)", min(iqr), max(iqr), diff(range(iqr)))
say("mean IQR by timepoint: %s",
    paste(sprintf("%s %.3f", levels(targets$time), tapply(iqr, targets$time, mean)), collapse = "  "))

# Residual spread SHOULD be tiny. If IQR still tracks timepoint after
# normalisation, that is real biology leaking through (a huge coordinated
# response widens the distribution), not a technical failure.
a <- summary(aov(iqr ~ targets$time))[[1]]
say("IQR by timepoint: F = %.2f, p = %.4g", a[1,"F value"], a[1,"Pr(>F)"])
boxplot(iqr ~ targets$time, main = "array IQR by timepoint", xlab = "", ylab = "IQR")
say("")


# =============================================================================
# 2. SAMPLE-SAMPLE CORRELATION  --  an identity check
#
# Between-subject variation in the muscle transcriptome is large and stable
# (fibre-type proportion, training history, genotype). So a subject's own six
# biopsies should resemble each other more than anyone else's. If they do not,
# a tube was mislabelled -- the one error that is invisible to every downstream
# test and fatal to all of them.
#
# Spearman, not Pearson: rank-based, so a handful of extreme probes cannot
# manufacture a correlation.
# =============================================================================
cm  <- cor(expr, method = "spearman")
ord <- order(targets$subject, targets$period, targets$time)

par(mar = c(2,2,3,1))
image(cm[ord, ord], axes = FALSE, col = hcl.colors(64, "YlOrRd", rev = TRUE),
      main = "sample-sample Spearman correlation, ordered by subject\n(expect bright 6x6 blocks on the diagonal)")

hc <- hclust(as.dist(1 - cm), method = "average")
par(mar = c(9,4,3,1))
plot(hc, labels = targets$sample_id, cex = .35, main = "average-linkage clustering (1 - Spearman r)",
     xlab = "", sub = "")

say("--- 2. sample identity ---")
subj <- as.character(targets$subject)
same_subj <- outer(subj, subj, "==")
ut <- upper.tri(cm)

# The naive comparison, and why it is the WRONG one here. The exercise response
# is so large that a subject's own PRE and 24hr arrays are less alike than two
# different subjects' 24hr arrays. Pooled over timepoints, this statistic
# measures the exercise response, not sample identity -- it can come out
# backwards on perfectly clean data.
say("pooled over all timepoints (CONFOUNDED -- see below):")
say("  within-subject  r = %.4f", mean(cm[same_subj & ut]))
say("  between-subject r = %.4f", mean(cm[!same_subj & ut]))

# The right comparison: like with like. Hold timepoint fixed, so the only thing
# separating a within-subject pair from a between-subject pair is whose muscle
# it came from.
say("stratified by timepoint (the valid comparison):")
for (tp in levels(targets$time)) {
  k <- targets$time == tp
  sub_cm <- cm[k, k]; ss <- outer(subj[k], subj[k], "=="); u <- upper.tri(sub_cm)
  say("  %-5s (%2d arrays)  within-subject r = %.4f   between-subject r = %.4f   diff %+.4f",
      tp, sum(k), mean(sub_cm[ss & u]), mean(sub_cm[!ss & u]),
      mean(sub_cm[ss & u]) - mean(sub_cm[!ss & u]))
}

# Nearest neighbour, also restricted to the same timepoint. Within a timepoint a
# subject has at most one other array (the other condition), so hitting it by
# chance is roughly 1 in 30 -- a demanding test, and a specific one: it names
# the arrays that fail rather than reporting an average.
cm_nn <- cm
cm_nn[outer(as.character(targets$time), as.character(targets$time), "!=")] <- -Inf
diag(cm_nn) <- -Inf
nn      <- colnames(cm_nn)[apply(cm_nn, 2, which.max)]
nn_subj <- targets$subject[match(nn, targets$sample_id)]
hit     <- nn_subj == targets$subject
say("within timepoint, nearest neighbour is the same subject: %d of %d (%.0f%%; chance ~3%%)",
    sum(hit), ncol(expr), 100*mean(hit))

top5 <- apply(cm_nn, 2, function(x) head(order(x, decreasing = TRUE), 5))
in5  <- vapply(seq_len(ncol(expr)), function(j)
          any(targets$subject[top5[, j]] == targets$subject[j]), logical(1))
say("same subject in the top 5 same-timepoint neighbours: %d of %d (%.0f%%)",
    sum(in5), ncol(expr), 100*mean(in5))
say("")


# =============================================================================
# 3. PCA  --  what is the largest coordinated pattern in the data?
#
# 95 samples described by 21,448 numbers each = a cloud of 95 points in
# 21,448-dimensional space. PCA rebuilds the axes: PC1 is the direction of
# greatest spread; PC2 the direction of greatest REMAINING spread that is
# perpendicular to PC1; and so on. Perpendicular means their variances do not
# overlap, so the percentages add.
#
# t(expr): prcomp wants samples as rows. Our matrix is genes x samples.
# scale. = FALSE: all values are already log2 on one common scale, and we WANT
#   high-variance genes to dominate -- those are the responding genes. Scaling
#   would give 21,000 flat low-expressed probes equal weight and bury the signal.
# Centring (the default) is kept: without it PC1 is just "average brightness".
#
# PCA is told NOTHING about the design. If it recovers the timepoints anyway,
# that is an unsupervised method independently confirming the experiment worked.
# =============================================================================
pca <- prcomp(t(expr), scale. = FALSE)
pv  <- 100 * pca$sdev^2 / sum(pca$sdev^2)

check("PCA ran on samples, not genes", nrow(pca$x) == ncol(expr))
check("PCA rows are the sample ids", identical(rownames(pca$x), as.character(targets$sample_id)))

say("--- 3. PCA ---")
say("variance explained: %s", paste(sprintf("PC%d %.1f%%", 1:6, pv[1:6]), collapse = "  "))

par(mfrow = c(1,1), mar = c(4,4,3,1))
barplot(pv[1:10], names.arg = paste0("PC", 1:10), col = "grey80",
        ylab = "% variance", main = "scree")

par(mfrow = c(2,3), mar = c(4,4,3,1))
for (v in c("time","condition","subject","period","sequence")) {
  f <- factor(targets[[v]])
  plot(pca$x[,1], pca$x[,2], col = as.integer(f), pch = 19, cex = .9,
       xlab = sprintf("PC1 (%.1f%%)", pv[1]), ylab = sprintf("PC2 (%.1f%%)", pv[2]),
       main = paste("coloured by", v))
  if (nlevels(f) <= 8) legend("topright", levels(f), col = seq_len(nlevels(f)), pch = 19, cex = .7)
}
par(mfrow = c(1,1))

# Loadings name the program. A PC is a weighted recipe over genes; the weights
# ARE the loadings. Reading the extremes tells you which biology this axis is.
tail_of <- function(pc, dec) {
  x <- head(sort(pca$rotation[, pc], decreasing = dec), N_LOAD)
  data.frame(PC = pc, end = if (dec) "positive" else "negative",
             ID = names(x), symbol = gene_label[names(x)],
             loading = round(unname(x), 4), row.names = NULL)
}
loadings <- do.call(rbind, lapply(1:3, function(i)
  rbind(tail_of(i, TRUE), tail_of(i, FALSE))))
write.csv(loadings, file.path(REPORT, "03_pca_loadings.csv"), row.names = FALSE)

for (i in 1:2) {
  say("PC%d positive end: %s", i,
      paste(head(loadings$symbol[loadings$PC==i & loadings$end=="positive"], 12), collapse=", "))
  say("PC%d negative end: %s", i,
      paste(head(loadings$symbol[loadings$PC==i & loadings$end=="negative"], 12), collapse=", "))
}
say("")


# =============================================================================
# 4. VARIANCE PARTITIONING  --  how much of each axis does each design factor own?
#
# This is the number that justifies the crossover. Fit PC scores on one factor
# at a time and read R^2: the share of that axis the factor accounts for.
#
# R^2 is inflated by degrees of freedom -- subject has 16 levels (15 df) and can
# fit noise simply by having enough parameters. So adjusted R^2 is reported
# beside it, and PC1-by-subject is checked against a permutation null: shuffle
# the subject labels, refit, repeat. That null is what "explaining nothing"
# actually looks like for a factor with this many levels.
# =============================================================================
vars <- c("time","condition","subject","period","sequence"); NPC <- 5

r2 <- sapply(vars, function(v) sapply(1:NPC, function(i) {
  a <- summary(aov(pca$x[, i] ~ targets[[v]]))[[1]]; a[1,"Sum Sq"] / sum(a[,"Sum Sq"]) }))
adj <- sapply(vars, function(v) sapply(1:NPC, function(i)
  summary(lm(pca$x[, i] ~ targets[[v]]))$adj.r.squared))
rownames(r2) <- rownames(adj) <- paste0("PC", 1:NPC)

say("--- 4. variance partitioning ---")
say("R^2 (%%) of each PC explained by each factor:")
say("%s", paste(capture.output(print(round(100*r2, 1))), collapse = "\n"))
say("adjusted R^2 (%%), same cells -- the honest version:")
say("%s", paste(capture.output(print(round(100*adj, 1))), collapse = "\n"))

obs  <- r2["PC1","subject"]
null <- replicate(N_PERM, { a <- summary(aov(pca$x[,1] ~ sample(targets$subject)))[[1]]
                            a[1,"Sum Sq"]/sum(a[,"Sum Sq"]) })
say("PC1 by subject: observed R^2 = %.3f; permutation null 95%% range %.3f - %.3f; percentile %.1f%%",
    obs, quantile(null, .025), quantile(null, .975), 100*mean(null < obs))
hist(null, breaks = 40, col = "grey85", xlim = range(c(null, obs)),
     main = "PC1 variance explained by subject: permutation null", xlab = "R^2")
abline(v = obs, col = "#D55E00", lwd = 2)
say("")


# =============================================================================
# 5. RESIDUAL PCA  --  what is left once the designed effect is removed?
#
# Strip condition x time, then look again. Whatever structure survives is
# nuisance: subject, plate position, or something unmodelled. Large structure
# here that is NOT subject is a confounder you have not accounted for.
# =============================================================================
res   <- residuals(lm(t(expr) ~ targets$grp))
pca_r <- prcomp(res)
pvr   <- 100 * pca_r$sdev^2 / sum(pca_r$sdev^2)

say("--- 5. residual PCA (condition x time removed) ---")
say("variance explained: %s", paste(sprintf("rPC%d %.1f%%", 1:5, pvr[1:5]), collapse = "  "))
for (v in vars)
  say("  rPC1-3 adj R^2 by %-9s : %s", v,
      paste(sprintf("%5.1f%%", 100*sapply(1:3, function(i)
        summary(lm(pca_r$x[, i] ~ targets[[v]]))$adj.r.squared)), collapse = " "))
par(mfrow = c(1,2))
for (v in c("subject","period")) {
  f <- factor(targets[[v]])
  plot(pca_r$x[,1], pca_r$x[,2], col = as.integer(f), pch = 19, cex = .9,
       xlab = sprintf("rPC1 (%.1f%%)", pvr[1]), ylab = sprintf("rPC2 (%.1f%%)", pvr[2]),
       main = paste("residual PCA, coloured by", v))
}
par(mfrow = c(1,1))
say("")


# =============================================================================
# 6. PLATE GEOMETRY  --  the confound that cannot be modelled away
#
# All 96 arrays sat on ONE plate, so there is no batch effect to worry about.
# But each subject occupies six consecutive wells with T1 in the first three and
# T2 in the last three, identically for every subject. Plate position is
# therefore ALIASED with period: a spatial gradient across the plate and a true
# period effect would produce the same numbers, and no model can tell them apart.
#
# Condition is NOT aliased, because sequence varies between subjects (10 A-first,
# 6 B-first). That is what protects the treatment contrast.
#
# We cannot separate position from period. We CAN ask whether a gradient exists
# at all -- if none does, the aliasing costs nothing.
# =============================================================================
plate_row <- match(substr(targets$well, 1, 1), LETTERS)     # A-H -> 1-8
plate_col <- as.integer(substr(targets$well, 2, 3))         # 01-12

say("--- 6. plate geometry ---")
say("one plate, barcode %s; rows A-H, columns 1-12", unique(targets$barcode))
for (nm in c("row","col")) {
  pos <- if (nm == "row") plate_row else plate_col
  s <- summary(lm(med ~ pos))
  say("array median vs plate %s: slope %+.4f log2/step, R^2 = %.3f, p = %.3g",
      nm, coef(s)[2,1], s$r.squared, coef(s)[2,4])
  for (i in 1:3) {
    si <- summary(lm(pca$x[, i] ~ pos))
    say("  PC%d vs plate %-3s : R^2 = %.3f, p = %.3g", i, nm, si$r.squared, coef(si)[2,4])
  }
}
par(mfrow = c(1,2))
boxplot(med ~ plate_row, main = "array median by plate row", xlab = "row (A-H)", ylab = "median log2")
boxplot(pca$x[,1] ~ plate_col, main = "PC1 by plate column", xlab = "column", ylab = "PC1")
par(mfrow = c(1,1))
say("")


# =============================================================================
# 7. THE TWO NAMED ODDITIES  --  carried over unresolved from the Summer analysis
# =============================================================================
say("--- 7. flagged samples and probes ---")

# (a) S07_T1_PRE_A: an immediate-early signature in a RESTING biopsy. If real,
#     that biopsy is not a baseline. Leading hypothesis is excision-to-freeze
#     time -- handling stress induces these genes ex vivo within minutes.
ieg <- c("NR4A3","EGR1","NR4A1","ATF3","JUNB","FOS","MYC")
ieg_id <- names(gene_label)[gene_label %in% ieg]
if (length(ieg_id)) {
  z <- t(scale(t(expr[ieg_id, , drop = FALSE])))          # per-gene z across arrays
  ieg_score <- colMeans(z, na.rm = TRUE)
  o <- order(ieg_score, decreasing = TRUE)
  say("immediate-early score, top 8 arrays (mean z over %s):", paste(ieg, collapse="/"))
  for (j in head(o, 8))
    say("  %-16s %-5s z = %+.2f", targets$sample_id[j], as.character(targets$time[j]), ieg_score[j])
  pre <- targets$time == "Pre"
  say("rank of S07_T1_PRE_A among the %d PRE arrays: %d",
      sum(pre), which(targets$sample_id[pre][order(ieg_score[pre], decreasing=TRUE)] == "S07_T1_PRE_A"))
  boxplot(ieg_score ~ targets$time, main = "immediate-early score by timepoint",
          xlab = "", ylab = "mean z")
  points(as.integer(targets$time)[targets$sample_id=="S07_T1_PRE_A"],
         ieg_score[targets$sample_id=="S07_T1_PRE_A"], col = "#D55E00", pch = 19, cex = 1.4)
} else say("  (immediate-early genes not found in the annotation -- skipped)")

# (b) TC0900008219.hg.1: swings ~5 -> ~16 -> ~7 log2 at Pre/3hr/24hr in EVERY
#     subject. An 11-log2 (2,000-fold) induction is not plausible biology; a
#     probe behaving like this is usually cross-hybridising or a control.
probe <- "TC0900008219.hg.1"
if (probe %in% rownames(expr)) {
  say("%s (%s) mean log2 by timepoint: %s", probe, gene_label[probe],
      paste(sprintf("%s %.2f", levels(targets$time),
                    tapply(expr[probe, ], targets$time, mean)), collapse = "  "))
  boxplot(expr[probe, ] ~ targets$time, main = paste(probe, "-", gene_label[probe]),
          xlab = "", ylab = "log2 signal")
}

rng <- apply(expr, 1, function(x) diff(range(x)))
say("transcripts with range > 10 log2 units: %d", sum(rng > 10))
say("their labels: %s", paste(head(gene_label[order(rng, decreasing=TRUE)], 15), collapse=", "))
say("")


# ---- 8. save -----------------------------------------------------------------
scores <- data.frame(targets[, c("sample_id","subject","period","time","condition","well")],
                     round(pca$x[, 1:6], 3), row.names = NULL)
write.csv(scores, file.path(REPORT, "03_pca_scores.csv"), row.names = FALSE)
writeLines(.log, file.path(REPORT, "03_qc_summary.txt"))

cat("\nStage 03 complete. Figures, summary, scores and loadings are in 01_reports/.\n")
