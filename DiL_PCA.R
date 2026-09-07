# Sample Identity 
Correlation_Matrix <- cor(expr, method = "spearman")
order <- order(targets$subject, targets$period, targets$time)
image(Correlation_Matrix[order, order], axes = FALSE, main = "sample-sample correlation, ordered by subject")

hc <- hclust(as.dist(1 - Correlation_Matrix), method = "average")
plot(hc, labels = targets$array, cex = 0.35, main = "", xlab = "", sub = "")


# What is driving the variability?
# PCA: Principle Component Analysis
pca <- prcomp(t(expr), scale. = FALSE) # These are already log2 and on a common scale
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

for (v in c("time","condition","subject","period","sequence")) {
  f <- factor(targets[[v]])
  plot(pca$x[,1], pca$x[,2], col = as.integer(f), pch = 19,
       xlab = paste0("PC1 (", pv[1], "%)"), ylab = paste0("PC2 (", pv[2], "%)"),
       main = paste("coloured by", v))
  legend("topright", levels(f), col = seq_along(levels(f)), pch = 19, cex = .6, ncol = 2)
}

# variance in each PC attributable to subject
sapply(1:5, function(i) {
  a <- summary(aov(pca$x[,i] ~ targets$subject))[[1]]
  round(a[1,2] / sum(a[,2]), 3)
})

# CRITICAL: assert alignment before you trust any of this
stopifnot(identical(rownames(pca$x), targets$sample_id))

vars <- c("time", "condition", "subject", "period", "sequence")
npc  <- 5

r2 <- sapply(vars, function(v) {
  sapply(seq_len(npc), function(i) {
    a <- summary(aov(pca$x[, i] ~ targets[[v]]))[[1]]
    a[1, "Sum Sq"] / sum(a[, "Sum Sq"])
  })
})
rownames(r2) <- paste0("PC", seq_len(npc))
round(100 * r2, 1)

adj_r2 <- sapply(vars, function(v) {
  sapply(seq_len(npc), function(i) {
    summary(lm(pca$x[, i] ~ targets[[v]]))$adj.r.squared
  })
})
rownames(adj_r2) <- paste0("PC", seq_len(npc))
round(100 * adj_r2, 1)

identical(rownames(pca$x), targets$sample_id)   # must be TRUE
length(pca$x[, 1])                              # must be 95, not 21448
table(table(targets$subject))          # must be: 5 appearing once, 6 appearing 15 times


lab(pca$rotation[, 1])           # PC1 positive end
lab(pca$rotation[, 2], dec = FALSE)   # PC2 negative end

out <- order(pca$x[, 2])[1:2]
targets[out, c("sample_id", "well", "subject", "period", "time", "condition")] # These are the two dots at the both of the PCA plots

