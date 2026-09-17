# =============================================================================
# 02_annotate.R  --  DiLeucine transcriptomics (Fall 2026), Stage 02
#
# One job: give every transcript cluster ID a human-readable name.
#
# The matrix rownames are Affymetrix transcript cluster IDs (TC0900008219.hg.1).
# They are internal manufacturing identifiers -- they say nothing about biology.
# Without a symbol table you can compute a p-value but you cannot read a result,
# and you certainly cannot read a PCA loading.
#
# Source: the Bioconductor annotation package for this array. Chosen over the
# TAC comparison exports because it is versioned and citable -- you can state
# in the methods which annotation build named your genes, which you cannot do
# for an ad-hoc software export.
#
# In : 01_Preprocessing/02_data/01_expression.rds   (rownames only)
# Out: 01_Preprocessing/02_data/02_annotation.csv
#      01_Preprocessing/02_data/sessionInfo_02.txt
#
# First run needs:
#   install.packages("BiocManager")
#   BiocManager::install("clariomshumanhttranscriptcluster.db")
# =============================================================================

library(here)

DATA <- here("01_Preprocessing", "02_data")
PKG  <- "clariomshumanhttranscriptcluster.db"

if (!requireNamespace(PKG, quietly = TRUE))
  stop("Annotation package missing. Run:\n",
       '  install.packages("BiocManager")\n',
       '  BiocManager::install("', PKG, '")\n',
       "If the IDs do not map (see the coverage check below), try the non-HT\n",
       "build instead: clariomshumantranscriptcluster.db", call. = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(PKG, character.only = TRUE)
})
db <- get(PKG)

d   <- readRDS(file.path(DATA, "01_expression.rds"))
ids <- rownames(d$expr)
cat(sprintf("transcript clusters to annotate: %d\n", length(ids)))


# ---- 1. map ------------------------------------------------------------------
# A transcript cluster can overlap more than one gene, so a key can return
# several values. multiVals = "first" takes one and keeps the table rectangular
# (one row per cluster, which is what every downstream join assumes). How often
# that happened is counted below rather than hidden.

pull <- function(column)
  unname(suppressMessages(
    AnnotationDbi::mapIds(db, keys = ids, keytype = "PROBEID",
                          column = column, multiVals = "first")))

ann <- data.frame(
  ID        = ids,
  symbol    = pull("SYMBOL"),
  gene_name = pull("GENENAME"),
  entrez    = pull("ENTREZID"),
  stringsAsFactors = FALSE
)
stopifnot(identical(ann$ID, ids), nrow(ann) == length(ids))


# ---- 2. how good is the mapping? --------------------------------------------
mapped <- !is.na(ann$symbol) & ann$symbol != ""
cat(sprintf("mapped to a gene symbol: %d of %d (%.1f%%)\n",
            sum(mapped), nrow(ann), 100 * mean(mapped)))

# Sanity floor. Clariom S is a curated, well-annotated array; a large majority
# of clusters should resolve. A very low rate does not mean "poorly annotated
# array", it means the keytype or the package build is wrong -- stop rather
# than carry an empty column through four more stages.
if (mean(mapped) < 0.50)
  stop(sprintf("only %.1f%% of IDs mapped. Wrong package build or keytype -- ",
               100 * mean(mapped)),
       "try clariomshumantranscriptcluster.db before going further.", call. = FALSE)

multi <- suppressMessages(
  AnnotationDbi::mapIds(db, keys = ids, keytype = "PROBEID",
                        column = "SYMBOL", multiVals = "list"))
n_multi <- sum(lengths(multi) > 1)
cat(sprintf("clusters matching more than one symbol: %d (%.1f%%) -- first kept\n",
            n_multi, 100 * n_multi / nrow(ann)))

dup <- sum(duplicated(ann$symbol[mapped]))
cat(sprintf("symbols appearing on more than one cluster: %d\n", dup))
# Not an error. Several clusters can measure the same gene. It matters later:
# a gene set test will see that gene more than once unless clusters are
# collapsed, which is a Stage 04 decision, not one to make silently here.


# ---- 3. save -----------------------------------------------------------------
write.csv(ann, file.path(DATA, "02_annotation.csv"), row.names = FALSE)
cat(sprintf("wrote 02_annotation.csv (%d rows)\n", nrow(ann)))

writeLines(c(sprintf("annotation package: %s %s", PKG,
                     as.character(utils::packageVersion(PKG))),
             capture.output(sessionInfo())),
           file.path(DATA, "sessionInfo_02.txt"))

cat("\nStage 02 complete.\n")
