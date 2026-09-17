# =============================================================================
# 01_input.R  --  DiLeucine transcriptomics Stage 01
#
# One job: turn two inputs into one aligned, asserted, correctly typed object

# In : 00_inputs/Raw_data.txt   21,448 transcript clusters x 95 arrays,
#                               log2 SST-RMA signal exported from TAC 4.0.
#                               Normalisation happened in TAC.
#      00_inputs/targets.csv    curated sample sheet. Built once from the
#                               Affymetrix .ARR UserAttributes (Subject,
#                               Time Point, Condition, Group) and checked
#                               against the array names. It is an INPUT, not
#                               a derivation -- who a sample is should not
#                               live inside a regex.
#
# Out: 01_Preprocessing/02_data/01_expression.rds        expr + targets
#      01_Preprocessing/02_data/01_alignment_checks.txt  the assertion log
#      01_Preprocessing/02_data/sessionInfo_01.txt
#
# Nothing in 00_inputs is ever written to.
# =============================================================================

library(here)

IN   <- here("00_inputs")
DATA <- here("01_Preprocessing", "02_data")
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)


# ---- 0. assertion machinery -------------------------------------------------
# Every structural claim this script makes gets written to a file. A pipeline
# that only prints its checks to the console has not actually recorded them.

.log <- character()
say <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  .log <<- c(.log, msg)
  cat(msg, "\n", sep = "")
  invisible(msg)
}
check <- function(label, ok) {
  ok <- isTRUE(ok)
  say("[%s] %s", if (ok) "PASS" else "FAIL", label)
  if (!ok) stop("failed check: ", label, call. = FALSE)
  invisible(TRUE)
}

say("Stage 01 run: %s", format(Sys.time()))
say("")


# ---- 1. the expression matrix ----------------------------------------------
# check.names = FALSE keeps the TAC column names byte-exact, which is what
# targets$array is matched on. Let R "fix" them and the join silently breaks.

raw <- read.delim(file.path(IN, "Raw_data.txt"), row.names = 1, check.names = FALSE)

# as.matrix() will coerce the whole thing to character if ONE column is not
# numeric, and every downstream number would then be wrong but not obviously so.
check("every column of Raw_data.txt is numeric", all(vapply(raw, is.numeric, logical(1))))

expr <- as.matrix(raw)
rm(raw) # remove raw
say("expression matrix: %d transcript clusters x %d arrays", nrow(expr), ncol(expr))

check("matrix fingerprint matches the Summer 2026 ingest (12843602.47)",
      abs(sum(expr) - 12843602.47) < 0.01)

check("no missing values", !anyNA(expr))
say("value range: %.2f to %.2f, mean %.2f", min(expr), max(expr), mean(expr))
say("")


# ---- 2. the curated sample sheet -------------------------------------------
# barcode is a 22-digit plate number. Left to itself R reads it as a double and
# silently renders it 5.50738454734002e+21 -- an identifier is a label, not a
# quantity, so it is read as text.
targets <- read.csv(file.path(IN, "targets.csv"), stringsAsFactors = FALSE,
                    colClasses = c(barcode = "character"))

check("targets.csv has 95 rows",        nrow(targets) == 95)
check("sample_id is unique",            !any(duplicated(targets$sample_id)))
check("array is unique",                !any(duplicated(targets$array)))
check("no missing sample metadata",
      !anyNA(targets[, c("subject", "period", "time", "condition", "sequence")]))
say("")


# ---- 3. cross-check the sheet against the array names ----------------------
# targets.csv came from the .ARR attributes. The array names are a completely
# independent record of the same facts. Re-derive from the names and demand
# agreement: two sources that were never copied from each other.

pat <- "^[0-9]+_([A-H][0-9]{2})_(S[0-9]+)(T[12])(Pre|3hr|24hr)([AB])\\.sst-rma"
m   <- regmatches(targets$array, regexec(pat, targets$array))
check("every array name parses", all(lengths(m) == 6))

from_name <- data.frame(
  well      = vapply(m, `[`, "", 2),
  subject   = vapply(m, `[`, "", 3),
  period    = vapply(m, `[`, "", 4),
  time      = vapply(m, `[`, "", 5),
  condition = vapply(m, `[`, "", 6)
)

for (v in names(from_name))
  check(sprintf("%s: .ARR attribute agrees with the array name", v),
        identical(from_name[[v]], targets[[v]]))
say("")


# ---- 4. align the matrix to the sheet --------------------------------------
# targets.csv is already ordered subject / period / time. Order the matrix TO
# the sheet rather than sorting both, so the canonical order lives in exactly
# one place -- the input file.

check("sheet and matrix describe the same 95 arrays",
      setequal(targets$array, colnames(expr)))

expr <- expr[, targets$array]
check("matrix columns now follow targets.csv", identical(colnames(expr), targets$array))

# From here on, columns are readable sample ids (S01_T1_PRE_B), not TAC strings.
colnames(expr) <- targets$sample_id
check("matrix columns are the sample ids", identical(colnames(expr), targets$sample_id))
say("")


# ---- 5. factors, with levels stated out loud -------------------------------
# R orders factor levels alphabetically. Left alone, time becomes
# 24hr / 3hr / Pre and every plot axis and contrast silently reverses.

targets$subject   <- factor(targets$subject)
targets$period    <- factor(targets$period,    levels = c("T1", "T2"))
targets$time      <- factor(targets$time,      levels = c("Pre", "3hr", "24hr"))

# B is the reference, so a positive log2 fold change always means "higher in A".
# Fix the direction once, here, and never think about it again.
targets$condition <- factor(targets$condition, levels = c("B", "A"))
targets$sequence  <- factor(targets$sequence,  levels = c("B", "A"))

grp_built <- paste(targets$condition, targets$time, sep = "_")
check("rebuilt grp agrees with the grp column in targets.csv",
      identical(grp_built, targets$grp))
targets$grp <- factor(grp_built,
                      levels = c("B_Pre", "B_3hr", "B_24hr", "A_Pre", "A_3hr", "A_24hr"))
say("")


# ---- 6. the design, stated as facts ----------------------------------------
say("--- design ---")
say("subjects: %d", nlevels(targets$subject))

n_per_subject <- table(table(targets$subject))
say("arrays per subject: %s",
    paste(sprintf("%s arrays x %d subjects", names(n_per_subject), n_per_subject),
          collapse = "; "))

cells <- table(targets$condition, targets$time)
say("condition x time cell counts:")
say("%s", paste(capture.output(print(cells)), collapse = "\n"))

# 95 of a planned 96. S19 has no A_Pre array, so that one cell is n = 15.
check("every cell is n = 16 except A_Pre at n = 15",
      all(cells == 16 | (rownames(cells)[row(cells)] == "A" &
                         colnames(cells)[col(cells)] == "Pre" & cells == 15)))

seq_tab <- table(targets$sequence[!duplicated(targets$subject)])
say("sequence allocation: %s", paste(names(seq_tab), seq_tab, sep = "-first = ", collapse = ", "))

# Sequence came out imbalanced, so condition and period are correlated in this
# sample. phi is the size of that association, and |phi| is the fraction of any
# period effect that leaks into a simple treatment effect. It is the reason
# `period` belongs in the model at Stage 02.
p1  <- table(targets$condition, targets$period)
phi <- (p1[1,1]*p1[2,2] - p1[1,2]*p1[2,1]) / sqrt(prod(rowSums(p1)) * prod(colSums(p1)))
say("condition-period association: phi = %.3f  (bias in a simple effect = %.2f x period effect)",
    phi, abs(phi))
say("")


# ---- 7. control probesets --------------------------------------------------
# A guard, not a filter. TAC's gene-level export already drops the AFFX
# controls, so this removes nothing here -- but it would catch a future export
# that included them, and it documents that they are absent rather than
# leaving a reader to assume it.

n_affx <- sum(grepl("^AFFX", rownames(expr)))
say("AFFX control probesets found: %d", n_affx)
if (n_affx > 0) expr <- expr[!grepl("^AFFX", rownames(expr)), ]

# No expression-level filter is applied at this stage. Filtering changes the
# multiple-testing burden, so it has to be justified against a specific
# contrast (independent filtering, Bourgon 2010) -- that belongs at Stage 02,
# not here, where it would silently set the denominator for everything.
say("carrying %d transcript clusters x %d arrays forward", nrow(expr), ncol(expr))
say("")


# ---- 8. save ----------------------------------------------------------------
check("final alignment holds", identical(colnames(expr), as.character(targets$sample_id)))

saveRDS(list(expr = expr, targets = targets), file.path(DATA, "01_expression.rds"))
say("wrote %s", file.path("01_Preprocessing", "02_data", "01_expression.rds"))

writeLines(.log,                          file.path(DATA, "01_alignment_checks.txt"))
writeLines(capture.output(sessionInfo()), file.path(DATA, "sessionInfo_01.txt"))

cat("\nStage 01 complete.\n")
