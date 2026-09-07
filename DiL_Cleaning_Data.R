here::here()
library(here)

dir_in   <- here("R", "Analysis", "00_Input")
dir_out  <- here("R", "Analysis", "01_Preprocessing", "01_Output")
dir_fig  <- here("R", "Analysis", "01_Preprocessing", "02_Figures")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

# DiLeucine transcriptomics — Stage 01: ingest and build sample sheet
#
# In : 00_Input/Raw_data.txt   (21,448 transcript clusters x 95 arrays,
#                               log2 SST-RMA signal exported from TAC 4.0)
# Out: 01_Output/targets.csv   (derived sample sheet)
#      01_Output/01_clean.rds  (expression matrix + targets, aligned)

# Step 1: Read the raw data file and make changes
raw <- read.delim(file.path(dir_in, "Raw_data.txt"),
                  row.names = 1, check.names = FALSE) #check.names = FALSE preserves your column names exactly

stopifnot(all(sapply(raw, is.numeric))) #The stopifnot guards as.matrix's mistake of having one chracter column
expr <- as.matrix(raw)

# Derive the sample sheet from the column names

pat <- "_([A-H][0-9]{2})_(S[0-9]+)(T[12])(Pre|3hr|24hr)([AB])\\.sst-rma"
m <- regmatches(colnames(expr), regexec(pat, colnames(expr)))
stopifnot(all(lengths(m) == 6))

targets <- data.frame(
  array     = colnames(expr),
  well      = sapply(m, `[`, 2),
  subject   = factor(sapply(m, `[`, 3)),
  period    = factor(sapply(m, `[`, 4)),
  time      = factor(sapply(m, `[`, 5), levels = c("Pre", "3hr", "24hr")),
  condition = factor(sapply(m, `[`, 6), levels = c("B", "A")),
  stringsAsFactors = FALSE
)

f <- targets[targets$period == "T1" & !duplicated(targets$subject), ]
targets$sequence <- factor(setNames(as.character(f$condition),
                                    as.character(f$subject))[as.character(targets$subject)])

targets$grp <- factor(paste(targets$condition, targets$time, sep = "_"))

stopifnot(nrow(targets) == ncol(expr),
          identical(targets$array, colnames(expr)))

write.csv(targets, file.path(dir_out, "targets.csv"), row.names = FALSE)
saveRDS(list(expr = expr, targets = targets), file.path(dir_out, "01_clean.rds"))

writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo_01.txt"))

# Checkpoints
dim(expr)                                              # 21448   95
table(targets$condition, targets$time)                 # 16 everywhere EXCEPT A/Pre = 15
table(targets$sequence[!duplicated(targets$subject)])  # A = 10, B = 6
sum(expr)                                              # 12843602.47

git init
