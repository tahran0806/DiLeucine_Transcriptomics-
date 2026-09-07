# DiLeucine_Transcriptomics

A project that runs transcriptomic statistical analyses on gene-level data from an acute resistance training study with dileucine supplementation versus placebo.
This stage-by-stage pipeline analyzes crossover microarray data in R. 

### About the project 

Sixteen participants completed an acute resistance training study that involved the consumption of either a placebo (PLA) or a di-leucine (DiL) supplement before a resistance training bout. Skeletal muscle biopsies were extracted from the vastus lateralis at three time points: PRE, 3hr POST, and 24hr POST.

Each participant completed both conditions, with each condition separated by a week. 

### Transcriptomic Tissue Preprocessing

RNA from VL tissue (~15 mg) was isolated using TRIzol. RNA suspended in DEPC-treated water was shipped on dry ice to a commercial laboratory for RNA integrity checks and transcriptome-wide analysis using the Clariom S Assay Human mRNA array. Raw .CEL files were analyzed using Transcriptome Analysis Console, and gene-level data were extracted from the software for downstream analysis in R.

The data went through a TAC pipeline before being exported from the device and then analyzed in R. 

Signal Space Transformation (SST) - Robust Multi-Chip Analysis (RMA).

What is SST-RMA?

It is the default preprocessing and normalization algorithm used by Thermo Fisher Scientific for modern whole-transcriptome microarrays, including Clariom D and Human Transcriptome Arrays (HTA). 

SST-RMA uses two main phases:
1. Signal Space Transformation: This applies a GC - content correlation to the raw data based on the sequences of the probes. It scales and transforms the raw fluorescence signal.
2. Robust Multi-Chip Analysis: This passes the above corrected data through normalization workflows that include:
  - Background correction to remove noise
  - Quantile normalization to make overall distribution of intensities identical across all sample arrays
  - Summarization to combine individual probe values into a single final expression value for each gene or transcript.

This improves inferential capacity: 
- The correction of fold change scores allows more accurate interpretation of true biological variations
- The final values after applying the SST-RMA are automatically log2-transformed
