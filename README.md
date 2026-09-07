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
  - 2a. Background correction to remove noise
      - Subtracts noise such as optical background from unwanted fluorescence and nonspecific hybridization.
      
  - 2b. Quantile normalization 
      - Forces every array to have an identical distribution of intensities across all sample arrays 
      
  - 2c. Median Polish Summarization
      - Collapses many probes into one value per transcript, using a robust fit across probes and arrays, as individual probes misbehave
      - Bad probes are outvoted rather than allowed to drag the summary

This improves inferential capacity: 
- The correction of fold change scores allows more accurate interpretation of true biological variations
- The final values after applying the SST-RMA are automatically log2-transformed

## `DiL_PCA`
This script works on understanding the variability in the data by implementing principal component analysis plots, looking at PC1(%) and PC2(%). The data matrix is composed of 95 samples, each describing around 21,500 numbers. To visualize across such high-dimensional degrees, a principal component analysis (PCA) crunches the data such that it looks for a single direction across all the 21,500 dimensions along which samples are spread out. This is considered PC1. The second principal component, PC2(%), stands perpendicular to PC1(%), which is the greatest direction of the remaining spread in the data.

The PCA plot looking at how samples are clustered by time (PRE, 3Hr post, 24Hr post) showcased:
* PC1 (37.1%) separates 24 hr from everything else.
* PC2 (8%) separates 3 hr downward from PRE, which clusters tightly in the upper left.

### PC1: 
PC1 is the delayed remodelling and inflammatory program. ECM (TNC, TIMP1, THBS1, SERPINE1, PRG4), mechanical stress (ANKRD1, ACTC1), and — importantly — markers of infiltrating immune cells: S100A8/S100A9 are calprotectin, expressed by neutrophils and monocytes, and CCL2 and SPP1 recruit and activate them.

### PC2 
PC2 is the acute transcriptional response. FOS, JUNB, EGR1, MYC, ATF3, NR4A1/3. ABRA (STARS), XIRP1, and CYR61 are mechanotransduction genes.


Additionally, the PCA plot colored by condition (PLA versus DiL) showed no compartmentalization, completely intermixed. However, treatment effects in this kind of study live in coordinated changes across specific pathways, of modest size — invisible to PCA by construction

Two samples sat pretty low on the PC2 scale. These were:

`S03_T2_3HR_B`
`S07_T1_PRE_A`

This was noted for downstream analysis of outliers.




## Terms to understand
**Probe**: A synthetic tool created by the manufacturer and attached to the physical microarray chip. This is, in most cases, an oligonucleotide that is present in thousands on a standard microarray chip that acts as `bait` to capture and bind to transcripts floating over the gene chip. They are designed using the known genetic sequence of an organism. They are engineered to be perfectly complementary to a very specific genetic target.

**Probe Level**: A single gene transcript is usually too long to be captured by just one short probe. Therefore, manufacturers design a probe set (a group of 10 to 25 different probes) that all target different parts of the exact same transcript.

**Transcript**: Naturally occurring RNA transcripts that are biologically transcribed in our cells from DNA. More specifically, a transcript is an actual messenger RNA (mRNA) or long non-coding RNA molecule that a cell produces when a gene is turned on (expressed). 

**Transcript Level**:  "Individual signal intensities from different probes in a set are mathematically crunched down into a single, final 'transcript-level' expression value.
