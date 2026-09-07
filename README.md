# DiLeucine_Transcriptomics

A project that runs transcriptomic statistical analyses on gene-level data from an acute resistance training study with dileucine supplementation versus placebo.
This stage-by-stage pipeline analyzes crossover microarray data in R. 

##About the project:
Sixteen participants completed an acute resistance training study that involved the consumption of either a placebo (PLA) or a di-leucine (DiL) supplement before a resistance training bout. Skeletal muscle biopsies were extracted from the vastus lateralis at three time points: PRE, 3hr POST, and 24hr POST.

Each participant completed both conditions, with each condition separated by a week. 

###Transcriptomic Tissue Preprocessing
RNA from VL tissue (~15 mg) was isolated using TRIzol. RNA suspended in DEPC-treated water was shipped on dry ice to a commercial laboratory for RNA integrity checks and transcriptome-wide analysis using the Clariom S Assay Human mRNA array. Raw .CEL files were analyzed using Transcriptome Analysis Console, and gene-level data were extracted from the software for downstream analysis in R.
