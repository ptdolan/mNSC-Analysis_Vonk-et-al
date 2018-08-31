# mNSC-Analysis_Vonk-et-al
RNAseq analysis of RNAseq data from differentiating murine NSCs

########
Figures: 
########

RPB_ChangePlot: Relative transcription at *gene* level. 

TPM_ChangePlot: Relative transcription at *transcript* level. PHN genes are highlighted. 

Diff_meanRPB_Scatter: Scatter of reads per base at the *gene* level.

Diff_meanTPM_Scatter: scatter of mean expression of *transcripts*. 

Smoothed Counts: the change in read counts between undifferentiated ("") and differentiated mNSCs ("D") on the x axis. 

all: All genes annotated as chaperone-, proteasome- or ubiquitin-related.

Chaps: Volcano plot of all Chaperone related genes. 'b' is a adjusted log2 fold-change value, it controls for experimental covariates (like age).

Chaps_grid: All Chaperone annotated genes, diveider by subset (Hsp90, Hsp70, sHSP, TRIC). 

PSM: Proteasome-related gene volcano. 

Ubiq: Ubiquitin related gene volcano. 

########
Sleuth App
########

Directory contains everything to view interactive results from RNAseq using sleuth and shiny in R. 

########
DataTables
########

Directory containing all results tables from fitted ststs models.
