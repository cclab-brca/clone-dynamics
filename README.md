# clone-dynamics

This Github directory includes information relevant to the following manuscript from the Nguyen and Caldas labs:

Nguyen LV et al. Dynamics and plasticity of human breast cancer single cell-derived clones.

The CODE folder contains the scripts used for data analysis. Specifically, DTime_Sim and Plots_Sim contain the code for the mathematical modelling presented in Figure 2A and Supplementary Figure 5. 

PDX_cln_scRNAseq_analysis.dmg contains the code for all the scRNAseq data analysis presented in the paper.

The PROCESSED_DATA folder contains the processed data tables required for the analysis. This includes the input dataframe (hist.Rda) required to run DTime_Sim and Plots_Sim, and mcid_QC_v2.rda which contains the metacell ID (mcid), cell type, lentiviral-based cellular barcode sequence (var_seq), and doubling time information pertaining to each single cell RNA profile. sc_annotation.csv indicates the cell ID suffixes corresponding to each xenograft sample and clones_info.csv provides the lentiviral barcode sequence corresponding to each of the clones highlighted in the manuscript.

Processed count matrices from both bulk and scRNAseq can be accessed through Zenodo: https://zenodo.org/records/10978990. Raw sequencing data will be submitted to an appropriate repository and made available upon publication.

FOR REVIEWERS: Encrypted folders can be accessed using a password, which is the first word of the Introduction section of the manuscript. Instructions to reviewers on how to access the raw sequencing data are provided in the manuscript submission (either in the Code availability statement and/or the Reporting Summary document).
