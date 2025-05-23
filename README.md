# clone-dynamics

This Github directory includes information relevant to the following manuscript from the Nguyen and Caldas labs:

Nguyen LV et al. Dynamics and plasticity of human breast cancer single cell-derived clones. Cell Reports 2025. doi: https://doi.org/10.1016/j.celrep.2025.115699

The CODE folder contains the scripts used for data analysis. Specifically, DTime_Sim and Plots_Sim contain the code for the mathematical modelling presented in Figure 2A and Supplementary Figure 5. 

The PDX_cln_scRNAseq_analysis folder contains the code for all the scRNAseq data analysis presented in the paper.

The PROCESSED_DATA folder contains the processed data tables required for the analysis. This includes the input dataframe (hist.Rda) required to run DTime_Sim and Plots_Sim, and mcid_QC_v2.rda which contains the metacell ID (mcid), cell type, lentiviral-based cellular barcode sequence (var_seq), and doubling time information pertaining to each single cell RNA profile. sc_annotation.csv indicates the cell ID suffixes corresponding to each xenograft sample and clones_info.csv provides the lentiviral barcode sequence corresponding to each of the clones highlighted in the manuscript.

Processed count matrices from both bulk and scRNAseq can be accessed through Zenodo: https://zenodo.org/records/14996776
