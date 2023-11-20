# A TCF4-dependent gene regulatory network confers resistance to immunotherapy in melanoma

This repository contains scripts used to analyse scRNAseq and spatial transcriptomic data.

## List of scripts / steps of the analysis

1. Preprocessing of all samples including DoubletFinder (GC_ALL.R)
2. Inference of copy number variations (CNV; HB_GC_ALL.R)
3. Analysis of the malignant compartment only (GC_ALL_MALIGNANT.R)
4. Analysis of the malignant compartment only in the treatment naive subset (MALIGNANT_BT.R)
4. Config file for running SCENIC 50x using Nextflow (harmony_scenic50x.vsn-pipelines)
5. Mouse tumour heterogeneity in vivo and in vitro (Mouse_heterogeneity_label_transfer.R)
6. Example of Visium (10x Genomics) data analysis on one sample (Visium_example.R)

## Location of the raw data

The data was deposited to the EGA portal with the study number: EGAS00001006488

### Location of the rds file of the malignant cells

https://doi.org/10.48804/GSAXBN
