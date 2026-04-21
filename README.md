# Spatial-Transcriptomics-Masther-Thesis
This repository contains the computational analysis developed for my Master's thesis on spatial transcriptomics in bladder cancer. The project focuses on understanding the cellular composition and spatial organization of the tumor microenvironment using high-resolution spatial gene expression data. More concretely, this involves analyzing post-treatment samples to compare tumor evolution across patients and treatment types.

The workflow includes data preprocessing, quality control, clustering, and cell type annotation, integrating tools from the Seurat ecosystem and publicly available marker databases (e.g., PanglaoDB). Special emphasis is placed on identifying distinct cellular populations and their spatial relationships within tumor tissue, which may provide insights into tumor heterogeneity and immune infiltration.

In addition, the repository explores clustering strategies (e.g., Louvain vs. Leiden) and marker-based annotation approaches to improve the robustness of cell type identification in spatial transcriptomics datasets.

Overall, this work aims to demonstrate how spatially resolved transcriptomic data can reveal biologically meaningful patterns in bladder cancer tissues, contributing to a better understanding of the tumor microenvironment.

# Repository Structure
Because this project is still ongoing, some scripts remain unfinished. The documents that start with **“Playing”** are exploratory; they were mainly used to test ideas and approaches and are therefore likely incomplete. Once the most suitable method was identified, we proceeded to create the **“Final”** scripts. Files that start with **“Final”** are therefore the ones most likely to be complete.

```
Spatial-Transcriptomics-Master-Thesis/
│
├── Code/                                # Scripts and data used for the analysis
│   ├── Annotations/                     # Annotation scripts and related resources
│   │
│   ├── DUTRENEO/                        # Spatial data from the DUTRENEO clinical trial
│   │   ├── RAW/                         # Raw data (not uploaded for privacy/safety reasons)
│   │   └── QC/                          # Preprocessed data (not uploaded for privacy/safety reasons)
│   │
│   ├── 1_QC_Visium.R                    # Quality control for Visium data
│   ├── 1_QC_Xenium.R                    # Quality control for Xenium data
│   ├── Final_QC_Xenium.R                # Final QC pipeline for Xenium
│   │
│   ├── Playing_*.Rmd                    # Exploratory scripts used while testing approaches
│   │
│   ├── Quality_controls_comparison.txt
│   ├── How_to_access_cluster_JIC.txt
│   └── Links_tutorials.txt
│
├── Papers (documentation)/              # Reference papers and supporting documentation
│
└── Thesis/                              # Thesis manuscript and LaTeX files
    ├── figures/                         # Figures used in the thesis
    ├── rho-class/                       # LaTeX thesis class/template
    ├── main.tex                         # Main LaTeX document
    ├── main.pdf                         # Compiled thesis
    └── rho.bib                          # Bibliography file
```
