# Data README

## Overview

This directory contains all the data files used in the analyses presented in this project.
Each file corresponds to a specific dataset or intermediate result, and serves as input or output for the various analysis steps, supporting full reproducibility and traceability.

## Data description


- **10k_vizinho_novo_df_jaccard.pkl**:  Pickle file containing a DataFrame of gene loci with at least one T6SS component gene; output from jaccard.py (in the source folder)..
- **FDEvolvedCargo5.xlsx**: Excel table listing protein IDs detected by HMM toxin models used in this study, including system classification and additional metadata.
- **Toxins_nomenclature 05-01-2024.xlsx**: Excel table tracking the renaming history of the models used in this study. External models include reference to the original publication..
- **figure_1E.tsv**: table of statistics comparing hits to HMM models from this resource to publicily available models. This data was used to generate figure 1E. For details on how this table was generated see 
- **from_rob.xlsx**: Excel table with annotations of toxic domains and T6SS component domains.
- **function2color.yaml**: YAML file specifying colors and ranking for thematic functions (e.g., peptidoglycan, nuclease, pore-forming), used to generate figures via scripts in the source folder.
- **jaccard_type_classification.yaml**:YAML file with manual mapping of gene neighborhood annotations based on Jaccard + Louvain clustering. Louvain communities were labeled (e.g., T6SSiii) using phylogenetic analysis of grouped proteins.
- **jaccard_group_classification.yaml**: YAML file mapping specific Jaccard clusters to broader functional groups.
- **md5.tsv**: MD5 checksums for all [alignments](../alns) and [HMM models](../hmm) available in this resource.
- **meta_s3.xlsx**: Metadata table for genomes used in the study; sourced from [Perez-Sepulveda et al., 2021](https://doi.org/10.1186/s13059-021-02536-3).
- **model2function.yaml**: YAML file mapping each HMM toxin model to a thematic function (e.g., peptidoglycan, nuclease, pore-forming) used in this study..
- **rename_model_name.yaml**:  YAML file used to rename selected models during the final analysis stage.
- **toxdist4.xlsx**: Excel table with references, statistics, and supporting information used to manually infer the function of toxin models.

