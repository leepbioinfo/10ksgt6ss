# Source README

## Overview

This directory contains all the **Python** scripts used for the analyses in this project.  
Each file generates one or more figures, loads a dataset or performs some analysis.

## Required datsets

To be able to reproduce the results available in this repository, input datasets should be retrieved from three different sources:

1. This repository's [data](../data) directory.
2. The [Zenodo repository](https://zenodo.org/records/16358274) associated with this work, were very large files  (> 500Gb) are provided.
3. The genome assemblies from BioProjects PRJEB35182 and PRJEB47910, available at [NCBI](https://www.ncbi.nlm.nih.gov/) website.

## Source Structure

- Code that requires changing file paths and/or download of publicily available data
  - `genomes.py`: this script loads all genome annotations, MMseqs2 protein clusters and merges the resulting table with domain annotations identified using from different Hidden Markov Model databases: [Pfam](https://www.ebi.ac.uk/interpro/), [CDD](https://www.ncbi.nlm.nih.gov/cdd), [TXSScan](https://github.com/macsy-models) and models by [Zhang *et al.* (2012)].(https://ftp.ncbi.nlm.nih.gov/pub/aravind/TOXIMM/)
  - `filter_t6ss_neighbors.py`: Script to remove results from the T6SSiii_tssH model (which retrieves many non-T6SS proteins) and fetch 10 genes upstream and downstream of the identified T6SS components.
  - `jaccard.py`: Script that uses the filtered data and classifies T6SS clusters using Jaccard's distance and the Louvain community detection algorithm.
  - `working_dfs.py`: Script that loads multiple dataframes and merges this data with manually curated annotations. These dataframes encompasses all data needed for figure generation and statistical analyses.
- Code that can be executed as soon as [the environment is setup](../README.md)
  - `fig_<figure number/item>.py`: Scripts used to generate the figures for the final manuscript. In the manuscript, some figures may have been edited with image software to enhance their quality and add other informative elements.

## Notes

- All code for the figures is organized so that each script can be run independently.
- **No code was modified after the results were finalized for publication.** These scripts represent exactly what was used to generate the figures and tables in the paper.

## Questions or Issues

If you encounter any problems reproducing the results or have questions about the code, please open an [issue](https://github.com/leepbioinfo/10ksgt6ss/issues).

---

*This repository is intended to maximize transparency and reproducibility of computational research. Feel free to reuse or adapt the scripts for your own work, citing this project where appropriate.*
