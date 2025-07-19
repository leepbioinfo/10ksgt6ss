# Source README

## Overview

This directory contains all the **Python** scripts used for the analyses in this project.  
Each file generates one or more figures, loads a dataset or performs some analysis.

## Source Structure

- `genomes.py`: this script loads all genome annotations and merges the resulting table with domain annotations identified using from different Hidden Markov Model databases: ([Pfam](https://www.ebi.ac.uk/interpro/), [CDD](https://www.ncbi.nlm.nih.gov/cdd), [TXSScan](https://github.com/macsy-models) and models by [Zhang *et al.* (2012)](https://ftp.ncbi.nlm.nih.gov/pub/aravind/TOXIMM/)
- `filter_t6ss_neighbors.py`: Script to filter results from the T6SSiii_tssH model (which retrieves many non-T6SS proteins) and fetch 10 genes upstream and downstream of the identified T6SS components.
- `jaccard.py`: Script that uses the filtered data and classifies T6SS clusters using Jaccard and Louvain algorithms..
- `working_dfs.py`: Script that generates multiple dataframes incorporating manually curated annotations; these final dataframes serve as the basis for figure generation and statistical analyses.
- `Fig_x.py`: Scripts used to generate the figures for the final manuscript; some figures may have been subsequently edited with image software to enhance their aesthetics.
- `requirements.txt`: List of all Python dependencies.

## How to Reproduce the Analyses

1. **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

2. **Run the source script:**
    ```bash
    python source.py
    ```
    - Some analyses may require specific data files in the `data/` subfolder.

3. **Output:**
    - Results (plots, tables, etc.) will be generated in the same directory or in an `output/` subfolder, as described in the `source.py`.

## Notes

- All code is organized so that each analysis can be reproduced independently.
- **No code was modified after the results were finalized for publication.** These scripts represent exactly what was used to generate the figures and tables in the paper.
- Utility scripts in `shared/` may be imported by multiple analyses for common tasks (e.g., parsing, plotting).

## Questions or Issues

If you encounter any problems reproducing the results or have questions about the code, please open an [issue](https://github.com/leepbioinfo/10ksgt6ss/issues) or contact me at [nicastrogg@nih.gov].

---

*This repository is intended to maximize transparency and reproducibility in computational research. Feel free to reuse or adapt the scripts for your own work, citing this project where appropriate.*
