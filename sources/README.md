# Source README

## Overview

This directory contains all the `source.py` files used for the analyses in this project.  
Each source corresponds to a distinct analysis or dataset, and contains the exact scripts used to generate the results, ensuring full reproducibility.

## Source Structure


- **filter_t6ss_neighbors.py**: Script to filter results from the T6SSiii_tssH model (which retrieves many non-T6SS proteins) and fetch 10 genes upstream and downstream of the identified T6SS components.
- **jaccard.py**: Script that uses the filtered data and classifies T6SS clusters using Jaccard and Louvain algorithms..
- **working_dfs.py**: Script that generates multiple dataframes incorporating manually curated annotations; these final dataframes serve as the basis for figure generation and statistical analyses.
- **Fig_x.py**: Scripts used to generate the figures for the final manuscript; some figures may have been subsequently edited with image software to enhance their aesthetics.
- **requirements.txt**: List of all Python dependencies.

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

If you encounter any problems reproducing the results or have questions about the code, please open an [issue](https://github.com/yourusername/yourrepo/issues) or contact me at [your_email@example.com].

---

*This repository is intended to maximize transparency and reproducibility in computational research. Feel free to reuse or adapt the scripts for your own work, citing this project where appropriate.*
