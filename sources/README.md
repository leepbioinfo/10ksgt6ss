# README

## Overview

This repository contains all the source files (`source.py`) used for the analyses presented in [project/paper name].  
Each directory corresponds to a distinct analysis or dataset, and contains the exact scripts used to generate the results, ensuring full reproducibility.

## Directory Structure


- **analysis_X/source.py**: Main script used for a specific analysis or figure.
- **analysis_X/data/**: Input data and/or intermediate files.
- **shared/**: Utility functions or modules used by multiple analyses.
- **requirements.txt**: List of all Python dependencies.

## How to Reproduce the Analyses

1. **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

2. **Navigate to the directory of the analysis you want to reproduce:**
    ```bash
    cd analysis_1
    ```

3. **Run the source script:**
    ```bash
    python source.py
    ```
    - Some analyses may require specific data files in the `data/` subfolder.
    - Each `source.py` includes documentation at the top explaining the required inputs and expected outputs.

4. **Output:**
    - Results (plots, tables, etc.) will be generated in the same directory or in an `output/` subfolder, as described in each `README.md`.

## Notes

- All code is organized so that each analysis can be reproduced independently.
- **No code was modified after the results were finalized for publication.** These scripts represent exactly what was used to generate the figures and tables in the paper.
- Utility scripts in `shared/` may be imported by multiple analyses for common tasks (e.g., parsing, plotting).

## Questions or Issues

If you encounter any problems reproducing the results or have questions about the code, please open an [issue](https://github.com/yourusername/yourrepo/issues) or contact me at [your_email@example.com].

---

*This repository is intended to maximize transparency and reproducibility in computational research. Feel free to reuse or adapt the scripts for your own work, citing this project where appropriate.*
