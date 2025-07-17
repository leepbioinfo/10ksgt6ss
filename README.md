# Genome-Directed Study of *Salmonella* T6SS Effectors

[![DOI](https://img.shields.io/badge/DOI-10.1101/2024.09.27.615498-blue)](https://doi.org/10.1101/2024.09.27.615498)

This repository contains the data, scripts, and resources associated with our study on the diversity of Type VI Secretion System (T6SS) effectors in *Salmonella* and the discovery of a novel family of lipid-targeting antibacterial toxins.

## 📰 Citation

Nicastro, G.G., Sibinelli-Sousa, S., Hespanhol, J.T., et al. (2025)  
**Genome-directed study reveals the diversity of *Salmonella* T6SS effectors and identifies a novel family of lipid-targeting antibacterial toxins**. *bioRxiv*. https://doi.org/10.1101/2024.09.27.615498

---

## 🧬 Summary

- Analyzed **10,000 *Salmonella* genomes** to identify T6SS-encoded effectors.
- Identified **128 effector groups**, including **45 novel or highly divergent domains**.
- Biochemically validated **Tox-Act1**, a phospholipase effector secreted by T6SS.
- Tox-Act1 confers **competitive advantage** in the mouse gut and represents the **first lipid-targeting NlpC/P60 toxin** described in this context.

---

## 🔬 Repository Contents

- `hmm/` – HMM profiles used in this study
- `sources/` – Python scripts for T6SS detection, clustering, figure production and annotation
- `html/` – Sequence alignment files in HTML format, used in the interactive supplementary data.
- `alns/` –Raw alignments in FASTA format for each model generated in this study, along with corresponding HHpred results in both .hhr (text) and HTML formats.
- `pdbs/` –Mol* sessions with AlphaFold-predicted structures for each toxic domain model created in this study.
- `data/` – Data folder containing Python pickle objects, Excel tables, and YAML files used in both manual analyses and automated scripts from the `source` folder.

---

## 🌐 Online Resource

- **Effector & Structure Browser**:  
  https://leepbioinfo.github.io/10ksgt6ss  
  Includes novel STox families, alignments, and AlphaFold structures.

---

## 📦 Requirements

- Python 3.8+
- Dependencies:
  - `pandas 2.0+`, `Biopython`, `matplotlib`
  - External tools: `mmseqs2`, `jackhmmer`, `mafft`, `hh-suite`, `FoldSeek`, `AlphaFold`, `DALI`

Install via `environment.yml` or `requirements.txt`.

---

## 🚀 Reproducing the Analysis

1. **Prepare input data**: GFF and protein FASTA files from the 10KSG dataset.
2. **Detect T6SS clusters**:  
   `python scripts/identify_T6SS_clusters.py`
3. **Effector mining and classification**:  
   Use scripts and HMMs to identify candidates.
4. **Structure prediction and comparison**:  
   AlphaFold2, FoldSeek, and DALI-based remote homology.
5. **Manual curation**:  
   Cross-validate domain predictions with structural context and conserved neighborhoods.

---

## 🔍 Key Findings

- Identification of serovar-specific effector repertoires in *Salmonella*.
- STox_15 (Tox-Act1) is a **novel phospholipase effector** with NlpC/P60 fold.
- T6SS subtypes i1 and i3 associate with distinct effector classes (anti-eukaryotic vs. antibacterial).

---

## 📄 License

This repository is licensed under the **Creative Commons BY-NC-ND 4.0** license.

---

## ✉️ Contact

For questions or feedback, contact:

- Gianlucca Nicastro – `gianlucca.nicastro@nih.gov`
- Robson F. de Souza – `rfsouza@usp.br`
- Ethel Bayer-Santos – `ebayersantos@austin.utexas.edu`
