# Genome-Directed Study of *Salmonella* T6SS Effectors

[![DOI](https://img.shields.io/badge/DOI-10.1101/2024.09.27.615498-blue)](https://doi.org/10.1101/2024.09.27.615498)

This repository contains the data, scripts, and resources associated with our study on the diversity of Type VI Secretion System (T6SS) effectors in *Salmonella* and the discovery of a novel family of lipid-targeting antibacterial toxins.

## 🧬 Summary

- This work analyzed **10,000 *Salmonella* genomes** to identify T6SS-encoded effectors.
- A total of **128 effector families** were identified, including **47 novel** (shared hits <= 25%) or **highly divergent*** domains (shared hits <= 75%) (Figure 2B).
- Serovar-specific effector repertoires are present in *Salmonella*.
- T6SS subtypes i1 and i3 associate with distinct effector classes (anti-eukaryotic vs. antibacterial).
- STox_15 (Tox-Act1) is a **novel effector** of the NlpC/P60 fold.
- Biochemical assays demonstrated **Tox-Act1** has phospholipase activity
- *In vivo* assays with mutant strains showed that Tox-Act1 is secreted by the T6SS.
- Tox-Act1 confers **competitive advantage** in the mouse gut and represents the **first lipid-targeting NlpC/P60 toxin** described in this context.

--

## 🌐 Online Resource

- **Effector & Structure Browser**:  
  https://leepbioinfo.github.io/10ksgt6ss  
  Includes novel STox families, alignments, and AlphaFold structures.

---

## 📰 Citation

Nicastro, G.G., Sibinelli-Sousa, S., Hespanhol, J.T., et al. (2025)  
**Genome-directed study reveals the diversity of *Salmonella* T6SS effectors and identifies a novel family of lipid-targeting antibacterial toxins**. *bioRxiv*. https://doi.org/10.1101/2024.09.27.615498

---

### Note on large files

This repository contains files larger than 100MB. Downloading these files requires Git LFS support.

```bash
conda install git-lfs
git lfs install
git lfs pull
```

### Setup the environment

```bash
conda env create -n <env-name> -f environment.yml
conda activate <env-name>
```

These two commands will install and enable the software environment needed to run scripts in this repository.

---

## 📦 Requirements

See `environment.yml`.

## 🔬 Repository Contents

- `hmm/` – HMM profiles used in this study
- `sources/` – Python scripts for T6SS detection, clustering, figure production and annotation
- `html/` – Sequence alignment files in HTML format, used in the interactive supplementary data.
- `alns/` –Raw alignments in FASTA format for each model generated in this study, along with corresponding HHpred results in both .hhr (text) and HTML formats.
- `pdbs/` –Mol* sessions with AlphaFold-predicted structures for each toxic domain model created in this study.
- `data/` – Data folder containing Python pickle objects, Excel tables, and YAML files used in both manual analyses and automated scripts from the `source` folder.

---

## 📄 License

This repository is licensed under the **Creative Commons BY-NC-ND 4.0** license.

---

## ✉️ Contact

For questions or feedback, contact:

- Gianlucca Nicastro – `gianlucca.nicastro@nih.gov`
- Robson F. de Souza – `rfsouza@usp.br`
- Ethel Bayer-Santos – `ebayersantos@austin.utexas.edu`
