# Metabolic modeling links leaf anatomy to environment-specific benefits on the C3-C4 spectrum

A computational framework for incorporating anatomical parameters into constraint-based models of plant photosynthesis. This project explores how anatomical traits, specifically the Bundle Sheath to Mesophyll volume ratio ($V_{BS}/V_{M}$) and suberisation influence metabolic fluxes and photosyntetic optimality across C3, C3–C4 intermediate, and C4 photosynthetic types.

---

## 📘 Overview

This repository contains the code, data, and documentation for simulating and analyzing the effect of leaf anatomical constraints on metabolic flux distributions. The project builds upon **stoichiometric models of plant core metabolism**, integrating tissue-specific compartmentalization and flux coupling between mesophyll and bundle sheath cells to create a **two cell-type**.

---

## 🧠 Objectives

- **Model Construction**: Generate a two-cell (Mesophyll & Bundle Sheath) metabolic model derived from a generic plant core model.
- **Anatomical Integration**: Constrain metabolic fluxes using anatomical parameters like cell volume ratios and plasmodesmatal conductance.
- **Simulation**: Perform weighted parsimonious Flux Balance Analysis (pFBA) under varying environmental (Light, CO2, N) and anatomical conditions.
- **Thermodynamics**: Analyze the thermodynamic feasibility of non-canonical ("Unidentified") flux solutions using **Component Contribution**.
- **Phylogeny**: Map anatomical traits onto phylogenetic trees to understand the evolutionary context of C4 traits.

---

## 🧩 Project Structure

The analysis is divided into sequential Jupyter Notebooks:

### 1. Model Building & Simulation
* **`[1]_Model_Creation.ipynb`**: Construction of the C4 model from a generic C3 base, including biomass reaction curation and transporter definitions.
* **`[2]_Anatomical_Data_Analysis.ipynb`**: Statistical analysis of empirical $V_{BS}/V_{M}$ data across C3, C3-C4, and C4 species.
* **`[3]_Parameter_Scans.ipynb`**: Execution of large-scale parameter scans (Volume Ratio, Light, CO2) across multiple scenarios (C3, C4, Variable, Suberised/Unsuberised).

### 2. Flux Analysis & Visualization
* **`[4]_Fluxmaps.ipynb`**: Generation of vector-based (SVG) metabolic flux maps visualizing pathway usage across different photosynthetic modes.
* **`[5]_Photosyntetic_Choice.ipynb`**: Comparative analysis of ATP turnover, quantum yield, and prediction of optimal photosynthetic types based on anatomy.
* **`[6]_Suberised_vs_Unsuberised.ipynb`**: Investigation into the metabolic impact of the suberin layer in Bundle Sheath cells.
* **`[7]_Bundle_Sheath_Leakage.ipynb`**: Sensitivity analysis of CO2 and O2 leakage rates on photosynthetic efficiency.

### 3. Empirical & Theoretical Validation
* **`[8]_Phylogenetic_Analysis.ipynb`**: Phylogenetic mapping of anatomical constraints using NCBI taxonomy.
* **`[9]_Unidentified_Flux_Modes.ipynb`**: Thermodynamic analysis, including Minimal Driving Force,  ($\Delta G'$) of novel/unidentified decarboxylation pathways using `equilibrator-api`.

```
project_root/
│
├── Datasets/            # Raw and processed datasets (metabolites, fluxes, parameters, phylogenetics)
├── Models/              # SBML models 
├── Notebooks/           # Numbered analysis notebooks (see above)
└── Images/              # Generated plots and flux maps

```

---

## ⚙️ Installation

Clone this repository and install the required dependencies specified in **`requirements.txt`**.

## 👩‍🔬 Author

**Tiago M. Machado**  
PhD candidate in Theoretical Plant Biology and Data Science  
University of Cologne – CEPLAS Cluster of Excellence

Final code optimization, and visualization refinement were completed by **Alejandro León-Ramírez**

---

## 🧭 Acknowledgements

Supervised by **Prof. Nadine Töpfer**.  

Developed as part of the **CEPLAS** initiative on theoretical and computational plant systems biology (EXC-2048/1 — project ID390686111).

CRC TRR 341 (project ID 456082119).

---

## 🧾 Citation

If you use this work, please cite:


