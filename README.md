# CRISPRi Perturb-seq analysis in human cortical neurons

This repository contains the analytical code and downstream processing workflows used to analyze single-cell **CRISPRi Perturb-seq** data targeting **Autism Spectrum Disorder (ASD) risk genes** in H9 human embryonic stem cell (hESC)-derived cortical neurons.

---

## Overview

The analytical pipeline spans initial fastq/matrix processing, quality control, energy distance analysis, cell type annotation, trajectory analysis, gene program discovery via consensus matrix factorization, and differential expression analysis.

---

## Analysis Workflow

### 1. Upstream Data Processing
Initial raw sequencing processing, guide assignment, and cell-by-gene matrix generation are performed using the standardized upstream pipeline:
* **Pipeline Repository:** [Hon-lab/Perturb-Seq-Processing-Pipeline](https://github.com/Hon-lab/Perturb-Seq-Processing-Pipeline)

---

### 2. Energy Distance Analysis
To evaluate the phenotypic impact of individual guide perturbations relative to non-targeting controls:
* **Pipeline Repository:** [Chikara-Takeuchi/energy_dist_pipeline](https://github.com/Chikara-Takeuchi/energy_dist_pipeline)

---

### 3. Clustering & Cell Annotation
Filtering, clustering, and cell type annotation are performed using:
* `Filtering_and_Annotation.R`

---

### 4. Cell Proportion Analysis
Quantification and statistical testing of cellular composition shifts across target perturbations:
* `CellProportions.R`

---

### 5. Pseudotime & Trajectory Analysis
Reconstruction of differentiation trajectories to map perturbation effects along developmental progression:
* `Pseudotime.R`

---

### 6. Gene Program Analysis (cNMF)
Consensus Non-negative Matrix Factorization (cNMF) was utilized to identify gene programs. 
* **Method Reference:** [Kotliar et al., eLife 2019;8:e43803](https://elifesciences.org/articles/43803)
* **Official cNMF Tutorial:** [dylkot/cNMF PBMC Tutorial](https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_pbmc_example_data.ipynb)

Because original cNMF workflows were adapted and extended for custom downstream analysis, stepwise custom scripts are provided in execution order:

1. `Prepare_cNMF.sh` – Prepares inputs and submits execution scripts (internally invokes `factorize_cnmf.sh`).
2. `Combine_cNMF.py` – Consolidates output factorization iterations across parallel runs.
3. `SelectK_cNMF.ipynb` – Evaluates factor stability and error metrics to select the optimal number of programs ($K$).
4. `CompareScore_and_GeneOntology_cNMF.py` – Calculates program gene scores and executes Functional Annotation / Gene Ontology enrichment.
5. `GP_downstream.ipynb` – Vizualization.
6. `GP_human_development.ipynb` – Maps identified gene programs against human developmental reference dataset.
7. `GP_mixture_modelling.ipynb` – Performs mixture modeling on program usage matrix.

---

### 7. Differential Gene Expression
Differentially expressed genes (DEGs) across perturbed targets was conducted using **pySpade**:
* **DEG Tool Repository:** [Hon-lab/pySpade](https://github.com/Hon-lab/pySpade)

*Note: The same `pySpade` pipeline was applied to quantify differentially expressed **transposable elements (TEs)** quantified via the **scTE** tool.*

#### Downstream DEG Analyses:
* `CoDysregulatedPairs.R` – Identification of co-dysregulated lncRNA and protein coding gene pairs across perturbed targets.
* `DEGs.R` – Post-processing, visualization, and filtering of differentially expressed genes and TEs.

---
