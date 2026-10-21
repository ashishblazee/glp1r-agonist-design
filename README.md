# De Novo GLP1R Agonist Design

Computational design of de novo functional agonists targeting human Glucagon-Like Peptide-1 Receptor (**GLP1R**) for **Challenge 3 of the Anthropic × Adaptyv Protein Design Competition**.

Author: [Ashish](https://github.com/ashishblazee)

---

## Overview

Unlike conditional binder challenges, GLP1R requires **functional agonism** (switching on intracellular $G_s$ signaling) rather than passive extracellular binding. 

Analysis of active GLP1R cryo-EM coordinates (**PDB: 6X18**) reveals a two-domain activation mechanism:
* **The Extracellular Domain (ECD) Groove**: Anchors the C-terminal helix for high affinity. Pure ECD binders act as antagonists.
* **The 7TM Transmembrane Bundle**: The N-terminal trigger (His1–Asp9) inserts deep into the orthosteric pocket to switch on $G_s$ coupling.

This project implements a modular de novo framework: an unhindered N-terminal catalytic warhead connected via calibrated linkers (35–50 Å reach) to de novo structural scaffolds (3-helix bundles, $\alpha/\beta$ clamps, and VHH nanobodies).

---

## Key Features

* **Distance-Calibrated Linkers**: Built to bridge the measured 36.6 Å gap between the ECD pocket centroid and the 7TM orthosteric cavity.
* **De Novo Novelty Prescreen**: Engineered bio-isosteric substitutions preserving key catalytic contacts (His1, Glu3, Phe6, Asp9) while ensuring no continuous $\ge 6$-mer matches against known peptides.
* **Non-Redundant Portfolio**: 20 distinct designs with $<70\%$ pairwise sequence identity across 8 different structural topologies.
* **In Silico Orthogonal Consensus**: Evaluated with both **Meta ESMFold** and **Google DeepMind AlphaFold 3** (scaffold pLDDT up to 85.2, warhead pLDDT 71.5).

---

## Repository Structure

```text
├── scripts/
│   ├── measure_geometry.py          # 3D vector and distance calculations from PDB 6X18
│   ├── contact_analysis.py          # Atomic contact mapping (<4.0 Å)
│   ├── fold_candidate.py            # Automated Meta ESMFold prediction pipeline
│   ├── test_denovo_warheads.py      # Substring novelty check (<6 aa matching)
│   └── generate_diverse_batch_1.py  # Portfolio generator with <70% identity clustering
├── data/
│   ├── proteinbase_submission_batch_1_diverse.csv   # 20 portal-verified de novo designs
│   └── glp1r_candidate_library_ranked.csv          # 40-design ranked library with biophysical metrics
├── models/
│   ├── GLP1R_DES_A01_esmfold.pdb    # ESMFold predicted atomic coordinates
│   └── GLP1R_DES_A01_alphafold3.cif # AlphaFold 3 predicted atomic coordinates
└── README.md
```

---

## Quickstart

### Prerequisites
* Python 3.10+
* Standard libraries (`urllib`, `math`, `csv`, `os`)

### 1. Measure 6X18 Cryo-EM Target Geometry
```bash
python scripts/measure_geometry.py
```

### 2. Run In Silico Structure Prediction (ESMFold API)
```bash
python scripts/fold_candidate.py
```

### 3. Generate Non-Redundant Candidate Batches
```bash
python scripts/generate_diverse_batch_1.py
```

---

## References

1. **Active GLP1R Structure**: Zhang et al., *Nature* (2020) — PDB [6X18](https://www.rcsb.org/structure/6X18).
2. **Competition**: [Anthropic × Adaptyv Protein Design Competition](https://proteinbase.com).
3. **Structure Prediction**: AlphaFold 3 (Google DeepMind) & ESMFold (Meta AI).

---

## License
MIT License.
