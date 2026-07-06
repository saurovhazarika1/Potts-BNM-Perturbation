# Potts-BNM-Perturbation

> A statistical-mechanical framework for predicting mutation-induced allosteric rewiring from wild-type molecular dynamics simulations.

---

## Overview

**Potts-BNM-Perturbation** is a theoretical framework for inferring how local perturbations (e.g., point mutations) reshape global allosteric communication networks **without requiring mutant molecular dynamics (MD) simulations**.

The framework combines

- Molecular Dynamics (MD)
- Bayesian Networks (BN)
- Potts models
- Frustratometer mutation energetics
- Maximum Entropy statistical inference

to reconstruct mutation-induced equilibrium ensembles and infer perturbed allosteric communication networks.

---

## Motivation

Large-scale mutational analysis using molecular dynamics is computationally expensive. While methods such as Frustratometer estimate the energetic consequences of mutations, they do not directly predict how mutations reorganize global allosteric communication.

This project aims to bridge that gap by integrating local energetic perturbations with a statistical-mechanical model learned from wild-type simulations.

---

## Conceptual Workflow

```text
Wild-Type MD
      │
      ▼
Residue Energy Decomposition
      │
      ▼
Potts State Assignment
      │
      ▼
Bayesian Network Learning
      │
      ▼
Inverse Potts Model
      │
      ▼
Wild-Type Hamiltonian
      │
Mutation
      │
      ▼
Frustratometer
      │
      ▼
Mutation-Induced Local Energy Perturbation
      │
      ▼
Target Potts State Distribution
      │
      ▼
Maximum Entropy Inference
      │
      ▼
Mutant Ensemble
      │
      ▼
Inverse Potts Refit
      │
      ▼
Mutant Hamiltonian
      │
      ▼
Allosteric Network Rewiring
```

---

# Theoretical Foundation

The wild-type protein is represented by a Potts Hamiltonian

```math
H(\sigma)
=
-
\sum_i h_i(\sigma_i)
-
\sum_{(i,j)\in E_{BN}}
J_{ij}(\sigma_i,\sigma_j)
```

where

- \(h_i\) are local energetic fields,
- \(J_{ij}\) are effective cooperative couplings,
- \(E_{BN}\) is the Bayesian-network-derived interaction graph.

The corresponding equilibrium ensemble is

```math
P_{WT}(\sigma)
=
\frac{1}{Z}
e^{-\beta H(\sigma)}.
```

Mutation-derived energetic perturbations are converted into local constraints, and the mutant ensemble is inferred using the principle of maximum entropy.

---

# Current Methodology

The current framework consists of the following steps:

1. Learn a Bayesian-network-constrained Potts model from wild-type MD simulations.
2. Compute mutation-induced energetic perturbations using Frustratometer.
3. Map energetic perturbations onto the MD-derived Potts representation.
4. Infer the least-biased mutant ensemble using maximum entropy.
5. Refit the Potts Hamiltonian using the inferred mutant ensemble.
6. Quantify mutation-induced rewiring of the allosteric communication network.

---

# Current Status

🚧 **Active research project**

The theoretical framework is under active development.

Current work focuses on

- rigorous mapping between Frustratometer energetics and Potts-state occupancies,
- maximum entropy inference,
- mutation-induced Hamiltonian reconstruction,
- validation against experimentally characterized proteins.

---

# Validation Systems

The framework will initially be validated using small proteins with extensive experimental mutation data.

- **GB1 (Protein G B1 domain)**
- **PDZ2 domain**
- **WW domain**

before being applied to larger biologically relevant systems, including nuclear receptors.

---

# Planned Applications

- Mutation-induced allosteric rewiring
- Communication pathway analysis
- Identification of allosteric control residues
- Disease-associated mutation analysis
- Protein engineering
- Rational mutagenesis
- Drug discovery

---

# Repository Structure

```text
Potts-BNM-Perturbation/
│
├── theory/              # Mathematical derivations
├── notebooks/           # Exploratory analyses
├── src/                 # Source code
├── examples/            # Example workflows
├── validation/          # Benchmark proteins
├── figures/             # Figures
├── docs/                # Documentation
└── README.md
```

---

# Citation

This repository accompanies the ongoing development of a statistical-mechanical framework for mutation-induced allosteric rewiring. A preprint and peer-reviewed publication will be linked here upon release.

---

# License

This project is released under the MIT License.

---

# Author

**Saurov Hazarika**

Postdoctoral Researcher  
SciLifeLab / Stockholm University (incoming)  
City of Hope (former)

---

## Disclaimer

This repository contains an actively developing theoretical framework. Mathematical formulations, algorithms, and implementation details may evolve as the methodology is refined and experimentally validated.
