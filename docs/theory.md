# Theory

## Objective

Predict mutation-induced allosteric rewiring **without running mutant MD** by combining:

- Wild-type molecular dynamics
- Bayesian-network-constrained Potts models
- Frustratometer mutation energetics
- Maximum Entropy inference

### Workflow

```text
WT MD
  ↓
Potts-BNM
  ↓
Frustratometer
  ↓
Target marginal
  ↓
Maximum Entropy
  ↓
Mutant ensemble
  ↓
Inverse Potts
  ↓
Mutant Hamiltonian
```

## Wild-Type Potts Model

Each residue occupies one of three Potts states

```math
\sigma_i \in \{0,1,2\}.
```

The Hamiltonian is

```math
H_{WT}(\sigma)=
-\sum_i h_i(\sigma_i)
-\sum_{(i,j)\in E_{BN}}
J_{ij}(\sigma_i,\sigma_j).
```

The WT equilibrium distribution is

```math
P_{WT}(\sigma)=
\frac1Z
e^{-\beta H_{WT}(\sigma)}.
```

## Mutation Constraint

Frustratometer-derived energetic perturbations are converted into a target marginal

```math
P_m^{target}(a).
```

## Maximum Entropy

Infer the least-biased mutant ensemble

```math
P_{mut}
=
\arg\min_P
D_{KL}(P||P_{WT})
```

subject to

```math
\sum_\sigma P(\sigma)=1
```

and

```math
\sum_\sigma
P(\sigma)
\mathbf1[\sigma_m=a]
=
P_m^{target}(a).
```

The resulting distribution is

```math
P_{mut}(\sigma)
=
\frac1Z
P_{WT}(\sigma)
\exp\left[
\sum_a
\lambda_m(a)
\mathbf1[\sigma_m=a]
\right].
```

Finally, refit a new Potts model to obtain

```math
\Delta J_{ij}
=
J_{ij}^{mut}
-
J_{ij}^{WT}.
```
