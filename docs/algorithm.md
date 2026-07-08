# Algorithm

## Inputs

- WT Potts model
- Bayesian-network edges
- Frustratometer mutation energies
- WT AMBER residue energies

## Pipeline

1. Build the WT Potts model.
2. Compute Frustratometer mutation energetics.
3. Map mutation energetics onto WT AMBER residue energies.
4. Re-bin the perturbed energies.
5. Compute the target marginal

```math
P_m^{target}(a).
```

6. Solve the maximum-entropy problem

```math
P_{mut}
=
\arg\min
D_{KL}(P||P_{WT}).
```

7. Sample the mutant ensemble using Metropolis Monte Carlo.
8. Estimate

```math
P_i^{mut}(a),
\qquad
P_{ij}^{mut}(a,b).
```

9. Refit the Potts model.
10. Compute

```math
\Delta J_{ij}
=
J_{ij}^{mut}
-
J_{ij}^{WT}.
```
