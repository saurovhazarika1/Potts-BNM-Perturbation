# Maximum Entropy Derivation

## Optimization

```math
P_{mut}
=
\arg\min_P
D_{KL}(P||P_{WT})
```

where

```math
D_{KL}
=
\sum_\sigma
P(\sigma)
\ln
\frac{P(\sigma)}
{P_{WT}(\sigma)}.
```

Subject to

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

## Lagrangian

```math
\mathcal L
=
D_{KL}
+
\gamma
\left(
\sum_\sigma P(\sigma)-1
\right)
+
\sum_a
\lambda_m(a)
\left[
\sum_\sigma
P(\sigma)
\mathbf1[\sigma_m=a]
-
P_m^{target}(a)
\right].
```

Differentiate:

```math
\frac{\partial\mathcal L}
{\partial P(\sigma)}
=
\ln\frac{P(\sigma)}{P_{WT}(\sigma)}
+1+\gamma+
\sum_a
\lambda_m(a)
\mathbf1[\sigma_m=a]
=0.
```

Therefore,

```math
P_{mut}(\sigma)
=
\frac1Z
P_{WT}(\sigma)
\exp
\left[
\sum_a
\lambda_m(a)
\mathbf1[\sigma_m=a]
\right].
```

The multipliers are iteratively optimized until

```math
P_m^{mut}(a)
=
P_m^{target}(a).
```
