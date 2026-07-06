#!/usr/bin/env python3
from __future__ import annotations
import numpy as np
import pandas as pd

def build_neighbors(N, edges):
    neigh = [[] for _ in range(N)]
    for eid, (i, j) in enumerate(edges):
        neigh[i].append((j, eid, True))
        neigh[j].append((i, eid, False))
    return neigh

def local_energy(h, J, neigh, sigma, i, state):
    e = -h[i, state]
    for j, eid, is_i_side in neigh[i]:
        sj = sigma[j]
        e -= J[eid, state, sj] if is_i_side else J[eid, sj, state]
    return e

def run_mc(h, J, edges, n_steps=80000, burn_in=10000, thin=20, seed=3):
    rng = np.random.default_rng(seed)
    N, q = h.shape
    neigh = build_neighbors(N, edges)
    sigma = rng.integers(0, q, size=N)
    samples = []
    for step in range(burn_in + n_steps):
        i = int(rng.integers(0, N))
        old = int(sigma[i])
        new = int(rng.integers(0, q - 1))
        if new >= old:
            new += 1
        dE = local_energy(h, J, neigh, sigma, i, new) - local_energy(h, J, neigh, sigma, i, old)
        if dE <= 0 or rng.random() < np.exp(-dE):
            sigma[i] = new
        if step >= burn_in and (step - burn_in) % thin == 0:
            samples.append(sigma.copy())
    return np.asarray(samples, dtype=np.int16)

def main():
    rng = np.random.default_rng(7)
    N, q = 5, 3
    h = rng.normal(0.0, 0.25, size=(N, q))
    edges = np.array([[0,1],[1,2],[2,3],[3,4],[0,3]], dtype=int)
    J = rng.normal(0.0, 0.15, size=(len(edges), q, q))
    J[0] += 0.8 * np.eye(q)
    J[4] += np.array([[0.4,0.1,-0.1],[0.0,0.3,0.1],[-0.2,0.1,0.7]])

    samples = run_mc(h, J, edges)
    m = 0
    wt_state = samples[:, m]

    centers = np.array([-5.0, -3.0, -1.0])
    amber = centers[wt_state] + rng.normal(0.0, 0.35, size=len(wt_state))

    frus_wt = -1.5 + 0.2 * wt_state + rng.normal(0.0, 0.05, size=len(wt_state))
    delta_frus = 0.15 + 0.25 * wt_state + rng.normal(0.0, 0.05, size=len(wt_state))
    frus_mut = frus_wt + delta_frus

    df = pd.DataFrame({
        "frame": np.arange(len(wt_state)),
        "amber_energy": amber,
        "wt_state": wt_state,
        "frus_wt_energy": frus_wt,
        "frus_mut_energy": frus_mut
    })

    np.save("toy_h_wt.npy", h)
    np.save("toy_J_wt.npy", J)
    np.save("toy_wt_samples.npy", samples)
    pd.DataFrame(edges, columns=["i","j"]).to_csv("toy_edges.csv", index=False)
    df.to_csv("toy_frame_data.csv", index=False)
    print("Wrote toy data. Suggested bins: --bins -4.0 -2.0")

if __name__ == "__main__":
    main()
