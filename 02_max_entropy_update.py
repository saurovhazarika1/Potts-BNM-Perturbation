#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass
class Model:
    h: np.ndarray
    J: np.ndarray
    edges: np.ndarray
    neighbors: list

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h", required=True, type=Path)
    p.add_argument("--J", required=True, type=Path)
    p.add_argument("--edges", required=True, type=Path)
    p.add_argument("--target-marginal", required=True, type=Path)
    p.add_argument("--mutation-index", required=True, type=int)
    p.add_argument("--beta", default=1.0, type=float)
    p.add_argument("--outer-iters", default=25, type=int)
    p.add_argument("--eta", default=0.7, type=float)
    p.add_argument("--n-steps", default=30000, type=int)
    p.add_argument("--burn-in", default=5000, type=int)
    p.add_argument("--thin", default=10, type=int)
    p.add_argument("--seed", default=11, type=int)
    p.add_argument("--save-samples", action="store_true")
    p.add_argument("--output-prefix", default="maxent")
    return p.parse_args()

def load_model(args):
    h = np.load(args.h)
    J = np.load(args.J)
    edges = pd.read_csv(args.edges)[["i","j"]].to_numpy(int)
    N, q = h.shape
    neigh = [[] for _ in range(N)]
    for eid, (i,j) in enumerate(edges):
        neigh[i].append((j,eid,True))
        neigh[j].append((i,eid,False))
    return Model(h,J,edges,neigh)

def load_target(path, q):
    df = pd.read_csv(path)
    col = "target_probability" if "target_probability" in df.columns else "probability"
    target = np.zeros(q)
    for _, row in df.iterrows():
        target[int(row["state"])] = float(row[col])
    return target / target.sum()

def local_energy(model, sigma, i, state, m, lam, beta):
    e = -model.h[i, state]
    if i == m:
        e -= lam[state] / beta
    for j, eid, is_i_side in model.neighbors[i]:
        sj = sigma[j]
        e -= model.J[eid, state, sj] if is_i_side else model.J[eid, sj, state]
    return e

def run_mc(model, m, lam, beta, n_steps, burn_in, thin, rng, sigma0=None):
    N, q = model.h.shape
    sigma = rng.integers(0, q, size=N) if sigma0 is None else sigma0.copy()
    samples = []
    for step in range(burn_in + n_steps):
        i = int(rng.integers(0, N))
        old = int(sigma[i])
        new = int(rng.integers(0, q-1))
        if new >= old:
            new += 1
        dE = local_energy(model,sigma,i,new,m,lam,beta) - local_energy(model,sigma,i,old,m,lam,beta)
        if dE <= 0 or rng.random() < np.exp(-beta*dE):
            sigma[i] = new
        if step >= burn_in and (step-burn_in) % thin == 0:
            samples.append(sigma.copy())
    return np.asarray(samples, dtype=np.int16)

def marginals(samples, q):
    N = samples.shape[1]
    P = np.zeros((N,q))
    for a in range(q):
        P[:,a] = (samples == a).mean(axis=0)
    return P

def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    model = load_model(args)
    N, q = model.h.shape
    m = args.mutation_index
    target = load_target(args.target_marginal, q)

    lam = np.zeros(q)
    sigma0 = None
    hist = []
    for it in range(args.outer_iters):
        samples = run_mc(model,m,lam,args.beta,args.n_steps,args.burn_in,args.thin,rng,sigma0)
        sigma0 = samples[-1].astype(int)
        P = marginals(samples,q)
        err = target - P[m]
        lam += args.eta * err
        lam -= lam[0]
        row = {"iteration":it, "max_abs_error":float(np.max(np.abs(err)))}
        for a in range(q):
            row[f"target_{a}"] = target[a]
            row[f"sample_{a}"] = P[m,a]
            row[f"lambda_{a}"] = lam[a]
        hist.append(row)
        print(f"iter {it:03d}: max_error={row['max_abs_error']:.4f}, lambda={lam}")

    final = run_mc(model,m,lam,args.beta,args.n_steps,args.burn_in,args.thin,rng,sigma0)
    wt = run_mc(model,m,np.zeros(q),args.beta,args.n_steps,args.burn_in,args.thin,rng)
    Pmut = marginals(final,q)
    Pwt = marginals(wt,q)
    R = 0.5 * np.abs(Pmut-Pwt).sum(axis=1)

    prefix = Path(args.output_prefix)
    pd.DataFrame({"state":np.arange(q), "lambda":lam, "target_probability":target}).to_csv(prefix.with_suffix(".lambda.csv"), index=False)
    pd.DataFrame(hist).to_csv(prefix.with_suffix(".history.csv"), index=False)

    rows = []
    for i in range(N):
        row = {"residue_index":i, "response_R":R[i]}
        for a in range(q):
            row[f"P_wt_{a}"] = Pwt[i,a]
            row[f"P_mut_{a}"] = Pmut[i,a]
        rows.append(row)
    pd.DataFrame(rows).to_csv(prefix.with_suffix(".mutant_marginals.csv"), index=False)

    if args.save_samples:
        np.save(prefix.with_suffix(".samples.npy"), final)

    print("target:", target)
    print("final sampled:", Pmut[m])

if __name__ == "__main__":
    main()
