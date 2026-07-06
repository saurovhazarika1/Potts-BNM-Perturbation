#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--bins", required=True, nargs=2, type=float)
    p.add_argument("--alpha", default=1.0, type=float)
    p.add_argument("--eps", default=1e-8, type=float)
    p.add_argument("--output-prefix", default="target")
    return p.parse_args()

def assign_bins(x, bins):
    b1, b2 = bins
    s = np.zeros(len(x), dtype=int)
    s[(x >= b1) & (x < b2)] = 1
    s[x >= b2] = 2
    return s

def marginal(states, q=3):
    counts = np.array([(states == a).sum() for a in range(q)])
    return pd.DataFrame({"state": np.arange(q), "count": counts, "probability": counts / counts.sum()})

def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    amber = df["amber_energy"].to_numpy(float)
    wt_state = df["wt_state"].to_numpy(int)
    frus_wt = df["frus_wt_energy"].to_numpy(float)
    frus_mut = df["frus_mut_energy"].to_numpy(float)

    emin = amber.min()
    amber_shifted = amber - emin + args.eps
    delta_frus = frus_mut - frus_wt
    r_frus = delta_frus / (np.abs(frus_wt) + args.eps)

    amber_mut_shifted = amber_shifted * (1.0 + args.alpha * r_frus)
    amber_mut_shifted = np.maximum(amber_mut_shifted, args.eps)
    amber_mut = amber_mut_shifted + emin - args.eps
    target_state = assign_bins(amber_mut, args.bins)

    wt = marginal(wt_state).rename(columns={"count":"wt_count","probability":"wt_probability"})
    tg = marginal(target_state).rename(columns={"count":"target_count","probability":"target_probability"})
    summary = wt.merge(tg, on="state")

    out = df.copy()
    out["delta_frus"] = delta_frus
    out["r_frus"] = r_frus
    out["amber_mut"] = amber_mut
    out["target_state"] = target_state
    out["state_changed"] = wt_state != target_state

    prefix = Path(args.output_prefix)
    out.to_csv(prefix.with_suffix(".per_frame.csv"), index=False)
    summary.to_csv(prefix.with_suffix(".target_marginal.csv"), index=False)
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
