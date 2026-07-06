#!/usr/bin/env python3
from __future__ import annotations
import argparse, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samples", required=True, type=Path)
    p.add_argument("--edges", required=True, type=Path)
    p.add_argument("--h-wt", required=True, type=Path)
    p.add_argument("--J-wt", required=True, type=Path)
    p.add_argument("--l2", default=1e-3, type=float)
    p.add_argument("--maxiter", default=500, type=int)
    p.add_argument("--output-prefix", default="refit")
    return p.parse_args()

def unpack(theta,N,q,E):
    n = N*q
    return theta[:n].reshape(N,q), theta[n:].reshape(E,q,q)

def energy(configs,h,J,edges):
    Eout = np.zeros(configs.shape[0])
    for i in range(h.shape[0]):
        Eout -= h[i, configs[:,i]]
    for eid,(i,j) in enumerate(edges):
        Eout -= J[eid, configs[:,i], configs[:,j]]
    return Eout

def empirical(samples,edges,N,q):
    E = len(edges)
    eh = np.zeros((N,q))
    eJ = np.zeros((E,q,q))
    for i in range(N):
        for a in range(q):
            eh[i,a] = np.mean(samples[:,i] == a)
    for eid,(i,j) in enumerate(edges):
        for a in range(q):
            for b in range(q):
                eJ[eid,a,b] = np.mean((samples[:,i] == a) & (samples[:,j] == b))
    return eh,eJ

def obj_grad(theta,configs,edges,eh,eJ,N,q,E,l2):
    h,J = unpack(theta,N,q,E)
    En = energy(configs,h,J,edges)
    score = -En
    score -= score.max()
    w = np.exp(score)
    p = w/w.sum()

    mh = np.zeros_like(eh)
    mJ = np.zeros_like(eJ)
    for i in range(N):
        for a in range(q):
            mh[i,a] = p[configs[:,i] == a].sum()
    for eid,(i,j) in enumerate(edges):
        for a in range(q):
            for b in range(q):
                mJ[eid,a,b] = p[(configs[:,i] == a) & (configs[:,j] == b)].sum()

    # negative log likelihood up to constant
    logZ = np.log(np.exp(-En - (-En).max()).sum()) + (-En).max()
    emp_score = np.sum(h*eh) + np.sum(J*eJ)
    val = logZ - emp_score + 0.5*l2*np.sum(theta*theta)
    grad = np.concatenate([(mh-eh+l2*h).ravel(), (mJ-eJ+l2*J).ravel()])
    return val, grad

def main():
    args = parse_args()
    samples = np.load(args.samples)
    h_wt = np.load(args.h_wt)
    J_wt = np.load(args.J_wt)
    edges = pd.read_csv(args.edges)[["i","j"]].to_numpy(int)
    N,q = h_wt.shape
    E = len(edges)
    if N > 10:
        raise ValueError("This exact refit is only for toy systems.")

    eh,eJ = empirical(samples,edges,N,q)
    configs = np.array(list(itertools.product(range(q), repeat=N)), dtype=np.int16)
    theta0 = np.concatenate([h_wt.ravel(), J_wt.ravel()])
    hist = []
    def fun(theta):
        val,grad = obj_grad(theta,configs,edges,eh,eJ,N,q,E,args.l2)
        hist.append(val)
        return val,grad

    res = minimize(fun,theta0,jac=True,method="L-BFGS-B",options={"maxiter":args.maxiter})
    h_mut,J_mut = unpack(res.x,N,q,E)
    prefix = Path(args.output_prefix)
    np.save(prefix.name+"_h_mut.npy", h_mut)
    np.save(prefix.name+"_J_mut.npy", J_mut)

    rows=[]
    for eid,(i,j) in enumerate(edges):
        d = J_mut[eid] - J_wt[eid]
        rows.append({"edge_id":eid,"i":i,"j":j,"delta_J_frobenius":float(np.linalg.norm(d)),
                     "J_wt_norm":float(np.linalg.norm(J_wt[eid])),"J_mut_norm":float(np.linalg.norm(J_mut[eid]))})
    pd.DataFrame(rows).to_csv(prefix.with_suffix(".delta_J.csv"), index=False)
    pd.DataFrame({"iteration":np.arange(len(hist)), "objective":hist}).to_csv(prefix.with_suffix(".training_history.csv"), index=False)
    print("success:", res.success)
    print(pd.DataFrame(rows).to_string(index=False))

if __name__ == "__main__":
    main()
