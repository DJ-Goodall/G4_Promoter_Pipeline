#!/usr/bin/env python
# ============================================================================
# 10a_threshold_sweep.py   (env: g4sp)   --- topology calibration (Q3) ---
#
# Is the ~95%-parallel result a G4SP threshold artifact or robust? Sweep the
# precision setting and record how the parallel/antiparallel/hybrid/mixed
# fractions move. Works purely from the per-class probabilities already saved in
# g4sp_topology.csv (no model reload), reproducing G4SP get_threshold + the
# predict_G4 thresholded-argmax rule.
#
# Inputs:   results/tables/g4sp_topology.csv  (p_parallel/p_antiparallel/p_hybrid)
#           config: g4sp.model, g4sp.threshold_file; calibration.precision_grid
# Output:   results/tables/topology_precision_sweep.csv
#             precision, n, n_parallel, n_antiparallel, n_hybrid, n_mixed,
#             frac_parallel, frac_antiparallel, frac_hybrid, frac_mixed
# ============================================================================
import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None

# proba columns are class 0,1,2 in this order
PROBA_COLS = ["p_parallel", "p_antiparallel", "p_hybrid"]
LABELS = {0: "parallel", 1: "antiparallel", 2: "hybrid", -1: "mixed"}


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def per_class_thresholds(thr_file, model_key, precision):
    """G4SP get_threshold(): per class, the precision key closest to `precision`."""
    per_class = thr_file[model_key]
    out = []
    for top in range(3):
        d = per_class[top]
        key = min(d.keys(), key=lambda x: abs(x - precision))
        out.append(float(d[key]["threshold"]))
    return np.array(out, dtype=float)


def classify(probs, thresholds):
    """thresholds None -> plain argmax; else argmax if any prob>thr, else -1 mixed."""
    if thresholds is None or np.any(probs > thresholds):
        return int(np.argmax(probs))
    return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--g4sp", default="results/tables/g4sp_topology.csv")
    ap.add_argument("--out", default="results/tables/topology_precision_sweep.csv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    g4sp_cfg = cfg.get("g4sp", {})
    repo = cfg["paths"]["g4sp_repo"]
    model_dir = os.path.join(repo, "g4sp application code")
    model_key = os.path.splitext(os.path.basename(
        g4sp_cfg.get("model", "RandomForest (default).pkl")))[0]
    thr_name = g4sp_cfg.get("threshold_file", "(precision)optimized_threshold.pkl")
    grid = cfg.get("calibration", {}).get(
        "precision_grid", ["default", 0.70, 0.75, 0.80, 0.85, 0.90, 0.95])

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df = pd.read_csv(args.g4sp)
    P = df[PROBA_COLS].to_numpy(dtype=float)
    valid = ~np.isnan(P).any(axis=1)
    P = P[valid]
    n = len(P)

    thr_file = None
    thr_path = os.path.join(model_dir, thr_name)
    try:
        with open(thr_path, "rb") as fh:
            thr_file = pickle.load(fh)
    except Exception as exc:
        sys.stderr.write("WARNING: threshold file unloadable (%s); only 'default' "
                         "argmax row will be produced.\n" % exc)

    rows = []
    for prec in grid:
        is_default = (isinstance(prec, str) and prec.strip().lower() in ("default", "none", ""))
        if is_default:
            thresholds = None
            prec_label = "default"
        else:
            if thr_file is None:
                continue
            prec_f = float(prec)
            thresholds = per_class_thresholds(thr_file, model_key, prec_f)
            prec_label = "%.2f" % prec_f
        calls = np.array([classify(P[i], thresholds) for i in range(n)])
        counts = {k: int(np.sum(calls == v)) for k, v in
                  {"parallel": 0, "antiparallel": 1, "hybrid": 2, "mixed": -1}.items()}
        row = {"precision": prec_label, "n": n}
        for k in ["parallel", "antiparallel", "hybrid", "mixed"]:
            row["n_" + k] = counts[k]
            row["frac_" + k] = counts[k] / n if n else np.nan
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    sys.stderr.write("Threshold sweep done over %d precisions on %d peaks.\n"
                     % (len(rows), n))


if __name__ == "__main__":
    main()
