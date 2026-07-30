#!/usr/bin/env python
# ============================================================================
# 04_g4sp_topology.py   (env: g4sp)
#
# Headless batch topology prediction with G4ShapePredictor (Liew et al. 2024).
# The upstream tool ships a GUI; we reproduce its exact inference path (verified
# from "g4sp application code/G4ShapePredictor.py": seqs_converter + predict_G4 +
# get_threshold) so we can score thousands of PQS at once.
#
#   encode: {'A':1,'T':2,'C':3,'G':4,'N':0}; symmetric-pad to 100; trim from end.
#   model.predict_proba(X) on the cloned .pkl classifier (RandomForest default).
#   Topology call (two modes, set by g4sp.precision in config):
#     - "default": plain argmax (every valid PQS gets a class; G4SP's own default).
#     - 0.70..0.95: per-class precision-optimized thresholds from
#       (precision)optimized_threshold.pkl. A PQS keeps argmax only if some class
#       prob exceeds its threshold, else it is called -1 "mixed" (judgment
#       suspended). Higher precision => stricter => more "mixed".
#   Labels: 0 parallel, 1 antiparallel, 2 hybrid, -1 mixed, -2 invalid.
#
# Inputs:   results/tables/peak_pqs.csv  (peak_id, has_pqs, pqs_seq)
# Output:   results/tables/g4sp_topology.csv
#             peak_id, has_pqs, pqs_seq, g4sp_class, p_parallel, p_antiparallel,
#             p_hybrid, g4sp_confidence
#
# If the G4SP repo / model can't be loaded, the script still writes the CSV with
# empty class calls and exits 0 — step 05 then falls back to the pqsfinder loop
# heuristic so the pipeline always completes.
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

CLASS_LABEL = {0: "parallel", 1: "antiparallel", 2: "hybrid",
               -1: "mixed", -2: "invalid"}
PROBA_COL = {0: "p_parallel", 1: "p_antiparallel", 2: "p_hybrid"}
PAD_LEN = 100
ENC = {"A": 1, "T": 2, "C": 3, "G": 4, "N": 0}


def load_config(path):
    if yaml is None:
        raise RuntimeError("pyyaml not available to read config")
    with open(path) as fh:
        return yaml.safe_load(fh)


def encode_sequence(seq, pad_length=PAD_LEN):
    """Replicate G4SP's seqs_converter: integer encode, symmetric pad, trim end."""
    seq = str(seq).upper()
    v = [ENC[c] for c in seq]                       # assumes validated ATCGN
    while len(v) < pad_length:
        v.insert(0, 0)
        v.append(0)
    return v[:pad_length]                            # trims surplus from the end


def valid_sequence(seq):
    return isinstance(seq, str) and len(seq) > 0 and all(c in "ATCGN" for c in seq.upper())


def resolve_precision(precision):
    """Return a float precision in [0,1], or None for plain-argmax 'default' mode."""
    if isinstance(precision, bool):
        return None
    if isinstance(precision, (int, float)):
        return float(precision)
    if isinstance(precision, str) and precision.strip().lower() not in ("", "default", "none"):
        try:
            return float(precision)
        except ValueError:
            return None
    return None


def per_class_thresholds(thr_file, model_key, precision, classes):
    """Replicate G4SP get_threshold(): for each class, pick the precision key
    closest to the requested precision and take its 'threshold'. Returns a vector
    aligned to the model's predict_proba column order (model.classes_)."""
    per_class = thr_file[model_key]                  # {0:{prec:{'threshold':..}}, 1:..., 2:...}
    thr_by_label = {}
    for top in range(3):
        d = per_class[top]
        key = min(d.keys(), key=lambda x: abs(x - precision))
        thr_by_label[top] = float(d[key]["threshold"])
    return np.array([thr_by_label[int(c)] for c in classes], dtype=float)


def write_empty(out_path, df, reason):
    sys.stderr.write("WARNING: G4SP unavailable (%s).\n"
                     "         Writing empty topology calls; step 05 will use the "
                     "pqsfinder heuristic.\n" % reason)
    out = df[["peak_id", "has_pqs", "pqs_seq"]].copy()
    out["g4sp_class"] = np.array([np.nan] * len(out), dtype=object)
    out["p_parallel"] = np.nan
    out["p_antiparallel"] = np.nan
    out["p_hybrid"] = np.nan
    out["g4sp_confidence"] = np.nan
    out.to_csv(out_path, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--pqs", default="results/tables/peak_pqs.csv")
    ap.add_argument("--out", default="results/tables/g4sp_topology.csv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    g4sp_cfg = cfg.get("g4sp", {})
    repo = cfg["paths"]["g4sp_repo"]
    model_dir = os.path.join(repo, "g4sp application code")
    model_name = g4sp_cfg.get("model", "RandomForest (default).pkl")
    # accept either key name for the threshold file (back-compat with older config)
    thr_name = g4sp_cfg.get("threshold_file", g4sp_cfg.get("threshold",
                            "(precision)optimized_threshold.pkl"))
    precision = resolve_precision(g4sp_cfg.get("precision", "default"))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df = pd.read_csv(args.pqs)

    model_path = os.path.join(model_dir, model_name)
    if not os.path.exists(model_path):
        write_empty(args.out, df, "model not found at %s" % model_path)
        return
    try:
        with open(model_path, "rb") as fh:
            model = pickle.load(fh)
    except Exception as exc:                       # version-mismatch on unpickle, etc.
        write_empty(args.out, df, "could not unpickle model: %s" % exc)
        return

    classes = list(getattr(model, "classes_", [0, 1, 2]))

    # Build per-class threshold vector for precision mode; None => plain argmax.
    threshold_vec = None
    if precision is not None:
        thr_path = os.path.join(model_dir, thr_name)
        model_key = os.path.splitext(os.path.basename(model_name))[0]   # 'RandomForest (default)'
        try:
            with open(thr_path, "rb") as fh:
                thr_file = pickle.load(fh)
            threshold_vec = per_class_thresholds(thr_file, model_key, precision, classes)
            sys.stderr.write("G4SP: precision=%.2f -> per-class thresholds (by class %s) = %s\n"
                             % (precision, list(map(int, classes)),
                                np.round(threshold_vec, 4).tolist()))
        except Exception as exc:
            sys.stderr.write("WARNING: could not apply precision thresholds (%s); "
                             "falling back to argmax.\n" % exc)
            threshold_vec = None
    else:
        sys.stderr.write("G4SP: precision=default -> plain argmax (no 'mixed' class).\n")

    # Predict only rows with a valid PQS; leave others as NaN.
    has = df["has_pqs"].fillna(False).astype(bool)
    seqs = df["pqs_seq"].where(has, other="")
    is_valid = seqs.map(valid_sequence)
    pred_idx = np.where(has.values & is_valid.values)[0]

    df["g4sp_class"] = np.array([np.nan] * len(df), dtype=object)   # object dtype: holds strings
    for col in PROBA_COL.values():
        df[col] = np.nan
    df["g4sp_confidence"] = np.nan

    if len(pred_idx) > 0:
        X = np.array([encode_sequence(seqs.iloc[i]) for i in pred_idx])
        proba = np.asarray(model.predict_proba(X))         # (n, n_cls), cols in classes_ order
        for k, i in enumerate(pred_idx):
            row = proba[k]
            if threshold_vec is None or np.any(row > threshold_vec):
                cls = int(classes[int(np.argmax(row))])
            else:
                cls = -1                                   # mixed: no class clears its threshold
            df.at[df.index[i], "g4sp_class"] = CLASS_LABEL.get(cls, "mixed")
            for c_pos, c_lbl in enumerate(classes):
                if int(c_lbl) in PROBA_COL:
                    df.at[df.index[i], PROBA_COL[int(c_lbl)]] = float(row[c_pos])
            df.at[df.index[i], "g4sp_confidence"] = float(np.max(row))

    out = df[["peak_id", "has_pqs", "pqs_seq", "g4sp_class",
              "p_parallel", "p_antiparallel", "p_hybrid", "g4sp_confidence"]]
    out.to_csv(args.out, index=False)

    n_called = int(out["g4sp_class"].notna().sum())
    n_mixed = int((out["g4sp_class"] == "mixed").sum())
    sys.stderr.write("G4SP done. Classified %d / %d peaks (%d mixed; model=%s).\n"
                     % (n_called, len(out), n_mixed, model_name))


if __name__ == "__main__":
    main()
