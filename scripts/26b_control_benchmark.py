#!/usr/bin/env python
# ============================================================================
# 10b_control_benchmark.py   (env: g4sp)   --- topology calibration (Q3) ---
#
# Benchmark G4SP on a curated set of literature G4s with reported topology
# (data/known_topology_controls.tsv) using the SAME inference path as rule 04.
# Tells us whether G4SP is calibrated on canonical references (c-MYC/c-KIT ->
# parallel; TBA/telomere -> antiparallel/hybrid) before we trust the genome-wide
# parallel-dominance.
#
# Inputs:   calibration.controls_tsv; g4sp.model/threshold_file/precision
# Output:   results/tables/g4sp_control_benchmark.csv
#             name, reported_topology, predicted_topology, p_parallel,
#             p_antiparallel, p_hybrid, confidence, correct, note
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

PAD_LEN = 100
ENC = {"A": 1, "T": 2, "C": 3, "G": 4, "N": 0}
CLASS_LABEL = {0: "parallel", 1: "antiparallel", 2: "hybrid", -1: "mixed", -2: "invalid"}
PROBA_COL = {0: "p_parallel", 1: "p_antiparallel", 2: "p_hybrid"}


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def encode_sequence(seq, pad_length=PAD_LEN):
    seq = str(seq).upper()
    v = [ENC[c] for c in seq]
    while len(v) < pad_length:
        v.insert(0, 0)
        v.append(0)
    return v[:pad_length]


def valid_sequence(seq):
    return isinstance(seq, str) and len(seq) > 0 and all(c in "ATCGN" for c in seq.upper())


def resolve_precision(precision):
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
    per_class = thr_file[model_key]
    thr_by_label = {}
    for top in range(3):
        d = per_class[top]
        key = min(d.keys(), key=lambda x: abs(x - precision))
        thr_by_label[top] = float(d[key]["threshold"])
    return np.array([thr_by_label[int(c)] for c in classes], dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--out", default="results/tables/g4sp_control_benchmark.csv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    g4sp_cfg = cfg.get("g4sp", {})
    repo = cfg["paths"]["g4sp_repo"]
    model_dir = os.path.join(repo, "g4sp application code")
    model_name = g4sp_cfg.get("model", "RandomForest (default).pkl")
    thr_name = g4sp_cfg.get("threshold_file", "(precision)optimized_threshold.pkl")
    precision = resolve_precision(g4sp_cfg.get("precision", "default"))
    controls_tsv = cfg.get("calibration", {}).get(
        "controls_tsv", "data/known_topology_controls.tsv")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    ctrls = pd.read_csv(controls_tsv, sep="\t")

    model_path = os.path.join(model_dir, model_name)
    try:
        with open(model_path, "rb") as fh:
            model = pickle.load(fh)
    except Exception as exc:
        sys.stderr.write("ERROR: could not load G4SP model (%s). Writing empty benchmark.\n" % exc)
        ctrls.assign(predicted_topology=np.nan).to_csv(args.out, index=False)
        return

    classes = list(getattr(model, "classes_", [0, 1, 2]))
    threshold_vec = None
    if precision is not None:
        try:
            with open(os.path.join(model_dir, thr_name), "rb") as fh:
                thr_file = pickle.load(fh)
            model_key = os.path.splitext(os.path.basename(model_name))[0]
            threshold_vec = per_class_thresholds(thr_file, model_key, precision, classes)
        except Exception as exc:
            sys.stderr.write("WARNING: thresholds unavailable (%s); using argmax.\n" % exc)

    recs = []
    for _, r in ctrls.iterrows():
        seq = str(r["sequence"]).upper()
        rec = {"name": r["name"], "reported_topology": r["topology"],
               "note": r.get("note", "")}
        if not valid_sequence(seq):
            rec.update({"predicted_topology": "invalid", "p_parallel": np.nan,
                        "p_antiparallel": np.nan, "p_hybrid": np.nan,
                        "confidence": np.nan, "correct": False})
            recs.append(rec)
            continue
        proba = np.asarray(model.predict_proba([encode_sequence(seq)]))[0]
        if threshold_vec is None or np.any(proba > threshold_vec):
            cls = int(classes[int(np.argmax(proba))])
        else:
            cls = -1
        pred = CLASS_LABEL.get(cls, "mixed")
        for c_pos, c_lbl in enumerate(classes):
            if int(c_lbl) in PROBA_COL:
                rec[PROBA_COL[int(c_lbl)]] = float(proba[c_pos])
        rec.update({"predicted_topology": pred, "confidence": float(np.max(proba)),
                    "correct": bool(pred == str(r["topology"]))})
        recs.append(rec)

    out = pd.DataFrame(recs)[["name", "reported_topology", "predicted_topology",
                              "p_parallel", "p_antiparallel", "p_hybrid",
                              "confidence", "correct", "note"]]
    out.to_csv(args.out, index=False)
    acc = float(out["correct"].mean()) if len(out) else float("nan")
    sys.stderr.write("Control benchmark: %d sequences, accuracy = %.2f (precision=%s).\n"
                     % (len(out), acc, precision))


if __name__ == "__main__":
    main()
