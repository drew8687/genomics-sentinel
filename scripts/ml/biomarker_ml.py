#!/usr/bin/env python3
"""
GenomicsSentinel — ML Biomarker Discovery
Pipeline: Feature matrix from DE genes → RF + XGB + SVM →
          Optuna HPO → SHAP explanation → MLflow tracking
"""

import json
import logging
import warnings
from pathlib import Path

import mlflow
import mlflow.sklearn
import numpy as np
import optuna
import pandas as pd
import shap
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    roc_auc_score,
    roc_curve,
    matthews_corrcoef,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")
optuna.logging.set_verbosity(optuna.logging.WARNING)
logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
log = logging.getLogger("genomics-ml")

# ── Snakemake I/O ─────────────────────────────────────────
feature_matrix = snakemake.input["features"]
labels_file    = snakemake.input["labels"]
out_metrics    = snakemake.output["metrics"]
out_shap       = snakemake.output["shap_plot"]
out_roc        = snakemake.output["roc_plot"]

cfg    = snakemake.config
ml_cfg = cfg["ml"]

# ── Load data ─────────────────────────────────────────────
log.info("Loading feature matrix …")
X = pd.read_csv(feature_matrix, index_col=0)
y = pd.read_csv(labels_file, index_col=0).squeeze().map({"control": 0, "treated": 1})

log.info(f"  Samples: {X.shape[0]} | Features: {X.shape[1]}")
log.info(f"  Class distribution: {y.value_counts().to_dict()}")

# ── Stratified CV ─────────────────────────────────────────
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=ml_cfg["random_state"])


# ── Optuna HPO — XGBoost ─────────────────────────────────
def objective_xgb(trial: optuna.Trial) -> float:
    params = {
        "n_estimators":     trial.suggest_int("n_estimators", 100, 600),
        "max_depth":        trial.suggest_int("max_depth", 3, 10),
        "learning_rate":    trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
        "subsample":        trial.suggest_float("subsample", 0.6, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "reg_alpha":        trial.suggest_float("reg_alpha", 1e-4, 10, log=True),
        "reg_lambda":       trial.suggest_float("reg_lambda", 1e-4, 10, log=True),
        "eval_metric": "logloss",
        "use_label_encoder": False,
        "random_state": ml_cfg["random_state"],
    }
    model = XGBClassifier(**params)
    scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
    return scores.mean()


log.info("Running Optuna HPO for XGBoost (50 trials) …")
study_xgb = optuna.create_study(direction="maximize",
                                 sampler=optuna.samplers.TPESampler(seed=42))
study_xgb.optimize(objective_xgb, n_trials=50, show_progress_bar=False)
best_xgb_params = study_xgb.best_params | {
    "eval_metric": "logloss",
    "use_label_encoder": False,
    "random_state": ml_cfg["random_state"],
}
log.info(f"  Best XGB AUC: {study_xgb.best_value:.4f}")


# ── Train all models ──────────────────────────────────────
mlflow.set_experiment("GenomicsSentinel-ML")

models = {
    "RandomForest": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(
            n_estimators=ml_cfg["n_estimators"],
            random_state=ml_cfg["random_state"],
            n_jobs=-1,
        )),
    ]),
    "XGBoost": XGBClassifier(**best_xgb_params),
    "SVM_RBF": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(kernel="rbf", probability=True,
                    random_state=ml_cfg["random_state"])),
    ]),
}

all_metrics = {}
best_model_name, best_auc, best_clf = None, 0.0, None

for name, clf in models.items():
    with mlflow.start_run(run_name=name):
        auc_scores = cross_val_score(clf, X, y, cv=cv,
                                     scoring="roc_auc", n_jobs=-1)
        mcc_scores = cross_val_score(clf, X, y, cv=cv,
                                     scoring="matthews_corrcoef", n_jobs=-1)

        clf.fit(X, y)

        metrics = {
            "model": name,
            "roc_auc_mean": round(auc_scores.mean(), 4),
            "roc_auc_std":  round(auc_scores.std(), 4),
            "mcc_mean":     round(mcc_scores.mean(), 4),
        }
        all_metrics[name] = metrics

        mlflow.log_metrics({
            "roc_auc": metrics["roc_auc_mean"],
            "mcc":     metrics["mcc_mean"],
        })
        mlflow.sklearn.log_model(clf, artifact_path=name)

        log.info(f"  {name}: AUC = {metrics['roc_auc_mean']:.4f} ± "
                 f"{metrics['roc_auc_std']:.4f} | MCC = {metrics['mcc_mean']:.4f}")

        if metrics["roc_auc_mean"] > best_auc:
            best_auc, best_model_name, best_clf = metrics["roc_auc_mean"], name, clf

# ── ROC Curves ────────────────────────────────────────────
log.info("Generating ROC curve …")
fig, ax = plt.subplots(figsize=(7, 6))
colors = {"RandomForest": "#D85A30", "XGBoost": "#1D9E75", "SVM_RBF": "#378ADD"}

for name, clf in models.items():
    y_prob = clf.predict_proba(X)[:, 1]
    fpr, tpr, _ = roc_curve(y, y_prob)
    auc = all_metrics[name]["roc_auc_mean"]
    ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})",
            color=colors[name], linewidth=2)

ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5)
ax.set_xlabel("False Positive Rate", fontsize=11)
ax.set_ylabel("True Positive Rate", fontsize=11)
ax.set_title("ROC Curves — All Models", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(out_roc, dpi=300)
plt.close()

# ── SHAP Explanation — Best model ────────────────────────
log.info(f"Computing SHAP values for best model: {best_model_name} …")

clf_for_shap = best_clf
if hasattr(best_clf, "named_steps"):
    X_shap = best_clf.named_steps["scaler"].transform(X)
    clf_for_shap = best_clf.named_steps["clf"]
else:
    X_shap = X.values

explainer   = shap.TreeExplainer(clf_for_shap)
shap_values = explainer.shap_values(X_shap)

if isinstance(shap_values, list):
    shap_values = shap_values[1]

fig, ax = plt.subplots(figsize=(10, 8))
shap.summary_plot(
    shap_values, X,
    plot_type="dot",
    max_display=25,
    show=False,
    plot_size=None,
)
ax.set_title(f"SHAP Summary — {best_model_name}", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(out_shap, dpi=300, bbox_inches="tight")
plt.close()

# ── Save metrics ──────────────────────────────────────────
summary = {
    "best_model": best_model_name,
    "best_roc_auc": best_auc,
    "all_models": all_metrics,
    "optuna_xgb_trials": len(study_xgb.trials),
}
with open(out_metrics, "w") as f:
    json.dump(summary, f, indent=2)

log.info("✓ ML analysis complete.")
log.info(f"  Best model: {best_model_name} | AUC: {best_auc:.4f}")
