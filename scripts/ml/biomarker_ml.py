# ================================================================
# Découverte de biomarqueurs — Apprentissage automatique
# Auteur  : Driss El Oifi
# Date    : mai 2026
# Dataset : GSE152418 — COVID-19 Severe vs Healthy
# ================================================================
# Approche :
#   Les gènes différentiellement exprimés (DESeq2) servent
#   de features pour entraîner des modèles de classification.
#
#   On cherche le sous-ensemble minimal de gènes qui prédit
#   le mieux la sévérité COVID — ces gènes sont les biomarqueurs.
#
#   Pipeline :
#   1. Random Forest — robuste, interprétable via importance
#   2. XGBoost — gradient boosting, souvent plus précis
#   3. SVM RBF — bon sur données de haute dimension
#   4. Optuna — optimisation automatique des hyperparamètres
#   5. SHAP — explication des prédictions (pas une boîte noire)
#
#   Métrique principale : AUC-ROC (robuste aux classes déséquilibrées)
#   Métrique secondaire : MCC (Matthews Correlation Coefficient)
#
# Référence SHAP : Lundberg & Lee, NeurIPS 2017
# Référence Optuna : Akiba et al., KDD 2019
# ================================================================

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

logging.basicConfig(
    level  = logging.INFO,
    format = "%(levelname)s | %(message)s"
)
log = logging.getLogger("biomarker-ml")

# ----------------------------------------------------------------
# Chargement des données
# ----------------------------------------------------------------
# feature_matrix : gènes DE comme colonnes, échantillons comme lignes
# labels         : 0 = Healthy, 1 = Severe
# ----------------------------------------------------------------
feature_matrix = snakemake.input["features"]
labels_file    = snakemake.input["labels"]
out_metrics    = snakemake.output["metrics"]
out_shap       = snakemake.output["shap_plot"]
out_roc        = snakemake.output["roc_plot"]

cfg    = snakemake.config
ml_cfg = cfg["ml"]

log.info("Chargement de la matrice de features ...")
X = pd.read_csv(feature_matrix, index_col=0)
y = pd.read_csv(labels_file, index_col=0).squeeze().map(
    {"Healthy": 0, "Severe": 1}
)

log.info(f"  Échantillons : {X.shape[0]}")
log.info(f"  Features     : {X.shape[1]}")
log.info(f"  Distribution : {y.value_counts().to_dict()}")

# ----------------------------------------------------------------
# Validation croisée stratifiée
# ----------------------------------------------------------------
# Stratifiée = chaque fold conserve la proportion Healthy/Severe.
# Important ici car les classes sont déséquilibrées (17 vs 8).
# 5 folds = bon compromis biais/variance avec n=25 échantillons.
# ----------------------------------------------------------------
cv = StratifiedKFold(
    n_splits   = 5,
    shuffle    = True,
    random_state = ml_cfg["random_state"]
)

# ----------------------------------------------------------------
# Optimisation XGBoost avec Optuna
# ----------------------------------------------------------------
# Optuna utilise l'algorithme TPE (Tree-structured Parzen Estimator)
# pour explorer intelligemment l'espace des hyperparamètres.
# Beaucoup plus efficace qu'un GridSearch exhaustif.
# ----------------------------------------------------------------
def objectif_xgb(trial: optuna.Trial) -> float:
    params = {
        "n_estimators"     : trial.suggest_int("n_estimators", 100, 600),
        "max_depth"        : trial.suggest_int("max_depth", 3, 10),
        "learning_rate"    : trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
        "subsample"        : trial.suggest_float("subsample", 0.6, 1.0),
        "colsample_bytree" : trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "reg_alpha"        : trial.suggest_float("reg_alpha", 1e-4, 10, log=True),
        "reg_lambda"       : trial.suggest_float("reg_lambda", 1e-4, 10, log=True),
        "eval_metric"      : "logloss",
        "use_label_encoder": False,
        "random_state"     : ml_cfg["random_state"],
    }
    modele = XGBClassifier(**params)
    scores = cross_val_score(
        modele, X, y, cv=cv, scoring="roc_auc", n_jobs=-1
    )
    return scores.mean()

log.info("Optimisation Optuna XGBoost (50 essais) ...")
etude = optuna.create_study(
    direction = "maximize",
    sampler   = optuna.samplers.TPESampler(seed=42)
)
etude.optimize(objectif_xgb, n_trials=50, show_progress_bar=False)

meilleurs_params = etude.best_params | {
    "eval_metric"      : "logloss",
    "use_label_encoder": False,
    "random_state"     : ml_cfg["random_state"],
}
log.info(f"  Meilleur AUC XGBoost : {etude.best_value:.4f}")

# ----------------------------------------------------------------
# Entraînement des 3 modèles
# ----------------------------------------------------------------
modeles = {
    "RandomForest": Pipeline([
        ("normalisation", StandardScaler()),
        ("clf", RandomForestClassifier(
            n_estimators = ml_cfg["n_estimators"],
            random_state = ml_cfg["random_state"],
            n_jobs       = -1,
        )),
    ]),
    "XGBoost": XGBClassifier(**meilleurs_params),
    "SVM_RBF": Pipeline([
        ("normalisation", StandardScaler()),
        ("clf", SVC(
            kernel       = "rbf",
            probability  = True,
            random_state = ml_cfg["random_state"]
        )),
    ]),
}

toutes_metriques = {}
meilleur_nom, meilleur_auc, meilleur_clf = None, 0.0, None

mlflow.set_experiment("GenomicsSentinel-COVID")

for nom, clf in modeles.items():
    with mlflow.start_run(run_name=nom):
        auc_scores = cross_val_score(
            clf, X, y, cv=cv, scoring="roc_auc", n_jobs=-1
        )
        mcc_scores = cross_val_score(
            clf, X, y, cv=cv, scoring="matthews_corrcoef", n_jobs=-1
        )
        clf.fit(X, y)

        metriques = {
            "modele"       : nom,
            "roc_auc_mean" : round(auc_scores.mean(), 4),
            "roc_auc_std"  : round(auc_scores.std(),  4),
            "mcc_mean"     : round(mcc_scores.mean(), 4),
        }
        toutes_metriques[nom] = metriques

        mlflow.log_metrics({
            "roc_auc" : metriques["roc_auc_mean"],
            "mcc"     : metriques["mcc_mean"],
        })
        mlflow.sklearn.log_model(clf, artifact_path=nom)

        log.info(
            f"  {nom:15s} AUC = {metriques['roc_auc_mean']:.4f} "
            f"± {metriques['roc_auc_std']:.4f} | "
            f"MCC = {metriques['mcc_mean']:.4f}"
        )

        if metriques["roc_auc_mean"] > meilleur_auc:
            meilleur_auc = metriques["roc_auc_mean"]
            meilleur_nom = nom
            meilleur_clf = clf

# ----------------------------------------------------------------
# Courbes ROC
# ----------------------------------------------------------------
couleurs = {
    "RandomForest": "#D85A30",
    "XGBoost"     : "#1D9E75",
    "SVM_RBF"     : "#378ADD"
}

fig, ax = plt.subplots(figsize=(7, 6))
for nom, clf in modeles.items():
    y_prob = clf.predict_proba(X)[:, 1]
    fpr, tpr, _ = roc_curve(y, y_prob)
    auc = toutes_metriques[nom]["roc_auc_mean"]
    ax.plot(fpr, tpr,
            label     = f"{nom} (AUC={auc:.3f})",
            color     = couleurs[nom],
            linewidth = 2)

ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5)
ax.set_xlabel("Taux faux positifs", fontsize=11)
ax.set_ylabel("Taux vrais positifs", fontsize=11)
ax.set_title("Courbes ROC — COVID Severe vs Healthy",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(out_roc, dpi=300)
plt.close()

# ----------------------------------------------------------------
# Valeurs SHAP — explication du meilleur modèle
# ----------------------------------------------------------------
# SHAP (SHapley Additive exPlanations) décompose la prédiction
# de chaque échantillon en contributions par feature.
#
# Interprétation : si SHAP(IFI27, patient_5) = +0.8
# → IFI27 pousse la prédiction vers "Severe" pour ce patient
# ----------------------------------------------------------------
log.info(f"Calcul des valeurs SHAP — modèle : {meilleur_nom} ...")

clf_shap = meilleur_clf
if hasattr(meilleur_clf, "named_steps"):
    X_shap   = meilleur_clf.named_steps["normalisation"].transform(X)
    clf_shap = meilleur_clf.named_steps["clf"]
else:
    X_shap = X.values

explicateur  = shap.TreeExplainer(clf_shap)
valeurs_shap = explicateur.shap_values(X_shap)

if isinstance(valeurs_shap, list):
    valeurs_shap = valeurs_shap[1]

fig, ax = plt.subplots(figsize=(10, 8))
shap.summary_plot(
    valeurs_shap, X,
    plot_type   = "dot",
    max_display = 25,
    show        = False,
    plot_size   = None,
)
ax.set_title(
    f"Importance SHAP — {meilleur_nom}\nCOVID Severe vs Healthy",
    fontsize=13, fontweight="bold"
)
plt.tight_layout()
plt.savefig(out_shap, dpi=300, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------
# Export des métriques
# ----------------------------------------------------------------
resume = {
    "meilleur_modele" : meilleur_nom,
    "meilleur_auc"    : meilleur_auc,
    "tous_modeles"    : toutes_metriques,
    "optuna_essais"   : len(etude.trials),
}
with open(out_metrics, "w") as f:
    json.dump(resume, f, indent=2, ensure_ascii=False)

log.info(f"✓ Analyse ML terminée — meilleur modèle : {meilleur_nom} (AUC={meilleur_auc:.4f})")