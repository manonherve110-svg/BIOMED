from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_csv(*parts: str, sep: str = ",") -> pd.DataFrame:
    """Charge un CSV avec séparateur configurable."""
    path = os.path.join(BASE_DIR, *parts)
    print("Chargement :", path)
    return pd.read_csv(path, sep=sep)


def pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def add_node_unique(nodes: List[Dict[str, Any]], seen: set, node: Dict[str, Any]) -> None:
    nid = node["id"]
    if nid not in seen:
        nodes.append(node)
        seen.add(nid)


def load_all_data() -> Dict[str, pd.DataFrame]:
    """
    Centralise tous les chargements CSV.
    Ajuste ici si tu changes l'arborescence.
    """
    # Données fusionnées
    proteines = load_csv("data", "processed", "merged_proteines.csv")
    ligands   = load_csv("data", "processed", "merged_ligands.csv")

    # Optionnel : ligands déjà scorés pour le GBM
    ligands_scored_path = os.path.join(BASE_DIR, "data", "processed", "ligands_GPCR_scored.csv")
    ligands_scored = None
    if os.path.exists(ligands_scored_path):
        ligands_scored = pd.read_csv(ligands_scored_path, dtype=str)

    signaling = load_csv("data", "processed", "merged_signaling.csv")
    
    # Données brutes ligands
    chembl    = load_csv("data", "raw", "ligands", "chembl_gper_ligands.csv")
    pubchem   = load_csv("data", "raw", "ligands", "pubchem_gper_ligands.csv")
    
    # Données brutes protéines
    uniprot   = load_csv("data", "raw", "proteines", "uniprot_GPER_GPCR.csv", sep=";")
    alphafold = load_csv("data", "raw", "proteines", "AlphaFold_GPER_GPCR.csv", sep=";")
    
    # Expression (GBM / cerveau sain)
    tcga      = load_csv("data", "raw", "expression", "tcga_gbm_gper_clinical.csv")
    gtex      = load_csv("data", "raw", "expression", "gtex_brain_gper_expression.csv")
    
    # Prédictions de nouveaux ligands (séquence complète)
    predicted_ligands_full = load_csv("data", "processed", "prediction_ligands_chaine_entiere.csv")
    
    # Voies de signalisation
    biogrid   = load_csv("data", "raw", "voie_signalisation", "biogrid_gper_interaction.csv")
    kegg_path = load_csv("data", "raw", "voie_signalisation", "kegg_gper_pathways.csv")

    return {
        "proteines": proteines,
        "ligands": ligands,
        "ligands_scored": ligands_scored,
        "signaling": signaling,
        "chembl": chembl,
        "pubchem": pubchem,
        "uniprot": uniprot,
        "alphafold": alphafold,
        "tcga": tcga,
        "gtex": gtex,
        "biogrid": biogrid,
        "kegg": kegg_path,
        "predicted_ligands_full": predicted_ligands_full,
    }
