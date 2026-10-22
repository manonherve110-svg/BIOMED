from __future__ import annotations

from typing import Any, Dict, List


def get_graph_modes() -> List[Dict[str, Any]]:
    return [
        
        {
            "id": "known_ligands",
            "title": "Ligands du GPER connus",
            "desc": "LIguands de la base de donnée chembl / pubchem",
            "params": {"seed": "GPER", "k": 1, "rels": "BINDS_TO"},
        },
        {
            "id": "prediction_GPCR",
            "title": "Similarité de séquence (GPCR) -- Solution 1",
            "desc": "Graphe basé sur la similarité Jaccard des séquences d'acides aminés.",
            "params": {
                "seed": "GPER", 
                "k": 1, 
                "rels": "SEQUENCE_SIMILARITY" # <--- C'est ce mot-clé qui déclenchera le calcul
                },
        },
        {
            "id": "lcs_GPCR",
            "title": "Similarité de séquence (GPCR) -- Solution 2",
            "desc": "Graphe basé sur la plus longue sous-chaîne commune (blocs dynamiques).",
            "params": {
                "seed": "GPER",
                "k": 1,
                "rels": "LCS_SIMILARITY"
            },
        },
        {
            "id": "predicted_ligands",
            "title": "Prédiction ligands",
            "desc": "Arêtes PREDICTED_BINDS avec score = ...",
            "params": {"seed": "GPER", "k": 1, "rels": "PREDICTED_BINDS"},
        }
    ]
