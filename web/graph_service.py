from __future__ import annotations

from typing import Any, Dict, List, Tuple

import pandas as pd


from .data_loader import add_node_unique, pick_col

from src.graphe.similarite_jaccord_proteine import get_similar_proteines
from src.graphe.plus_longue_sous_chaine import get_similarity_dynamic_blocks
# importer le prédicteur de ligands
try:
    from src.ligands.ligand_predictor import gbm_compatibility_score as ligand_gbm_compatibility
except Exception:
    ligand_gbm_compatibility = None

def build_subgraph(
    data: Dict[str, pd.DataFrame],
    seed: str,
    k: int,
    rels: List[str],
    max_n: int = 200,
) -> Dict[str, Any]:
    """
    Construit un sous-graphe (nodes/edges) autour de seed.
    k est là pour plus tard (BFS multi-hop), pour l'instant on fait 1-hop.
    """
    biogrid = data["biogrid"]
    proteines = data["proteines"] 

    nodes: List[Dict[str, Any]] = []
    edges: List[Dict[str, Any]] = []
    seen: set = set()

    # Alias canonicalisation (ex: GPR30, GPER1 -> GPER)
    alias_map = {
        'GPR30': 'GPER',
        'GPER1': 'GPER',
        'GPER': 'GPER',
        'A8K0M3': 'GPER',
        'Q99527': 'GPER',
        'AAH82766': 'GPER'
    }

    def normalize_id(x: str) -> str:
        if not isinstance(x, str):
            return x
        token = x.strip()
        for a, c in alias_map.items():
            if token.upper() == a.upper():
                return c
        return token

    # Seed node (normalisé)
    seed = normalize_id(seed)
    add_node_unique(nodes, seen, {"id": seed, "label": seed, "shape": "hexagon"})

    # ---- BioGRID: INTERACTS_WITH ----
    if "INTERACTS_WITH" in rels:
        col_a = pick_col(biogrid, ["OFFICIAL_SYMBOL_A"])
        col_b = pick_col(biogrid, ["OFFICIAL_SYMBOL_B"])
        if not col_a or not col_b:
            return {"error": "Colonnes BioGRID non trouvées", "columns": list(biogrid.columns)}

        a_series = biogrid[col_a].astype(str)
        b_series = biogrid[col_b].astype(str)

        df = biogrid[
            a_series.str.contains(seed, case=False, na=False) |
            b_series.str.contains(seed, case=False, na=False)
        ]

        for _, r in df.iterrows():
            a_raw = str(r[col_a])
            b_raw = str(r[col_b])
            a = normalize_id(a_raw)
            b = normalize_id(b_raw)

            add_node_unique(nodes, seen, {
                "id": a, "label": a, "shape": "dot",
                "group": "Protein",
                "title": f"Protein (BioGRID symbol): {a_raw}"
            })
            add_node_unique(nodes, seen, {
                "id": b, "label": b, "shape": "dot",
                "group": "Protein",
                "title": f"Protein (BioGRID symbol): {b_raw}"
            })

            edges.append({
                "from": a,
                "to": b,
                "label": "INTERACTS_WITH",
                "score": None,
                "source": "BioGRID",
                "title": "BioGRID | score: N/A"
            })

    # ---- LIGANDS: BINDS_TO (from merged_ligands.csv) ----
    if "BINDS_TO" in rels:
        lig_df = data.get("ligands")
        if lig_df is not None:
            col_ligand = pick_col(lig_df, ["source_id", "molecule_chembl_id"])
            col_acc = pick_col(lig_df, ["Protein_Accession", "target_id"])

            if col_ligand and col_acc:
                # Filtrer pour la protéine courante (seed)
                mask = lig_df[col_acc].astype(str).apply(normalize_id) == seed
                df_filtered = lig_df[mask].head(max_n)

                col_type = pick_col(lig_df, ["activity_type", "standard_type"])
                col_value = pick_col(lig_df, ["activity_value", "standard_value"])
                col_units = pick_col(lig_df, ["activity_unit", "standard_units"])
                col_origin = pick_col(lig_df, ["origin"])
                col_smiles = pick_col(lig_df, ["canonical_smiles"])

                # Groupement par ligand pour une arête unique
                ligand_groups = {}
                for _, r in df_filtered.iterrows():
                    lig_id = str(r[col_ligand])
                    if lig_id not in ligand_groups:
                        ligand_groups[lig_id] = {
                            "assays": [],
                            "smiles": str(r[col_smiles]) if col_smiles and pd.notna(r[col_smiles]) else "",
                            "origins": set()
                        }
                    
                    stype  = str(r[col_type]) if col_type and pd.notna(r[col_type]) else ""
                    sval   = str(r[col_value]) if col_value and pd.notna(r[col_value]) else ""
                    sunit  = str(r[col_units]) if col_units and pd.notna(r[col_units]) else ""
                    origin = str(r[col_origin]) if col_origin and pd.notna(r[col_origin]) else "N/A"
                    
                    if stype or sval:
                        val_str = f"{stype} = {sval} {sunit}".strip()
                        ligand_groups[lig_id]["assays"].append(f"[{origin}] {val_str}")
                    
                    ligand_groups[lig_id]["origins"].add(origin)

                for lig_id, info in ligand_groups.items():
                    add_node_unique(nodes, seen, {"id": lig_id, "label": lig_id, "shape": "diamond", "group": "Ligand"})
                    
                    tooltipLines = [f"Ligand: {lig_id}"]
                    if info["assays"]:
                        tooltipLines.append("-" * 20)
                        tooltipLines.extend(info["assays"][:10])
                        if len(info["assays"]) > 10:
                            tooltipLines.append(f"... and {len(info['assays'])-10} more")
                    
                    if info["smiles"]:
                        tooltipLines.append("-" * 20)
                        tooltipLines.append(f"SMILES: {info['smiles'][:100]}...")

                    edges.append({
                        "from": lig_id,
                        "to": seed,
                        "label": "BINDS_TO",
                        "source": ", ".join(info["origins"]),
                        "title": "\n".join(tooltipLines)
                    })

    # ---- Prediction placeholder: PREDICTED_BINDS ----
    if "PREDICTED_BINDS" in rels:
        # 1. Prédiction basée sur les scores GBM (calculate_jaccard_scores)
        predicted_edges = calculate_jaccard_scores(data, seed)
        
        # 2. On prend le Top N pour ne pas surcharger le graphe
        top_predictions = predicted_edges[:10] 
        
        # 3. On les ajoute au graphe final
        seen_pred_edges: set = set()
        for pred in top_predictions:
            target_id = str(pred["target"])
            edge_key = (pred["source"], target_id)
            if edge_key in seen_pred_edges:
                continue
            seen_pred_edges.add(edge_key)

            # Ajouter le nœud cible s'il n'est pas déjà dans le graphe
            add_node_unique(nodes, seen, {
                "id": target_id, 
                "label": target_id, 
                "shape": "diamond", 
                "group": "Ligand",
                "color": "#4BD1C3",
            }) 
            
            # afficher le score directement sur l'arête et la colorer (Solution 1)
            sc = float(pred.get("score", 0.0) or 0.0) *100
            edge_label = f"{sc:.2f}%"
            # Couleur = Solution 1 (vert de la légende)
            edge_color = "#4BD1C3"
            edges.append({
                "from": pred["source"],
                "to": target_id,
                "label": edge_label,
                "score": sc,
                "color": {"color": edge_color},
                "title": f"PREDICTED_BINDS | Score: {sc:.3f} | voisins communs: {pred['common_neighbors']}"
            })

        # 4. Ajouter les nouveaux ligands issus du fichier de prédiction séquence complète (Solution 2)
        # ---- Solution 2 : nouveaux ligands prédits ----
        # ---- Solution 2 : nouveaux ligands prédits (GBM) ----
        df_new = data.get("predicted_ligands_full")
        if df_new is not None and not df_new.empty and ligand_gbm_compatibility:
            # Colonnes possibles pour ID et SMILES
            id_col_candidates = ["Molecule ChEMBL ID", "chembl_ID", "Compound Key", "source_id"]
            smiles_col_candidates = ["Smiles", "smiles", "canonical_smiles"]

            id_col = next((c for c in id_col_candidates if c in df_new.columns), None)
            smiles_col = next((c for c in smiles_col_candidates if c in df_new.columns), None)

            if not id_col or not smiles_col:
                print("Cols ID ou SMILES non trouvées dans predicted_ligands_full")
            else:
                used_ids: set = set()
                count = 0
                print("=== Solution 2: Top 10 Predicted Ligands (GBM) ===")
                
                for _, r in df_new.iterrows():
                    if count >= 10:  # limiter à 10 ligands
                        break

                    lig_id = str(r[id_col]).strip()
                    if not lig_id or lig_id in used_ids:
                        continue
                    used_ids.add(lig_id)

                    smiles = str(r[smiles_col]).strip()
                    score = 0.0

                    if smiles:
                        try:
                            # Appel du prédicteur GBM
                            res = ligand_gbm_compatibility(smiles)
                            score = float(res)  # résultat attendu entre 0 et 100
                        except Exception as e:
                            print(f"Erreur GBM pour {lig_id}: {e}")
                            score = 0.0

                    # Print debug
                    print(f"{lig_id} | GBM score: {score:.2f}% | SMILES: {smiles[:50]}...")

                    # Ajouter le nœud ligand (orange)
                    add_node_unique(
                        nodes,
                        seen,
                        {
                            "id": lig_id,
                            "label": lig_id,
                            "shape": "diamond",
                            "group": "Ligand",
                            "color": "#EF9C42",
                        },
                    )

                    # Ajouter l'arête seed -> ligand
                    edges.append(
                        {
                            "from": seed,
                            "to": lig_id,
                            "label": f"{score:.2f}%",  # en pourcentage
                            "score": score,
                            "color": {"color": "#EF9C42"},
                            "title": f"PREDICTED_LIGAND | GBM score: {score:.2f}% | SMILES: {smiles[:100]}...",
                        }
                    )

                    count += 1
                    
            # ---- TCGA : EXPRESSED_IN_GBM ----
    if "EXPRESSED_IN_GBM" in rels:
        tcga_data = data.get("tcga")
        if tcga_data is not None and not tcga_data.empty:
            mean_age = tcga_data["age_at_initial_pathologic_diagnosis"].mean()
            gender_counts = tcga_data["gender"].value_counts().to_dict()

            add_node_unique(
                nodes,
                seen,
                {
                    "id": "GBM",
                    "label": "Glioblastoma",
                    "shape": "triangle",
                    "color": "#e74c3c",
                    "clinical_info": {
                        "cohort_size": len(tcga_data),
                        "mean_age": round(mean_age, 1),
                        "male_count": gender_counts.get("MALE", 0),
                        "female_count": gender_counts.get("FEMALE", 0),
                    },
                },
            )

            current_node_ids = {n["id"] for n in nodes}
            for gene in ["GPER", "ESR1", "ESR2"]:
                gene_norm = normalize_id(gene)
                if gene_norm in current_node_ids and gene in tcga_data.columns:
                    mean_expr = tcga_data[gene].mean()
                    edges.append(
                        {
                            "id": f"tcga_{gene}_GBM",
                            "from": gene_norm,
                            "to": "GBM",
                            "label": f"{mean_expr:.2f}",
                            "title": f"Expression moyenne TCGA ({gene}): {mean_expr:.2f}",
                            "color": {"color": "#e74c3c"},
                        }
                    )

    # Vue détaillée : quelques patients GBM très exprimés pour GPER
    if "GBM_PATIENTS" in rels:
        tcga_data = data.get("tcga")

        if tcga_data is not None and not tcga_data.empty and seed in ["GPER", "GPER1"]:
            sample_id_col = tcga_data.columns[0]
            top_patients = (
                tcga_data.sort_values(by="GPER", ascending=False)
                .head(min(10, len(tcga_data)))
            )

            for _, row in top_patients.iterrows():
                sample_id = str(row[sample_id_col])
                expr_val = float(row["GPER"])
                age = row.get("age_at_initial_pathologic_diagnosis", None)
                gender = row.get("gender", None)
                sample_type = row.get("sample_type", None)

                add_node_unique(
                    nodes,
                    seen,
                    {
                        "id": sample_id,
                        "label": sample_id,
                        "group": "GBM_sample",
                        "shape": "box",
                    },
                )

                title_parts = [f"Expression GPER: {expr_val:.2f}"]
                if age is not None:
                    title_parts.append(f"Âge: {age}")
                if gender:
                    title_parts.append(f"Sexe: {gender}")
                if sample_type:
                    title_parts.append(f"Type: {sample_type}")

                edges.append(
                    {
                        "from": seed,
                        "to": sample_id,
                        # label = quantité de GPER pour ce patient
                        "label": f"{expr_val:.2f}",
                        "value": expr_val,
                        "title": " | ".join(title_parts),
                        "color": {"color": "#e74c3c"},
                    }
                )

    if "EXPRESSED_IN_BRAIN" in rels:
        gtex_data = data.get("gtex")
        if gtex_data is not None and not gtex_data.empty:
            # Colonnes : Expression_TPM, Brain_Region, Type, SampleID
            mean_expr = gtex_data["Expression_TPM"].mean()

            # Nœud résumé "HEALTHY_BRAIN"
            add_node_unique(
                nodes,
                seen,
                {
                    "id": "HEALTHY_BRAIN",
                    "label": "Cerveau sain",
                    "shape": "triangle",
                    "color": "#2ecc71",
                    "clinical_info": {
                        "cohort_size": len(gtex_data),
                        "mean_expr": round(mean_expr, 2),
                    },
                },
            )

            # Lien moyen seed -> cerveau sain (on considère le fichier comme expression de GPER)
            edges.append(
                {
                    "id": f"gtex_{seed}_BRAIN",
                    "from": seed,
                    "to": "HEALTHY_BRAIN",
                    "label": f"{mean_expr:.2f}",
                    "title": f"Expression moyenne (GTEx, GPER): {mean_expr:.2f} TPM",
                    "color": {"color": "#2ecc71"},
                }
            )

    if "SEQUENCE_SIMILARITY" in rels:

        similar_proteins = get_similar_proteines(proteines, seed_proteine=seed, top_n=20)
        
        for rank, item in enumerate(similar_proteins):
            source_id = str(item["source_id"])
            target_gene = str(item["gene"])
            # Label d'affichage : gène si disponible, sinon source_id
            label = target_gene if target_gene and target_gene != "nan" else source_id
            label = normalize_id(label)

            score = item["score"]
            
            # --- Choix visuels ---
            color = "#2ecc71" if score > 0.5 else "#1ABC9C"
            
            # 2. Ajout du Nœud — on utilise source_id comme identifiant ET label unique
            add_node_unique(nodes, seen, {
                "id": source_id,
                "label": source_id,
                "shape": "hexagon",   # Forme spécifique GPCR
                "color": color,
                "size": 20 + (score * 20), # Taille proportionnelle au score
                "title": f"Protein: {source_id} ({label})\nSequence Similarity: {score:.4f}"
            })
            
            # 3. Ajout du Lien
            # Ignorer les auto-liens (self-loop) — pas utile visuellement
            if source_id == seed or label == seed:
                continue

            # Les 3 premiers (meilleurs scores) ont un trait épais
            edge_width = 6 if rank < 3 else 1.5

            edges.append({
                "from": seed,
                "to": source_id,
                "label": f"{score * 100:.2f}%",
                "width": edge_width,
                "color": {"color": color, "opacity": 0.8},
                "title": f"Jaccard Score: {score * 100:.2f}%"
            })

    if "LCS_SIMILARITY" in rels:

        lcs_proteins = get_similarity_dynamic_blocks(proteines, seed_proteine=seed, top_n=20, min_block_size=4)

        for rank, item in enumerate(lcs_proteins):
            source_id = str(item["source_id"])
            target_gene = str(item["gene"])

            score_pct = item["score_percent"]  # valeur en %
            nb_blocks = item.get("nb_blocks", 0)
            longest = item.get("longest_block", 0)

            # Gradient de couleur selon le score : orange 
            color = "#e67e22" 

            # Ajout du nœud — identifié ET affiché par source_id
            add_node_unique(nodes, seen, {
                "id": source_id,
                "label": source_id,
                "shape": "hexagon",
                "color": color,
                "size": 15 + (score_pct * 0.3),
                "title": f"Protein: {source_id} ({target_gene})\nLCS Score: {score_pct:.2f}%\nBlocs: {nb_blocks} | Plus long: {longest} aa"
            })

            # Ignorer les auto-liens
            if source_id == seed:
                continue

            # Les 3 premiers (meilleurs scores) ont un trait épais
            edge_width = 6 if rank < 3 else (5 + (score_pct * 0.15) if score_pct > 20 else 1.5)
            edges.append({
                "from": seed,
                "to": source_id,
                "label": f"{score_pct:.2f}%",
                "width": edge_width,
                "color": {"color": color, "opacity": 0.8},
                "title": f"LCS Score: {score_pct:.2f}% | Blocs communs: {nb_blocks} | Plus long: {longest} aa"
            })

    return {"nodes": nodes, "edges": edges}


def resolve_entity(data: Dict[str, pd.DataFrame], node_id: str) -> Tuple[Dict[str, Any], int]:
    """
    Renvoie (payload, status_code)
    - Essaie merged_proteines.csv
    - Puis merged_ligands.csv
    - Puis chembl_gper_ligands.csv (molecule_chembl_id)
    - Sinon: fiche minimale (200) au lieu de 404 (UX)
    """
    proteines = data["proteines"]
    ligands   = data["ligands"]
    chembl    = data["chembl"]

    node_id = str(node_id)

    # Proteins: match strict puis contains sur toutes colonnes
    for col in proteines.columns:
        s = proteines[col].astype(str)
        hit = proteines[s.str.fullmatch(node_id, case=False, na=False)]
        if len(hit):
            return hit.iloc[0].to_dict(), 200

    for col in proteines.columns:
        s = proteines[col].astype(str)
        hit = proteines[s.str.contains(node_id, case=False, na=False)]
        if len(hit):
            return hit.iloc[0].to_dict(), 200

    # Ligands
    for col in ligands.columns:
        s = ligands[col].astype(str)
        hit = ligands[s.str.fullmatch(node_id, case=False, na=False)]
        if len(hit):
            return hit.iloc[0].to_dict(), 200

    for col in ligands.columns:
        s = ligands[col].astype(str)
        hit = ligands[s.str.contains(node_id, case=False, na=False)]
        if len(hit):
            return hit.iloc[0].to_dict(), 200

    # ChEMBL fallback
    if "molecule_chembl_id" in chembl.columns:
        df = chembl[chembl["molecule_chembl_id"].astype(str) == node_id]
        if len(df):
            return df.iloc[0].to_dict(), 200

    # Fiche minimale (au lieu de 404)
    return {
        "id": node_id,
        "type": "Entity",
        "source": "Graph node",
        "note": "Pas de détails trouvés dans merged_proteines / merged_ligands / chembl.",
    }, 200


def calculate_jaccard_scores(data: Dict[str, pd.DataFrame], seed: str) -> List[Dict[str, Any]]:
    """Fallback prediction: return ligands scored by GBM compatibility.

    Returns a list of dicts with keys: source (protein), target (ligand id), score, common_neighbors.
    Uses `data['ligands_scored']` when available, otherwise `data['ligands']`.
    """
    df = data.get("ligands_scored") if data.get("ligands_scored") is not None else data.get("ligands")
    results: List[Dict[str, Any]] = []
    if df is None or df.empty:
        return results

    # detect smiles column
    smiles_col = pick_col(df, ["canonical_smiles", "Smiles", "smiles", "smiles_raw"])
    id_col = pick_col(df, ["Molecule ChEMBL ID", "chembl_ID", "Compound Key", "source_id"])
    if not id_col:
        id_col = df.columns[0]

    has_gbm_col = "gbm_score" in df.columns
    
    for _, r in df.iterrows():
        lig_id = str(r.get(id_col, '')).strip() or None
        if not lig_id:
            continue
            
        score = None
        # Priorité 1 : colonne déjà présente (cas du fichier _scored.csv)
        if has_gbm_col:
            try:
                val = r.get("gbm_score")
                if pd.notna(val) and val != "":
                    score = float(val)
            except (ValueError, TypeError):
                pass
        
        # Priorité 2 : calcul si nécessaire
        if score is None and ligand_gbm_compatibility:
            smiles = str(r.get(smiles_col, '')).strip() if smiles_col else ''
            if smiles:
                try:
                    res = ligand_gbm_compatibility(smiles)
                    score = res.get('gbm_score')
                except Exception:
                    score = None

        if lig_id:
            results.append({
                'source': seed,
                'target': lig_id,
                'score': float(score) if score is not None else 0.0,
                'common_neighbors': 0
            })

    # sort desc by score
    results.sort(key=lambda x: x.get('score', 0.0), reverse=True)
    return results
