import pandas as pd
import os
from difflib import SequenceMatcher

# Chargement des données
relative_path = "data/processed/merged_proteines.csv"
absolute_path = os.path.abspath(relative_path)
df = pd.read_csv(absolute_path)

def get_similarity_dynamic_blocks(df, seed_proteine, top_n=20, min_block_size=4):
    """
    Calcule la similarité en trouvant itérativement les plus longs blocs communs,
    quelle que soit leur position et leur taille.

    min_block_size : Taille minimale pour qu'un bloc soit considéré (4 est bien pour éviter le bruit).
    """

    # --- 1. Gestion des Alias ---
    alias_map = {
        'GPR30': 'GPER', 'GPER1': 'GPER', 'GPER': 'GPER',
        'A8K0M3': 'GPER', 'A8KM03': 'GPER'
    }
    seed_key = alias_map.get(seed_proteine.upper(), seed_proteine)
    target_aliases = [k for k,v in alias_map.items() if v == seed_key]
    if seed_key not in target_aliases:
        target_aliases.append(seed_key)

    df_cible = df[df['gene'].isin(target_aliases)]
    if df_cible.empty:
        df_cible = df[df['source_id'] == seed_proteine]

    if df_cible.empty:
        print(f"Erreur : '{seed_proteine}' introuvable.")
        return []

    seq_cible_originale = df_cible["sequence"].iloc[0]
    len_cible = len(seq_cible_originale)
    
    ids_to_exclude = set(df_cible["source_id"].unique())
    ids_to_exclude.update(target_aliases)

    results = []

    print(f"Recherche dynamique de blocs pour {seed_key} ({len_cible} aa)...")
    count = 0 

    # --- 2. Boucle de comparaison ---
    for index, row in df.iterrows():
        
        # Filtres de base
        if row['source_id'] in ids_to_exclude: continue
        if row['gene'] == seed_key: continue
            
        seq_compare = row['sequence']
        if not isinstance(seq_compare, str): continue
        if seq_compare == seq_cible_originale: continue

        s1 = seq_cible_originale
        s2 = seq_compare
        
        total_match_length = 0
        blocks_found = []
        
        while True:
            matcher = SequenceMatcher(None, s1, s2, autojunk=False)
            match = matcher.find_longest_match(0, len(s1), 0, len(s2))
            
            if match.size < min_block_size:
                break
                
            total_match_length += match.size
            motif = s1[match.a : match.a + match.size]
            blocks_found.append(motif)
            
            # On enlève ce motif pour la prochaine itération
            s1 = s1.replace(motif, "", 1)
            s2 = s2.replace(motif, "", 1)
            
            if len(s1) < min_block_size or len(s2) < min_block_size:
                break
        
        # Calcul du pourcentage de ressemblance
        score_percent = (total_match_length / len_cible) * 100
        
        if score_percent > 10.0:  # seuil minimum %
            results.append({
                "gene": row['gene'],
                "source_id": row['source_id'],
                "score_percent": round(score_percent, 2),  
                "nb_blocks": len(blocks_found),
                "longest_block": len(blocks_found[0]) if blocks_found else 0
            })
            
        count += 1
        if count % 1000 == 0:
            print(f"   ... {count} protéines traitées")

    results = sorted(results, key=lambda x: x['score_percent'], reverse=True)
    return results[:top_n]

# --- TEST ---
#print(get_similarity_dynamic_blocks(df, "GPER1", 20, min_block_size=4))