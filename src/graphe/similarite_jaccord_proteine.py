import pandas as pd 
import os
from collections import Counter  # Import nécessaire pour gérer les doublons

# chargement des données -----------------------------------------------------------------------------------------------
relative_path = "data/processed/merged_proteines.csv"
absolute_path = os.path.abspath(relative_path)
df = pd.read_csv(absolute_path)

def coupage_seq(sequence, k):
    """
    Découpe une séquence d'acides aminés en sous-séquences de longueur k.
    """
    liste = []
    n = len(str(sequence)) # str() pour éviter les erreurs si NaN
    for i in range(n - k + 1):
        subseq = sequence[i:i+k]  
        liste.append(subseq)
    return liste

def sequence_similarity(seq1, seq2, k=3):
    """
    Calcule la similarité de Jaccard pondérée (avec doublons) 
    entre deux séquences en utilisant Counter.
    """
    list1 = coupage_seq(seq1, k)
    list2 = coupage_seq(seq2, k)
    
    if not list1 or not list2:
        return 0.0
        
    # On compte le nombre d'occurrences de chaque k-mer
    counts1 = Counter(list1)
    counts2 = Counter(list2)
    
    # Intersection : on prend le minimum des occurrences pour chaque k-mer commun
    intersection = counts1 & counts2
    intersection_sum = sum(intersection.values())
    
    # Union : on prend le maximum des occurrences pour chaque k-mer présent
    union = counts1 | counts2
    union_sum = sum(union.values())
    
    return intersection_sum / union_sum





def get_similar_proteines(df, seed_proteine, top_n):

    # Alias connus pour GPER
    alias_map = {
        'GPR30': 'GPER',
        'GPER1': 'GPER',
        'GPER': 'GPER',
        'A8K0M3': 'GPER' 
    }

    seed_key = alias_map.get(seed_proteine.upper(), seed_proteine)

    target_aliases = [k for k,v in alias_map.items() if v == seed_key]
    
    # On s'assure que le seed lui-même est dans la liste
    if seed_key not in target_aliases:
        target_aliases.append(seed_key)

    # 3. Récupération des lignes cibles pour avoir la séquence
    df_cible = df[df['gene'].isin(target_aliases)]

    if df_cible.empty:
        df_cible = df[df['source_id'] == seed_proteine]

    if df_cible.empty:
        print(f"{seed_proteine} non trouvé.")
        return []

    # On prend la première séquence trouvée comme référence
    seq_cible = df_cible["sequence"].iloc[0]

    ids_to_exclude = set(df_cible["source_id"].unique())
    
    ids_to_exclude.update(target_aliases) 
    # ------------------------------

    results = []

    for index, row in df.iterrows():
        
        if row['source_id'] in ids_to_exclude:
            continue
            
        if row['gene'] == seed_key:
            continue

        seq = row['sequence']
        
        if seq == seq_cible:
            continue

        score = sequence_similarity(seq_cible, seq, 3)

        if score >= 0.05:
            results.append({
                "gene": row['gene'],
                "source_id": row['source_id'],
                "score": score
            })

    results = sorted(results, key=lambda x: x['score'], reverse=True)
    return results[:top_n]
#print(get_similar_proteines(df, "GPER1", 20))
