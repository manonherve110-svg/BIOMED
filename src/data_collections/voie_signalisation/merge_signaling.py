import pandas as pd
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(SCRIPT_DIR, "../../../data/raw/voie_signalisation")
PROCESSED_DIR = os.path.join(SCRIPT_DIR, "../../../data/processed")
OUTPUT_FILE = os.path.join(PROCESSED_DIR, "merged_signaling.csv")

def main():
    print("Chargement des données de signalisation...")
    
    # 1. Chargement des fichiers
    try:
        df_biogrid = pd.read_csv(os.path.join(RAW_DIR, "biogrid_gper_interaction.csv"))
        df_kegg_pathways = pd.read_csv(os.path.join(RAW_DIR, "kegg_gper_pathways.csv"))
        df_kegg_gpcrs = pd.read_csv(os.path.join(RAW_DIR, "kegg_all_gpcrs.csv"))
    except FileNotFoundError as e:
        print(f" Erreur : Fichier manquant - {e}")
        return

    # 2. Traitement BioGRID (Interactions physiques)
    print(" Traitement des interactions BioGRID...")
    # On normalise pour avoir une structure (GeneA, GeneB, Type, Source)
    biogrid_final = df_biogrid.rename(columns={
        'OFFICIAL_SYMBOL_A': 'entity_a',
        'OFFICIAL_SYMBOL_B': 'entity_b',
        'EXPERIMENTAL_SYSTEM': 'relation_type'
    })
    biogrid_final['source'] = 'BioGRID'
    biogrid_final = biogrid_final[['entity_a', 'entity_b', 'relation_type', 'source', 'PUBMED_ID']]

    # 3. Traitement KEGG Pathways (Appartenance aux voies)
    print("🗺️ Traitement des voies KEGG...")
    # On sait que hsa:2852 est GPER1
    df_kegg_pathways['entity_a'] = 'GPER1'
    kegg_path_final = df_kegg_pathways.rename(columns={
        'pathway_name': 'entity_b',
        'pathway_id': 'metadata'
    })
    kegg_path_final['relation_type'] = 'in_pathway'
    kegg_path_final['source'] = 'KEGG'
    kegg_path_final = kegg_path_final[['entity_a', 'entity_b', 'relation_type', 'source', 'metadata']]

    # 4. Traitement KEGG GPCRs (Classification)
    print(" Classification des GPCRs (KEGG)...")
    gpcr_class_final = df_kegg_gpcrs.rename(columns={
        'symbol': 'entity_a',
        'family': 'entity_b',
        'class': 'metadata'
    })
    gpcr_class_final['relation_type'] = 'belongs_to_family'
    gpcr_class_final['source'] = 'KEGG'
    gpcr_class_final = gpcr_class_final[['entity_a', 'entity_b', 'relation_type', 'source', 'metadata']]

    # 5. Fusion finale
    print(" Fusion de toutes les relations...")
    df_final = pd.concat([biogrid_final, kegg_path_final, gpcr_class_final], ignore_index=True)

    # 6. Sauvegarde
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    df_final.to_csv(OUTPUT_FILE, index=False)
    
    print(f"Terminé ! {len(df_final)} relations extraites.")
    print(f" Fichier sauvegardé")

if __name__ == "__main__":
    main()
