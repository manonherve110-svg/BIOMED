import pandas as pd
import os
import io
import requests

# --- CONFIGURATION ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../../data/raw/expression/tcga_gbm_gper_clinical.csv")

# URLs UCSC Xena
URL_EXPRESSION = "https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HiSeqV2.gz"
URL_CLINICAL = "https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/GBM_clinicalMatrix"

TARGET_GENES = ["GPER1", "GPER", "GPR30", "ESR1", "ESR2"] 

def fetch_xena_data(url):
    print(f"Téléchargement ...")
    try:
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            comp = 'gzip' if url.endswith('.gz') else 'infer'
            return pd.read_csv(io.BytesIO(r.content), sep='\t', compression=comp, index_col=0)
    except Exception as e:
        print(f"Erreur : {e}")
    return None

def main():
    # 1. EXPRESSION
    df_exp = fetch_xena_data(URL_EXPRESSION)
    if df_exp is None: return

    # Filtre des gènes
    genes_found = [g for g in TARGET_GENES if g in df_exp.index]
    if not genes_found:
        print(" Gènes non trouvés.")
        return
    
    # Transposition (Lignes = Patients)
    df_gper = df_exp.loc[genes_found].T
    
    print(" Harmonisation des identifiants (Sample -> Patient)...")
    df_gper['patient_id'] = df_gper.index.map(lambda x: x[:-3] if len(x) > 12 else x)

    # 2. CLINIQUE
    df_clin = fetch_xena_data(URL_CLINICAL)
    if df_clin is None: return

    # 3. FUSION INTELLIGENTE
    print(" Fusion avec les données de survie...")
    df_final = pd.merge(df_gper, df_clin, left_on='patient_id', right_index=True, how='inner')

    # 4. NETTOYAGE FINAL

    cols_map = {
        'OS.time': 'Survival_Days', 
        'OS': 'Survival_Status',
        'age_at_initial_pathologic_diagnosis': 'Age',
        'gender': 'Gender'
    }
    df_final = df_final.rename(columns=cols_map)

    # On sélectionne les colonnes finales
    wanted_cols = genes_found + ['Survival_Days', 'Survival_Status', 'Age', 'Gender', 'sample_type']
    available_cols = [c for c in wanted_cols if c in df_final.columns]
    
    df_final = df_final[available_cols]

    # Sauvegarde
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df_final.to_csv(OUTPUT_FILE)
    
    print(f" Terminé ! ")



if __name__ == "__main__":
    main()