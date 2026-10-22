import pandas as pd
import requests
import os

# --- CONFIGURATION ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../../data/raw/expression/gtex_brain_gper_expression.csv")

TARGET_SYMBOL = "GPER1"

# API Officielle GTEx Portal v2
API_BASE = "https://gtexportal.org/api/v2"

def get_gencode_id(symbol):
    """ Récupère l'ID Ensembl """
    print(f"🔍 Recherche de l'ID Ensembl pour {symbol}...")
    url = f"{API_BASE}/reference/gene"
    params = {"geneId": symbol, "pageSize": 1}
    try:
        r = requests.get(url, params=params)
        if r.status_code == 200:
            data = r.json()
            if data and 'data' in data and len(data['data']) > 0:
                gencode_id = data['data'][0]['gencodeId']
                print(f"   ✅ ID trouvé : {gencode_id}")
                return gencode_id
    except Exception as e:
        print(f"❌ Erreur recherche ID : {e}")
    return None

def fetch_expression_data(gencode_id):
    """ Récupère les données groupées par tissu """
    print("📡 Récupération des données d'expression (API GTEx)...")
    url = f"{API_BASE}/expression/geneExpression"
    params = {
        "gencodeId": gencode_id,
        "datasetId": "gtex_v8",
        "pageSize": 10000 
    }
    try:
        r = requests.get(url, params=params)
        if r.status_code == 200:
            data = r.json()
            return data.get('data', [])
    except Exception as e:
        print(f"❌ Erreur connexion : {e}")
    return []

def main():
    # 1. ID Gène
    gencode_id = get_gencode_id(TARGET_SYMBOL)
    if not gencode_id: return

    # 2. Données
    raw_data = fetch_expression_data(gencode_id)
    if not raw_data:
        print("⚠️ Aucune donnée reçue.")
        return

    # 3. DataFrame
    df = pd.DataFrame(raw_data)
    
    # 4. Filtrer pour le CERVEAU
    print("🧠 Filtrage pour les tissus cérébraux...")
    # On vérifie d'abord si la colonne existe
    if 'tissueSiteDetailId' not in df.columns:
        print("❌ Colonne 'tissueSiteDetailId' absente. Colonnes dispos :", df.columns.tolist())
        return

    df_brain = df[df['tissueSiteDetailId'].str.contains("Brain", case=False, na=False)].copy()
    
    if df_brain.empty:
        print("❌ Aucun échantillon de cerveau trouvé.")
        return

    print(f"   ℹ️ {len(df_brain)} régions cérébrales trouvées (données agrégées).")

    # 5. "EXPLOSION" DES DONNÉES (La correction est ici)
    # La colonne 'data' contient des listes [val1, val2...]. On les transforme en lignes.
    print("💥 Dépliage des listes de valeurs (Explode)...")
    df_exploded = df_brain.explode('data')

    # 6. Nettoyage
    # On garde la valeur ('data') et le tissu ('tissueSiteDetailId')
    df_final = df_exploded[['data', 'tissueSiteDetailId']].copy()
    
    df_final = df_final.rename(columns={
        'data': 'Expression_TPM',
        'tissueSiteDetailId': 'Brain_Region'
    })
    
    # On ajoute une colonne Type
    df_final['Type'] = 'Healthy (GTEx)'
    
    # On crée des faux SampleID car GTEx ne les donne pas ici (c'est juste pour avoir une structure compatible)
    df_final['SampleID'] = [f"GTEX_Sim_{i}" for i in range(len(df_final))]

    # 7. Sauvegarde
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df_final.to_csv(OUTPUT_FILE, index=False)
    
    print(f"✅ SUCCÈS ! {len(df_final)} échantillons récupérés.")
    print(f"📂 Fichier : {os.path.abspath(OUTPUT_FILE)}")
    
    avg_expr = df_final['Expression_TPM'].mean()
    print(f"\n--- Diagnostic ---")
    print(f"   Moyenne Expression (TPM) : {avg_expr:.4f}")

if __name__ == "__main__":
    main()