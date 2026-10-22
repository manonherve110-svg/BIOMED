import pandas as pd
import requests
import time
import os

# --- CHEMINS ---
script_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(script_dir, "../../data/raw/ligands/pubchem_gper_ligands.csv")
OUTPUT_FILE = os.path.join(script_dir, "../../data/raw/ligands/pubchem_gper_ligands_smiles.csv")

BATCH_SIZE = 50

def recup_smile(list_cid):
    """ Récupère les SMILES via l'API PubChem """
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
    cids_str = ",".join(list_cid) #
    url = f"{base_url}/{cids_str}/property/IsomericSMILES/JSON"
    
    try:
        reponse = requests.get(url)
        if reponse.status_code == 200:
            data = reponse.json()
            if 'PropertyTable' in data:
                properties = data['PropertyTable']['Properties']

                return {str(item['CID']): item.get('SMILES', '') for item in properties}
        return {}
    except Exception as e:
        print(f" Erreur connexion : {e}")
        return {}

def main():
    print(f" Lecture de : {os.path.basename(INPUT_FILE)}")
    try:
        df = pd.read_csv(INPUT_FILE)
    except FileNotFoundError:
        print(" Fichier introuvable.")
        return


    print(" Préparation des identifiants (Mode Texte)...")
    
    # On convertit d'abord en numérique pour gérer les erreurs, puis on vire les NaNs
    temp_series = pd.to_numeric(df['Compound_CID'], errors='coerce').dropna()

    unique_cids_str = temp_series.astype(int).astype(str).unique().tolist()
    
    print(f" {len(unique_cids_str)} molécules uniques à chercher.")

    cid_to_smiles = {}
    

    
    # --- ETAPE 2 : RECUPERATION ---
    print(" Démarrage de la récupération...")
    for i in range(0, len(unique_cids_str), BATCH_SIZE):
        batch = unique_cids_str[i:i + BATCH_SIZE]
        results = recup_smile(batch)
        cid_to_smiles.update(results)
        time.sleep(0.2)

    df['cid_join_key'] = pd.to_numeric(df['Compound_CID'], errors='coerce')
    df['cid_join_key'] = df['cid_join_key'].fillna(-1).astype(int).astype(str)
    df.loc[df['cid_join_key'] == "-1", 'cid_join_key'] = None




    # --- ETAPE 4 : MAPPING ---
    print(" Fusion des données...")
    df['SMILES'] = df['cid_join_key'].map(cid_to_smiles)
    

    # Nettoyage
    df.drop(columns=['cid_join_key'], inplace=True)

    nb_smiles = df['SMILES'].notna().sum()
    print(f" Résultat final : {nb_smiles} SMILES trouvés.")

    df.to_csv(OUTPUT_FILE, index=False)
    print(f" Fichier sauvegardé : {OUTPUT_FILE}")

if __name__ == "__main__":
    main()