import requests 
import pandas as pd
import os


KEY_ACCES = "...."
GENE_SYMBOL = "GPER1"
OUTPUT_FILE = "../../../data/raw/voie_signalisation/biogrid_gper_interaction.csv"


def fetch_pathway(gene_symbole) :
    """ recuperation donnée voie_signalisation de biogrid  --> retourne dataframe des pathway"""

    url = "https://webservice.thebiogrid.org/interactions/"

    params = {
        "accesskey": KEY_ACCES,
        "searchNames": "true",  # Chercher par nom de gène
        "geneList": gene_symbole,
        "taxId": 9606,          # Humain
        "format": "json",       # Format de réponse
        "includeInteractors": "true" # Inclure les partenaires
    }

    print("conexion api ....")


    try:
        reponse = requests.get(url, params=params)
        if reponse.status_code == 200:
            data = reponse.json()

            pathway = list(data.values())
            return pd.DataFrame(pathway)
        else:
            print(f"Erreur BioGRID : {reponse.status_code} (Vérifiez votre clé API)")
            return pd.DataFrame()

    except Exception as e:
        print(f" Erreur : {e}")
        return pd.DataFrame()


def main():
    if KEY_ACCES == "....." or KEY_ACCES == "VOTRE_CLE_API_ICI":
        print("Veuillez renseigner votre clé API BioGRID.")
        return


    df = fetch_pathway(GENE_SYMBOL)

    if not df.empty:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.normpath(os.path.join(script_dir, OUTPUT_FILE))

        # Création du dossier si il n'existe pas
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        cols = ['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B', 'EXPERIMENTAL_SYSTEM', 'PUBMED_ID']

        df_clean = df[cols] if set(cols).issubset(df.columns) else df
        
        df_clean.to_csv(save_path, index=False)
        print(f" {len(df)} interactions BioGRID trouvées.")
        print(f" Sauvegardé : {save_path}")
    else:
        print("Rien trouvé.")

if __name__ == "__main__":
    main()
