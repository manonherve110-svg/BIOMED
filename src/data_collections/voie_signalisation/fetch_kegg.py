import requests
import pandas as pd
import re
import os


BRITE_ID = "br:hsa04030"
OUTPUT_FILE = "../../../data/raw/voie_signalisation/kegg_all_gpcrs.csv"

def fetch_gpcr_family():
    """
    Récupère la liste hiérarchique des GPCRs humains depuis KEGG BRITE.
    """
    url = f"http://rest.kegg.jp/get/{BRITE_ID}"
    print(f" Téléchargement de la classification GPCR ({BRITE_ID})...")
    
    try:
        response = requests.get(url)
        if response.status_code != 200:
            print(f" Erreur API : {response.status_code}")
            return []
        
        # Le fichier est un texte hiérarchique (A, B, C, D...)
        # A = Super Famille
        # B = Classe (ex: Class A: Rhodopsin-like)
        # C = Sous-famille
        # D = Le Gène (Ce qu'on veut !)
        
        gpcrs = []
        current_class = "Unknown"
        current_family = "Unknown"
        
        for line in response.text.split('\n'):

            # Détection de la Classe (B)
            if line.startswith('B  '):
                current_class = line[3:].strip()

            # Détection de la Famille (C)
            elif line.startswith('C    '):
                current_family = line[5:].strip()

                
            # Détection du Gène (D) -> C'est la ligne de donnée
            elif line.startswith('D      '):
                # Format typique : "D      2852  GPER1; G protein-coupled estrogen receptor 1"
                content = line[7:].strip()
                
                # On utilise une expression régulière (Regex) pour séparer l'ID du reste
                # On cherche des chiffres au début
                match = re.match(r'^(\d+)\s+(.*)', content)
                if match:
                    gene_id = match.group(1) # ex: 2852
                    description = match.group(2) # ex: GPER1; G protein...
                    
                    # Séparation du symbole (GPER1) et du nom complet
                    if ";" in description:
                        symbol, full_name = description.split(';', 1)
                    else:
                        symbol = description
                        full_name = description
                    
                    gpcrs.append({
                        'gene_id': f"hsa:{gene_id}", # On remet le format hsa: pour être compatible
                        'symbol': symbol.strip(),
                        'name': full_name.strip(),
                        'class': current_class,
                        'family': current_family
                    })
        
        return gpcrs

    except Exception as e:
        print(f" Erreur de script : {e}")
        return []

def main():
    gpcr_list = fetch_gpcr_family()
    
    if gpcr_list:
        df = pd.DataFrame(gpcr_list)
        
        # Sauvegarde
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(script_dir, OUTPUT_FILE)
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        df.to_csv(save_path, index=False)
        print(f"{len(df)} GPCRs récupérés.")
        print(f"Sauvegardé ici : {save_path}")
        
    else:
        print("Aucun GPCR trouvé.")

if __name__ == "__main__":
    main()

