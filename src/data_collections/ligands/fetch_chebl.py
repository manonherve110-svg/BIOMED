from chembl_webresource_client.new_client import new_client
import pandas as pd


gpcr_proches = ['CHEMBL5872']  # CHEMBL5872 = GPER# interface
activity = new_client.activity # EC50 , KI , IC50  .....
molecule = new_client.molecule # ID , smiles , 

print(f"connection API ... ")



# 1. activité --> pandas consomme tout itérateur
print("1. Récupération des activités...")
import itertools

all_activities = []

for t_id in gpcr_proches:
    print(f"Récupération des activités pour {t_id} ...")
    res = activity.filter(target_chembl_id=t_id)
    df_temp = pd.DataFrame(res)
    df_temp['gpcr_target'] = t_id  # garde trace du GPCR
    all_activities.append(df_temp)

# On concatène tout
df_activity = pd.concat(all_activities, ignore_index=True)


# 2. RÉCUPÉRATION DES STRUCTURES (BATCHING)
unique_ids = df_activity['molecule_chembl_id'].unique().tolist()

smiles_map = {} 


# finction de batching
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))  # [ID1..ID50], [ID51..ID100],.... : on en a 50


# On boucle pour récupérer les données quand c'est possible
for group in chunker(unique_ids, 50):
    try:
        res_mols = molecule.filter(molecule_chembl_id__in=list(group))\
                        .only(['molecule_chembl_id', 'molecule_structures']) # requete API
        
        for mol in res_mols:
            m_id = mol['molecule_chembl_id']
            if mol['molecule_structures'] and 'canonical_smiles' in mol['molecule_structures']:
                smiles_map[m_id] = mol['molecule_structures']['canonical_smiles']
            else:
                smiles_map[m_id] = None # On note  que c'est vide

    except Exception as e:
        print(f"   Erreur sur un paquet (pas grave, on continue) : {e}")





# 3. MAPPING
print("3. Intégration des structures...")
df_activity['canonical_smiles'] = df_activity['molecule_chembl_id'].map(smiles_map)


# 4. NETTOYAGE
print("4. Nettoyage final (On garde tout)...")
colonnes_voulues = [
    'molecule_chembl_id',
    'target_chembl_id', 
    'canonical_smiles', 
    'standard_type', 
    'standard_value', 
    'standard_units', 
    'assay_description'
]

# On garde les colonnes existantes
cols = [c for c in colonnes_voulues if c in df_activity.columns]
df_clean = df_activity[cols]

print(f"--- TERMINÉ : {len(df_clean)} activités récupérées ---")

if not df_clean.empty:   
    output_path = "data/raw/ligands/chembl_gper_ligands.csv"
    df_clean.to_csv(output_path, index=False)
    print(f" SUCCÈS : Tout a été sauvegardé dans {output_path}")
else:
    print(" ÉCHEC : Tableau vide.")