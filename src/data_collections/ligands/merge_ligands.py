import pandas as pd
from rdkit import Chem
import os


script_dir = os.path.dirname(os.path.abspath(__file__))

FILE_PUBCHEM = os.path.join(script_dir, "../../../data/raw/ligands/pubchem_gper_ligands.csv")
FILE_CHEMBL  = os.path.join(script_dir, "../../../data/raw/ligands/chembl_gper_ligands.csv")

OUTPUT_FILE  = os.path.join(script_dir, "../../../data/processed/merged_ligands.csv")


def get_canonical_smiles(smiles_str):
    """ 
    Nettoie et standardise le SMILES grâce à RDKit. 
    Permet de comparer les molécules entre elles.
    """
    if pd.isna(smiles_str) or str(smiles_str).strip() == "":
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles_str))
        if mol:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except:
        pass
    return None


def main():
    print(" Chargement des données...")

    # PUBCHEM
    try:
        df_pub = pd.read_csv(FILE_PUBCHEM)
        
        df_pub_clean = pd.DataFrame({
            'source_id': df_pub['Compound_CID'],
            'Protein_Accession': df_pub['Protein_Accession'],
            'origin': 'PubChem',
            'activity_type': df_pub['Activity_Type'],
            'activity_value': df_pub['Activity_Value'],
            'activity_unit': 'nM', 
            'smiles_raw': df_pub['SMILES']
        })
    except FileNotFoundError:
        print(f" Erreur : Impossible de trouver {FILE_PUBCHEM}")
        return



    # CHEMBL
    try:
        df_chem = pd.read_csv(FILE_CHEMBL) 

        df_chem_clean = pd.DataFrame({
            'source_id': df_chem['molecule_chembl_id'],
            'origin': 'ChEMBL',
            'activity_type': df_chem['standard_type'],
            'activity_value': df_chem['standard_value'],
            'activity_unit': df_chem['standard_units'],
            'smiles_raw': df_chem['canonical_smiles']
        })
    except FileNotFoundError:
        print(f" Erreur : Impossible de trouver {FILE_CHEMBL}")
        return



    #  FUSION 
    print(" Fusion des deux datasets...")
    df_merged = pd.concat([df_pub_clean, df_chem_clean], ignore_index=True)

    print("  Standardisation des structures chimiques avec RDKit (patience)...")
    df_merged['canonical_smiles'] = df_merged['smiles_raw'].apply(get_canonical_smiles)

    # Nettoyage des SMILES invalides
    before_drop = len(df_merged)
    df_merged = df_merged.dropna(subset=['canonical_smiles'])
    after_drop = len(df_merged)
    print(f"   - Molécules valides : {after_drop} (suppression de {before_drop - after_drop} SMILES invalides)")

    # DEDUPLICATION
    print("  Suppression des doublons...")
    df_merged = df_merged.sort_values(by=['canonical_smiles', 'origin'])
    
    before_dedup = len(df_merged)
    df_merged = df_merged.drop_duplicates(subset=['canonical_smiles'], keep='first')
    after_dedup = len(df_merged)
    print(f"   - Molécules uniques : {after_dedup} (suppression de {before_dedup - after_dedup} doublons)")

    # SAUVEGARDE 
    df_merged.to_csv(OUTPUT_FILE, index=False)
    print(f"Fichier fusionné et dédoublonné sauvegardé : {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    main()