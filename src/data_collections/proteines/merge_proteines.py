import pandas as pd
import os


script_dir = os.path.dirname(os.path.abspath(__file__))



FILE_ALPHA = os.path.join(script_dir, "../../../data/raw/proteines/AlphaFold_GPER_GPCR.csv")
FILE_UNIPROT = os.path.join(script_dir, "../../../data/raw/proteines/uniprot_GPER_GPCR.csv")

OUTPUT = os.path.join(script_dir ,"../../../data/processed/merged_proteines.csv")



def main():
    print("chargement de donnée ...")

#alpha
    try:
        df_alpha = pd.read_csv(FILE_ALPHA, sep=';')

        df_alpha_clean = pd.DataFrame({
            'source_id' : df_alpha['uniprotAccession'],
            'name' : df_alpha['proteinFullNames'],
            'gene' : df_alpha['gene'],
            'description': df_alpha['uniprotDescription'],
            'organismScientificName' : df_alpha['organismScientificName'],
            'globalMetricValue' : df_alpha['globalMetricValue'],
            'sequenceStart' : df_alpha['sequenceStart'],
            'sequenceEnd' : df_alpha['sequenceEnd'],
            'sequence' : df_alpha['sequence'],
            'Recherche_perso' : df_alpha['Recherche_perso']
        })

    except FileNotFoundError:
        print(f" Erreur : Impossible de trouver {FILE_ALPHA}")
        return

#uniprot
    try:
        df_uniprot = pd.read_csv(FILE_UNIPROT, sep=';')

        df_uniprot_clean = pd.DataFrame({
            'source_id' : df_uniprot['Entry'],
            'chembl' : df_uniprot['ChEMBL'],
            'name' : df_uniprot['Protein names'],
            'gene' : df_uniprot['Gene Names'],
            'description': None,  # Non présent dans uniprot_GPER_GPCR.csv
            'organismScientificName' : df_uniprot['Organism'],
            'globalMetricValue' : None, # Non présent dans uniprot_GPER_GPCR.csv
            'sequence' : df_uniprot['Sequence'],
            'Interacts with' : df_uniprot['Interacts with'],
            'mass' : df_uniprot['Mass'],
            'Binding site' : df_uniprot['Binding site'],
            'Region' : df_uniprot['Region'],
            'Domain [CC]' : df_uniprot['Domain [CC]'],
            'SwissLipids' : df_uniprot['SwissLipids'],
            'Recherche_perso' : df_uniprot['Recherche_perso']

        })

    except FileNotFoundError:
        print(f" Erreur : Impossible de trouver {FILE_UNIPROT}")
        return


    # fusion
    print("concatenation des dataframe ....")
    df_merged = pd.concat([df_alpha_clean, df_uniprot_clean], ignore_index=True)

    # normalisation des alias (ex: GPR30, GPER1, GPER -> GPER)
    alias_map = {
        'GPR30': 'GPER',
        'GPER1': 'GPER',
        'GPER': 'GPER',
        'A8K0M3': 'GPER'
    }

    import re

    def normalize_gene_field(val_name, val_gene):
        # try gene first, then name; return canonical if alias found, else original gene (or name)
        source = ''
        if isinstance(val_gene, str) and val_gene.strip():
            source = val_gene
        elif isinstance(val_name, str) and val_name.strip():
            source = val_name
        else:
            return ''

        # split tokens by common separators
        tokens = re.split(r'[;,/|\s()\[\]-]+', source)
        for t in tokens:
            t_up = t.strip().upper()
            # match ignoring case
            for alias, canonical in alias_map.items():
                if t_up == alias.upper():
                    return canonical

        # no known alias found -> return original gene if present, else first token
        if isinstance(val_gene, str) and val_gene.strip():
            return val_gene
        return tokens[0] if tokens else ''

    df_merged['gene_normalized'] = df_merged.apply(lambda r: normalize_gene_field(r.get('name', ''), r.get('gene', '')), axis=1)

    # Si le source_id est un alias connu (ex: A8K0M3), forcer la normalisation
    try:
        mask = df_merged['source_id'].isin(alias_map.keys())
        if mask.any():
            df_merged.loc[mask, 'gene_normalized'] = df_merged.loc[mask, 'source_id'].map(lambda x: alias_map.get(x, x))
    except Exception:
        pass


    # suppression des doublons
    print ("suppression des doublons ...")
    # conserver un seul enregistrement pour les cas canoniques (ex: GPER)
    # on garde le premier et on supprime les autres quand 'gene_normalized' est identique
    keep_indices = []
    seen = set()
    for idx, gene in df_merged['gene_normalized'].items():
        key = (gene.strip(), df_merged.at[idx, 'organismScientificName'] if 'organismScientificName' in df_merged.columns else '')
        if gene and gene.upper() in [v.upper() for v in alias_map.values()]:
            if key in seen:
                continue
            seen.add(key)
            keep_indices.append(idx)
        else:
            keep_indices.append(idx)

    df_merged = df_merged.loc[keep_indices]

    # assurer également unicité par source_id
    df_merged = df_merged.drop_duplicates(subset=['source_id'], keep='first')

    # en fichier csv 
    df_merged.to_csv(OUTPUT, index=False)
    print(f"Fichier fusionné et dédoublonné sauvegardé : {os.path.abspath(OUTPUT)}")

if __name__ == "__main__":
    main()