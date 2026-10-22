# BIOMED: Identification de Ligands GPER dans les Glioblastomes

### Contexte Scientifique

Le glioblastome (GBM, grade 4) est la tumeur cérébrale la plus agressive chez l'adulte, avec une survie médiane de 15 mois. Des études récentes suggèrent un rôle protecteur des œstrogènes via le récepteur membranaire GPER. Cette plateforme intègre des données multi-omiques pour identifier des ligands candidats GPER potentiellement thérapeutiques.

---

## Quick Start

### Prérequis
- Python 3.8+
- pip/conda

### Installation

```bash
# Sur macOS/Linux
python3 -m venv venv
source venv/bin/activate

# Sur Windows
python -m venv venv
venv\Scripts\activate

# Installer les dépendances sur Windows
pip install -r requirements.txt

# Installer les dépendances sur Mac
pip install -r requirements-mac.txt
```

### Lancer l'application web

```bash
python3 -m web.app
```

L'application sera accessible sur `http://127.0.0.1:5000`

---

## Fonctionnalités Principales

### 1. **Collecte de Données Multi-Sources**
- **Ligands** : ChEMBL (validés) + PubChem (screening)
- **Protéines** : UniProt + structures AlphaFold
- **Expression** : TCGA GBM + GTEx Brain
- **Interactions** : BioGRID + voies KEGG
- **Voies de Signalisation** : KEGG G-couplées

### 2. **Analyses Computationnelles**
- **Similarité structurale Jaccard** : Identification GPCR proches de GPER
- **Blocs conservés dynamiques** : Détection de motifs fonctionnels (DRY, CWxP)
- **Scoring GBM** : Évaluation compatibilité ligands (MW, LogP, BBB)
- **Construction graphe** : Réseau multi-relation typed

### 3. **Interface Web Interactive**
- **Visualisation réseau** : Graphe interactif via vis.js
- **Exploration données** : Consultation CSV dynamique
- **Filtrage avancé** : Types de relations, profondeur BFS
- **Export** : Données et graphe

---

##  Structure du Projet

```
biomed/
├── src/                           # Code source métier
│   ├── data_collections/         # Scripts collecte données
│   │   ├── expression/           # TCGA, GTEx
│   │   ├── ligands/              # ChEMBL, PubChem
│   │   ├── proteines/            # UniProt, AlphaFold
│   │   └── voie_signalisation/   # BioGRID, KEGG
│   ├── graphe/                   # Analyses graphe
│   │   ├── similarite_jaccord_proteine.py    # Jaccard k-mer
│   │   └── plus_longue_sous_chaine.py        # Blocs conservés
│   ├── ligands/
│   │   └── ligand_predictor.py   # Score GBM
│   └── visualisation/
│       └── gper.ipynb            # Notebook exploration
├── web/                          # Application Flask
│   ├── app.py                    # Serveur principal
│   ├── graph_service.py          # Logique graphe
│   ├── graph_modes.py            # Modes visualisation
│   ├── data_loader.py            # Chargement données
│   ├── static/                   # CSS, JS, images
│   └── templates/                # HTML (menu, graphe, csv)
├── data/
│   ├── raw/                      # Données sources (APIs, téléchargements)
│   └── processed/                # Données nettoyées/fusionnées
|
├── requirements.txt              # Dépendances Python
└── README.md                     # Ce fichier
```

---

##  Technologies

| Composant | Technologie | Raison |
|-----------|-------------|--------|
| **Backend** | Flask + Python 3 | Écosystème scientifique riche |
| **Graphe** | NetworkX | Analyse algorithmes (clustering, BFS) |
| **Frontend** | HTML5 + vis.js + JS | Visualisation réseau performante |
| **Données** | pandas + CSV | Reproductibilité, versionnage Git |
| **Environnement** | venv + pip | Isolation dépendances |

---


## Utilisation

### Accéder à l'Application Web

1. **Menu** : Navigation et informations projet
2. **Graphe** : 
   - Sélectionner nœud graine (GPER, GPCR ou ligand)
   - Filtrer par type de relation
   - Ajuster profondeur BFS et limite nœuds
   - Cliquer pour détails
3. **CSV** : Consulter données brutes
   - Filtrer par catégorie (ligands, protéines, expression)
   - Export possible

---

##  Algorithmes de Scoring

### Score GBM Compatibility

Évalue la probabilité qu'un composé soit un ligand actif GPER en contexte GBM :

$$S_{\text{GBM}} = \sigma(6 \times (S_{MW} \times S_{\text{LogP}} - 0.5)) \times 100$$

Où :
- $S_{MW}$ = score poids moléculaire (optimum 300-500 Da, BBB)
- $S_{\text{LogP}}$ = score lipophilie (optimum 1-5, traversabilité)
- $\sigma$ = fonction sigmoïde (compression non-linéaire)

**Interprétation** :
- **80-100%** : Ligand très compatible → priorité 1
- **60-79%** : Compatible → validation recommandée
- **<60%** : Peu compatible → BBB challengé

Voir le rapport sur arche
### Similarité Structurale Jaccard

Identification GPCR proches de GPER via k-mers (k=3) :

$$\text{Jaccard} = \frac{\sum \min(c_1(k_i), c_2(k_i))}{\sum \max(c_1(k_i), c_2(k_i))}$$

Seuil : similarité > 0.15 (conservation minimale 7-TM)

---

##  Logique d'Intégration

```
Données Brutes (APIs/fichiers)
  ↓
[src/data_collections/] Nettoyage & fusion
  ↓
[data/processed/] CSV normalisés
  ↓
[web/graph_service.py] Construction graphe
  ↓
[web/app.py] API REST
  ↓
[templates/] Interface Web Interactive
```
---

##  Limitations et Perspectives

### Limitations Actuelles
- Score MW/LogP simplifié (vs XLOGP3 propriétaire)
- Pas de modélisation 3D/docking (phase future)
- Expression TCGA sur 153 samples (future : stratification sous-types)
- Similarité Jaccard insensible à fonction protéique fine

### Développements Futurs (Court-Moyen Terme)
-  **Node2vec embeddings** : Similarité latente graphe
-  **Link prediction** : PageRank pour découvrir interactions
-  **Docking rapide** : Smina scoring vs AlphaFold GPER
-  **Single-cell RNA-seq** : Hétérogénéité intra-tumorale
-  **PAINS filter** : Triage promiscuous binders

---

##  Contribution

Ce projet a été développé dans le contexte d'une formation en bioinformatique par Ambre EVANO, Manon Hervé et Ingrid Gompel
