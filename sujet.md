# Identification de nouveaux ligands du récepteur GPER dans les glioblastomes via un graphe de connaissances et l’analyse des GPCR

## Contexte

Le **glioblastome (GBM, grade 4 IDH wt)** est la tumeur cérébrale la plus fréquente et agressive chez l’adulte, avec une survie médiane limitée à 15 mois malgré la chirurgie et la radio-chimiothérapie. Des études suggèrent un effet protecteur des œstrogènes, dont l’action passe via les récepteurs nucléaires ERα/ERβ et le récepteur membranaire **GPER**.

Une forte expression de GPER est associée à une meilleure survie, notamment chez les femmes, mais ses ligands endogènes dans le GBM et ses mécanismes de signalisation restent peu connus.

> **Ligands :** Molécules capables de se lier à un récepteur et de moduler son activité (agonistes ou antagonistes). Identifier de nouveaux ligands pour GPER est crucial pour comprendre sa signalisation et développer des approches thérapeutiques ciblées.

## Objectif

L'objectif principal est d'identifier de nouveaux ligands spécifiques de GPER dans le GBM et d'explorer leurs effets. Le projet s'articule autour des axes suivants :

1.  **Classification :** Positionner GPER au sein de la famille des GPCR (Récepteurs couplés aux protéines G).
2.  **Construction de graphe :** Créer un graphe de connaissances intégrant structure, activité et fonction.
3.  **Visualisation :** Développer des outils pour analyser les interactions.
4.  **(OPTIONNEL) Prédiction :** Exploiter le graphe pour prédire de nouveaux ligands via des représentations vectorielles (embeddings) et des algorithmes de prédiction de liens.

## Méthodologie et sources de données

Le projet repose sur l'intégration de données multi-omiques et pharmacologiques :

*   **Séquences et structures :** UniProt, PDB, AlphaFold.
*   **GPCR et ligands connus :** GPCRdb, IUPHAR/BPS, ChEMBL, BindingDB.
*   **Expression génique / protéomique :** TCGA, GTEx, DepMap.
*   **Interactions et voies de signalisation :** STRING, BioGRID, KEGG, Reactome.
*   **Propriétés chimiques et docking :** PubChem, ZINC database.

### Construction du graphe multi-sources

Le graphe sera structuré comme suit :

*   **Nœuds :** Récepteurs, ligands, gènes, voies de signalisation.
*   **Arêtes :** Interactions physiques, corrélations expression-activité, relations structure-fonction.
*   **Analyse :** Utilisation de méthodes de clustering, détection de communautés et identification de nœuds clés (hubs).
*   **Approche prédictive (Optionnelle) :** Utilisation d'embeddings et d'algorithmes de liens pour suggérer de nouvelles interactions ligand-récepteur.

## Résultats attendus

*   Mise en place d'un **graphe de connaissances complet** reliant GPER aux autres GPCR dans le contexte spécifique du GBM.
*   Identification de **ligands candidats** et formulation d'hypothèses sur leurs effets thérapeutiques potentiels.
*   Développement d'un **outil interactif de visualisation** permettant l'exploration et la communication des données.

## Impact

Ce projet pourrait ouvrir la voie à de nouvelles stratégies diagnostiques et thérapeutiques ciblant GPER dans le glioblastome. Il combine bioinformatique, biologie des systèmes et pharmacologie des GPCR pour générer un impact biomédical concret et améliorer la compréhension de cette pathologie complexe.