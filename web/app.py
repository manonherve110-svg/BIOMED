from __future__ import annotations

from flask import Flask, jsonify, render_template, request

from .data_loader import load_all_data
from .graph_modes import get_graph_modes
from .graph_service import build_subgraph, resolve_entity


def create_app() -> Flask:
    app = Flask(__name__)
    data = load_all_data()

    # Pages
    @app.route("/")
    def menu():
        return render_template("menu.html")

    @app.route("/csv")
    def csv():
        return render_template("csv.html")

    @app.route("/graphe")
    def graphe():
        return render_template("graphe.html")

    # API: modes
    @app.route("/api/graph_modes")
    def graph_modes():
        return jsonify(get_graph_modes())

    # API: subgraph
    @app.route("/api/subgraph")
    def api_subgraph():
        seed = request.args.get("seed", "GPER").strip()
        k = int(request.args.get("k", 1))
        rels = [r.strip() for r in request.args.get("rels", "INTERACTS_WITH").split(",") if r.strip()]
        max_n = int(request.args.get("max_n", 200))

        payload = build_subgraph(data=data, seed=seed, k=k, rels=rels, max_n=max_n)

        # si build_subgraph renvoie une erreur structurée
        if "error" in payload:
            return jsonify(payload), 400

        return jsonify(payload)

    # API: entity
    @app.route("/api/entity/<node_id>")
    def api_entity(node_id: str):
        payload, status = resolve_entity(data, node_id)
        return jsonify(payload), status

    # API: données CSV dynamique
    @app.route("/api/data")
    def api_data():
        category = request.args.get("category", "ligands").lower()
        source = request.args.get("source", "all").lower()
        
        # Mapping catégorie -> données
        if category == "ligands":
            if source == "pubchem":
                df = data["ligands"][data["ligands"]["origin"] == "PubChem"]
            elif source == "chembl":
                df = data["ligands"][data["ligands"]["origin"] == "ChEMBL"]
            else:
                df = data["ligands"]
        elif category == "gper":
            if source == "uniprot":
                df = data["uniprot"][data["uniprot"]["Recherche_perso"] == "GPER"]
            elif source == "alphafold":
                df = data["alphafold"][data["alphafold"]["Recherche_perso"] == "GPER"]
            else:
                df = data["proteines"][data["proteines"]["Recherche_perso"] == "GPER"]
        elif category == "gpcr":
            if source == "uniprot":
                df = data["uniprot"][data["uniprot"]["Recherche_perso"] == "GPCR"]
            elif source == "alphafold":
                df = data["alphafold"][data["alphafold"]["Recherche_perso"] == "GPCR"]
            else:
                df = data["proteines"][data["proteines"]["Recherche_perso"] == "GPCR"]
        elif category == "tcga":
            df = data["tcga"]
        else:
            df = data.get(category, data["ligands"])
        
        return jsonify(df.fillna("").to_dict(orient="records"))

    # API: sources disponibles par catégorie
    @app.route("/api/sources")
    def api_sources():
        return jsonify({
            "ligands": ["Tous", "PubChem", "ChEMBL"],
            "gper": ["Tous", "UniProt", "AlphaFold"],
            "gpcr": ["Tous", "UniProt", "AlphaFold"],
            "tcga": ["TCGA"]
        })

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
