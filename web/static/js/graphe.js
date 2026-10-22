
(function () {
  function setStatus(msg) {
    const el = document.getElementById("statusLine");
    if (el) el.textContent = msg || "";
  }

  async function loadGraph(params) {
    const panel = document.getElementById("proteinPanel");
    try {
      setStatus("Chargement…");
      const qs = new URLSearchParams(params).toString();
      const res = await fetch(`/api/subgraph?${qs}`);
      const text = await res.text();

      if (!res.ok) {
        panel.innerHTML = `<pre style="white-space:pre-wrap;color:#b91c1c;">API ERROR ${res.status}\n\n${text}</pre>`;
        setStatus("Erreur API");
        return;
      }

      const payload = JSON.parse(text);

      window.__KG__.nodes.clear();
      window.__KG__.edges.clear();
      window.__KG__.nodes.add(payload.nodes || []);
      window.__KG__.edges.add(payload.edges || []);

      setStatus(`OK — ${window.__KG__.nodes.length} nœuds, ${window.__KG__.edges.length} arêtes`);
      window.__KG__.network.fit({ animation: true });
    } catch (e) {
      panel.innerHTML = `<pre style="white-space:pre-wrap;color:#b91c1c;">JS ERROR\n\n${e}</pre>`;
      setStatus("Erreur JS");
    }
  }

  async function loadEntity(nodeId) {
    const panel = document.getElementById("proteinPanel");

    const res = await fetch(`/api/entity/${encodeURIComponent(nodeId)}`);
    const text = await res.text();

    if (!res.ok) {
      panel.innerHTML = `
        <div class="text-danger"><strong>Aucune fiche trouvée</strong> pour : ${nodeId}</div>
        <pre style="white-space:pre-wrap;margin-top:8px;">HTTP ${res.status}\n${text}</pre>
      `;
      return;
    }

    const feat = JSON.parse(text);
    const rows = Object.entries(feat).map(([k, v]) =>
      `<tr><th style="width:45%">${k}</th><td>${v ?? ""}</td></tr>`
    ).join("");

    panel.innerHTML = `
      <div class="mb-2"><strong>ID :</strong> ${nodeId}</div>
      <table class="table table-sm"><tbody>${rows}</tbody></table>
    `;
  }

  async function loadGraphModes() {
    const box = document.getElementById("graphModes");
    box.innerHTML = `<div class="text-muted">Chargement…</div>`;

    const res = await fetch("/api/graph_modes");
    if (!res.ok) {
      box.innerHTML = `<div class="text-danger">Impossible de charger /api/graph_modes</div>`;
      return;
    }

    const modes = await res.json();
    box.innerHTML = "";

    const maxNInput = document.getElementById("maxNInput");
    for (const m of modes) {
      const btn = document.createElement("button");
      btn.className = "btn btn-outline-dark text-start";
      btn.innerHTML = `
        <div style="font-weight:600">${m.title}</div>
        <div class="small text-muted">${m.desc || ""}</div>
      `;
      btn.addEventListener("click", () => {
        const params = Object.assign({}, m.params || {});
        if (maxNInput) params.max_n = parseInt(maxNInput.value || 50, 10);
        loadGraph(params);
      });
      box.appendChild(btn);
    }
  }

  function initNetwork() {
    const container = document.getElementById("mynetwork");
    const nodes = new vis.DataSet([]);
    const edges = new vis.DataSet([]);

    const options = {
      edges: { smooth: { enabled: true, type: "dynamic" } },
      interaction: { dragNodes: true },
      physics: {
        barnesHut: {
          centralGravity: 0.3,
          damping: 0.09,
          gravitationalConstant: -80000,
          springConstant: 0.001,
          springLength: 250
        },
        enabled: true,
        stabilization: { enabled: true, iterations: 800, updateInterval: 50 }
      }
    };

    const network = new vis.Network(container, { nodes, edges }, options);
    window.__KG__ = { network, nodes, edges, loadGraph };

    network.on("click", function (params) {
      const panel = document.getElementById("proteinPanel");

      // 1. SI ON CLIQUE SUR UN NŒUD
      if (params.nodes.length > 0) {
        const nodeId = params.nodes[0];

        // --- Intercepter le Glioblastome ---
        if (nodeId === "GBM") {
          const nodeData = window.__KG__.nodes.get(nodeId);
          if (nodeData && nodeData.clinical_info) {
            const info = nodeData.clinical_info;
            panel.innerHTML = `
                      <h5 class="mb-3" style="color: #e74c3c;">Glioblastoma Multiforme</h5>
                      <ul class="list-group list-group-flush">
                          <li class="list-group-item"><b>Base de données :</b> TCGA</li>
                          <li class="list-group-item"><b>Taille de la cohorte :</b> ${info.cohort_size} patients</li>
                          <li class="list-group-item"><b>Âge moyen :</b> ${info.mean_age} ans</li>
                          <li class="list-group-item"><b>Répartition :</b> ${info.male_count} Hommes / ${info.female_count} Femmes</li>
                      </ul>
                      <div class="mt-3 text-muted small">
                          Ces données cliniques valident l'ancrage tumoral des protéines affichées sur le graphe.
                      </div>
                  `;
          }
          return;
        }

        // --- Intercepter le Tissu Sain ---
        if (nodeId === "HEALTHY_BRAIN") {
          panel.innerHTML = `
                  <h5 class="mb-3" style="color: #2ecc71;">Cerveau Sain</h5>
                  <ul class="list-group list-group-flush">
                      <li class="list-group-item"><b>Base de données :</b> GTEx Portal</li>
                      <li class="list-group-item"><b>Régions :</b> Substantia nigra, Cortex, etc.</li>
                  </ul>
              `;
          return;
        }

        // --- Pour les Protéines et Ligands (Appel de ton API existante) ---
        loadEntity(nodeId);
        return;
      }

      // 2. SI ON CLIQUE SUR UNE ARÊTE
      if (params.edges.length > 0 && params.nodes.length === 0) {
        const edgeId = params.edges[0];
        const clickedEdge = edges.get(edgeId);

        if (clickedEdge) {
          let htmlContent = `<h5 class="mb-3">Détails de l'interaction</h5>`;
          htmlContent += `<ul class="list-group list-group-flush">`;
          htmlContent += `<li class="list-group-item"><b>Source :</b> ${clickedEdge.from}</li>`;
          htmlContent += `<li class="list-group-item"><b>Cible :</b> ${clickedEdge.to}</li>`;
          htmlContent += `<li class="list-group-item"><b>Type :</b> ${clickedEdge.label || "Inconnu"}</li>`;

          if (clickedEdge.score !== undefined && clickedEdge.score !== null) {
            htmlContent += `<li class="list-group-item">
                      <b>Score Jaccard :</b> 
                      <span class="badge bg-success fs-6">${clickedEdge.score}</span>
                  </li>`;
          }
          if (clickedEdge.title) {
            htmlContent += `<li class="list-group-item text-muted"><small>${clickedEdge.title}</small></li>`;
          }
          htmlContent += `</ul>`;
          panel.innerHTML = htmlContent;
        }
      } else {
        panel.innerHTML = '<div class="text-muted">Clique sur un nœud ou une arête pour afficher ses infos.</div>';
      }
    });
  }


  function bindUI() {
    document.getElementById("btnBuild").addEventListener("click", () => {
      const seed = (document.getElementById("seedInput").value || "").trim() || "GPER";
      const k = parseInt(document.getElementById("kInput").value || "1", 10);
      const max_n = parseInt(document.getElementById("maxNInput").value || "50", 10);

      const rels = [];
      const predEl = document.getElementById("relPRED");
      if (document.getElementById("relPPI").checked) rels.push("INTERACTS_WITH");
      if (document.getElementById("relBINDS").checked) rels.push("BINDS_TO");
      if (predEl && predEl.checked) rels.push("PREDICTED_BINDS");

      // Ajout des nouvelles cases :
      if (document.getElementById("relTCGA") && document.getElementById("relTCGA").checked) {
        rels.push("EXPRESSED_IN_GBM");       // résumé (noeud GBM)
      }
      const tcgaDetailEl = document.getElementById("relTCGA_DETAIL");
      if (tcgaDetailEl && tcgaDetailEl.checked) {
        rels.push("GBM_PATIENTS");           // patients individuels reliés à GPER
      }
      if (document.getElementById("relGTEX") && document.getElementById("relGTEX").checked) {
        rels.push("EXPRESSED_IN_BRAIN");
      }

      loadGraph({ seed, k, rels: rels.join(","), max_n });
    });

    document.getElementById("btnFit").addEventListener("click", () => {
      window.__KG__.network.fit({ animation: true });
    });
  }

  // Boot
  document.addEventListener("DOMContentLoaded", () => {
    // vis-network est chargé via defer; DOMContentLoaded suffit en général.
    // Si jamais vis n'est pas défini (rare), c’est que l’ordre de scripts a été modifié.
    if (typeof vis === "undefined") {
      const panel = document.getElementById("proteinPanel");
      panel.innerHTML = `<pre style="white-space:pre-wrap;color:#b91c1c;">Erreur: vis-network non chargé.</pre>`;
      return;
    }

    initNetwork();
    bindUI();
    loadGraphModes();
    // initial load respects the default maxNInput value
    const initialMax = parseInt((document.getElementById("maxNInput") || { value: 50 }).value, 10) || 50;
    loadGraph({ seed: "GPER", k: 1, rels: "INTERACTS_WITH", max_n: initialMax });
  });
})();
