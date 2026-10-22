import re
import math
import os
from typing import Dict, Any, Optional


ATOM_WEIGHTS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Br': 79.904, 'I': 126.90
}


def _count_atoms_smiles(smiles: str) -> Dict[str, int]:
    if not smiles or not isinstance(smiles, str):
        return {}
    s = smiles
    counts = {}
    for two in ("Cl", "Br"):
        n = s.count(two)
        if n:
            counts[two] = n
            s = s.replace(two, '')
    atoms = re.findall(r"[A-Z][a-z]?", s)
    for a in atoms:
        counts[a] = counts.get(a, 0) + 1
    return counts


def estimate_molecular_weight(smiles: str) -> Optional[float]:
    counts = _count_atoms_smiles(smiles)
    if not counts:
        return None
    mw = 0.0
    for atom, n in counts.items():
        mw += ATOM_WEIGHTS.get(atom, 0.0) * n
    return mw


def estimate_logP_simple(smiles: str) -> Optional[float]:
    counts = _count_atoms_smiles(smiles)
    if not counts:
        return None
    hydrophobic = counts.get('C', 0) + counts.get('Cl', 0) + counts.get('Br', 0) + counts.get('I', 0) + counts.get('F', 0)
    polar = counts.get('O', 0) + counts.get('N', 0) + counts.get('S', 0) + counts.get('P', 0)
    diff = hydrophobic - polar
    logp = max(-5.0, min(8.0, diff * 0.6))
    return round(logp, 3)


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def bbb_compatibility_score(smiles: str) -> float:
    counts = _count_atoms_smiles(smiles)
    if not counts:
        return 0.0
    mw = estimate_molecular_weight(smiles) or 0.0
    hbd = counts.get('O', 0) + counts.get('N', 0)
    hba = hbd
    logp = estimate_logP_simple(smiles) or 0.0

    mw_score = 1.0 - max(0.0, (mw - 450.0) / 350.0)
    mw_score = max(0.0, min(1.0, mw_score))

    hbd_score = 1.0 - max(0.0, (hbd - 3.0) / 10.0)
    hbd_score = max(0.0, min(1.0, hbd_score))

    hba_score = 1.0 - max(0.0, (hba - 7.0) / 15.0)
    hba_score = max(0.0, min(1.0, hba_score))

    logp_pref = 1.0 - abs((logp - 2.0) / 4.0)
    logp_pref = max(0.0, min(1.0, logp_pref))

    combined = mw_score * hbd_score * hba_score * logp_pref
    return float(round(_sigmoid((combined - 0.5) * 6.0), 3))


def gbm_compatibility_score(smiles: str) -> float:
    """Retourne un score GBM entre 0 et 100 basé sur MW et logP."""
    if not smiles:
        return 0.0

    # MW et logP
    mw = estimate_molecular_weight(smiles) or 0.0
    logp = estimate_logP_simple(smiles) or 0.0

    # Scores individuels normalisés
    mw_score = max(0.0, min(1.0, 1.0 - (mw - 450.0) / 350.0))
    logp_score = max(0.0, min(1.0, 1.0 - abs(logp - 2.0) / 4.0))

    # Combinaison simple
    combined = mw_score * logp_score

    # Transformation sigmoïde pour étirer la valeur entre 0 et 1
    final_score = 1.0 / (1.0 + math.exp(-6.0 * (combined - 0.5)))

    # Retour en pourcentage
    return round(final_score * 100, 1)


def _detect_smiles_column(columns):
    candidates = ['Smiles', 'smiles', 'SMILES', 'canonical_smiles', 'smiles_raw']
    for c in columns:
        if c in candidates:
            return c
    lowmap = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    return None


if __name__ == '__main__':
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=None, help='Path to ligands_GPCR.csv')
    parser.add_argument('--nrows', type=int, default=None, help='Limit rows (for quick tests)')
    parser.add_argument('--output', default=None, help='Output CSV path')
    args = parser.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    input_path = args.input or os.path.join(repo_root, 'data', 'processed', 'ligands_GPCR.csv')
    output_path = args.output or os.path.join(repo_root, 'data', 'processed', 'ligands_GPCR_scored.csv')

    if not os.path.exists(input_path):
        print('Input not found:', input_path)
        raise SystemExit(1)

    df = pd.read_csv(input_path, dtype=str, nrows=args.nrows)
    df = df.fillna('')
    smiles_col = _detect_smiles_column(list(df.columns))
    if not smiles_col:
        print('No SMILES column detected in', input_path)
        print('Columns:', df.columns.tolist())
        raise SystemExit(1)

    df['gbm_score_percent'] = df[smiles_col].apply(gbm_compatibility_score)
    out = pd.concat([df.reset_index(drop=True), est_df.reset_index(drop=True)], axis=1)
    out.to_csv(output_path, index=False)
    print('Scored', len(out), 'rows ->', output_path)
