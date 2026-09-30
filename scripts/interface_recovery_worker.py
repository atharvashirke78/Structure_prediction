#!/usr/bin/env python
"""
interface_recovery_worker.py -- standalone heavy/light interface-contact
recovery computation, matching this project's published Table 5 methodology.

Runs in a SEPARATE Python environment from PyMOL, using gemmi (mmCIF I/O)
and Arpeggio (pdbe-arpeggio, interatomic interaction fingerprinting).
Invoked via subprocess from benchmark_master_pipeline.py.

Requires:
  1. conda create -n arpeggio-analysis -c conda-forge python=3.10 numpy
     pandas scipy biopython openbabel
  2. pip install pdbe-arpeggio   (into that same env)
  3. Both the reference and predicted structures must already have their
     heavy/light chains labeled "H" and "L" (i.e. run this AFTER the chain-
     identification step already in benchmark_master_pipeline.py, using the
     same h_chain/l_chain PyMOL selections saved out to temporary H/L-labeled
     PDB files -- see the calling convention in benchmark_master_pipeline.py).

Usage:
    python interface_recovery_worker.py <ref_HL_pdb> <pred_HL_pdb>

Prints one line of JSON to stdout:
    {"fraction_recovered": 0.79, "ref_n_pairs": 34, "pred_n_pairs": 30, "n_shared": 27}
or, on failure:
    {"error": "..."}
"""
import sys
import os
import json
import tempfile
import subprocess

SPECIFIC_TYPES = {"hbond", "weak_hbond", "polar", "weak_polar", "ionic",
                   "aromatic", "hydrophobic", "carbonyl", "cation_pi",
                   "halogen_bond", "metal_complex", "amide_pi", "donor_pi",
                   "halogen_pi", "pi_pi"}


def clip_fv_and_convert(in_pdb, out_cif, max_resi=120):
    import gemmi
    st = gemmi.read_structure(in_pdb)
    st.setup_entities()
    model = st[0]
    for chain in model:
        to_remove = [i for i, res in enumerate(chain) if res.seqid.num > max_resi or res.seqid.num < 1]
        for i in reversed(to_remove):
            del chain[i]
    doc = st.make_mmcif_document()
    doc.write_file(out_cif)


def run_arpeggio(cif_path, workdir):
    result = subprocess.run(
        ["pdbe-arpeggio", "-s", "/L//", "-o", workdir, cif_path],
        capture_output=True, text=True, timeout=90
    )
    json_path = cif_path.replace(".cif", ".json")
    if not os.path.exists(json_path):
        return None, result.stderr[-500:]
    with open(json_path) as f:
        data = json.load(f)
    return data, None


def extract_interface_pairs(data):
    """Return set of (heavy_resi, light_resi) pairs with at least one
    chemically specific interaction type (excludes generic/loose vdw contacts)."""
    pairs = set()
    for entry in data:
        if entry["interacting_entities"] != "INTER":
            continue
        types = set(entry["contact"])
        if not (types & SPECIFIC_TYPES):
            continue
        bgn, end = entry["bgn"], entry["end"]
        h_resi, l_resi = None, None
        for atom in (bgn, end):
            if atom["auth_asym_id"] == "H":
                h_resi = atom["auth_seq_id"]
            elif atom["auth_asym_id"] == "L":
                l_resi = atom["auth_seq_id"]
        if h_resi is not None and l_resi is not None:
            pairs.add((h_resi, l_resi))
    return pairs


def main():
    ref_pdb, pred_pdb = sys.argv[1], sys.argv[2]
    workdir = tempfile.mkdtemp(prefix="ifprint_")
    try:
        ref_cif = os.path.join(workdir, "ref.cif")
        pred_cif = os.path.join(workdir, "pred.cif")
        clip_fv_and_convert(ref_pdb, ref_cif)
        clip_fv_and_convert(pred_pdb, pred_cif)

        ref_data, err = run_arpeggio(ref_cif, workdir)
        if ref_data is None:
            print(json.dumps({"error": f"ref_arpeggio_failed: {err}"}))
            sys.exit(1)
        pred_data, err = run_arpeggio(pred_cif, workdir)
        if pred_data is None:
            print(json.dumps({"error": f"pred_arpeggio_failed: {err}"}))
            sys.exit(1)

        ref_pairs = extract_interface_pairs(ref_data)
        pred_pairs = extract_interface_pairs(pred_data)
        intersect = ref_pairs & pred_pairs
        fraction = len(intersect) / len(ref_pairs) if ref_pairs else None

        print(json.dumps({
            "fraction_recovered": round(fraction, 4) if fraction is not None else None,
            "ref_n_pairs": len(ref_pairs),
            "pred_n_pairs": len(pred_pairs),
            "n_shared": len(intersect),
        }))
    except Exception as e:
        print(json.dumps({"error": str(e)[:300]}))
        sys.exit(1)
    finally:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
