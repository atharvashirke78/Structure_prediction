"""IMGT-defined CDR-H3 RMSD, Table 1 protocol, with uniform SS normalisation.

Framework fit: beta-strand Calpha, residues 1-120, heavy chain (as Table 1).
Loop RMSD: ANARCI/IMGT CDR-H3 (positions 105-117) only, transform=0, cycles=0.
Uniform preprocessing: explicit hydrogens removed and cmd.dss() re-run on BOTH
reference and prediction before selection, so PyMOL's strand assignment is
computed from heavy atoms identically for every model (one structure,
ABodyBuilder3/3HC3, otherwise receives no strand assignment at all).
"""
import pymol, json, os, sys
pymol.finish_launching(["pymol","-qc"])
from pymol import cmd

BASE = "/home/indiskaai_laptop_3/claude-science-data/Structure Prediction Paper/All models fv structures"
sel = json.load(open("chk/h3_selections.json"))
MODELS = ["ABodyBuilder3","AbFold","AlphaFold2","Boltz2","Chai1","ESMFold2","Ibex","IgFold","IntelliFold","OpenFold3"]

def path_for(m, t):
    for e in (".pdb", ".cif"):
        p = os.path.join(BASE, m, t + e)
        if os.path.exists(p): return p

out = []
for t, info in sorted(sel.items()):
    rp = path_for("Reference", t)
    for m in MODELS:
        pp = path_for(m, t)
        row = {"Model": m, "PDB_ID": t}
        try:
            cmd.reinitialize()
            cmd.load(rp, "ref"); cmd.load(pp, "pred")
            cmd.remove("hydro")
            cmd.dss("ref"); cmd.dss("pred")
            rch, pch = info["Reference"]["chain"], info[m]["chain"]
            fr = cmd.super(f"pred and chain {pch} and resi 1-120 and name CA and ss S",
                           f"ref and chain {rch} and resi 1-120 and name CA and ss S")
            res = cmd.align(f"pred and chain {pch} and resi {info[m]['resi']} and name CA",
                            f"ref and chain {rch} and resi {info['Reference']['resi']} and name CA",
                            transform=0, cycles=0)
            row.update({"H3_IMGT_RMSD": round(res[0],4), "H3_n_aligned": res[1],
                        "framework_rmsd": round(fr[0],4), "framework_n": fr[1],
                        "n_ref_h3": info["Reference"]["n"], "n_pred_h3": info[m]["n"],
                        "Status": "Success"})
        except Exception as e:
            row.update({"Status": "ERR", "Notes": str(e)[:150]})
        out.append(row)
    print("done", t, file=sys.stderr)

json.dump(out, open("h3_imgt_v2.json","w"), indent=1)
print("TOTAL", len(out), file=sys.stderr)
