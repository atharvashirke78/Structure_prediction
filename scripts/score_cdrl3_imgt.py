"""IMGT CDR-L3 RMSD, Table 1 protocol, uniform SS normalisation, with a documented
fallback: if the beta-strand-only framework selection is empty on either structure
(observed for 7/58 Ibex light chains), fall back to an unrestricted-Calpha framework
fit (resi 1-120, name CA, no ss filter) on both structures for that pair, rather than
dropping the case. This mirrors the audit done for the ABodyBuilder3/3HC3 CDR-H3 case.
"""
import pymol, json, os, sys
pymol.finish_launching(["pymol","-qc"])
from pymol import cmd

BASE = "/home/indiskaai_laptop_3/claude-science-data/Structure Prediction Paper/All models fv structures"
sel = json.load(open("chk/l3_selections.json"))
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
            l3_ref, l3_pred = info["Reference"]["resi"], info[m]["resi"]
            n_ref_ss = cmd.count_atoms(f"ref and chain {rch} and resi 1-120 and name CA and ss S and not resi {l3_ref}")
            n_pred_ss = cmd.count_atoms(f"pred and chain {pch} and resi 1-120 and name CA and ss S and not resi {l3_pred}")
            fallback = (n_ref_ss == 0 or n_pred_ss == 0)
            # Framework fit always excludes the CDR-L3 residues themselves, so the loop
            # cannot influence the frame in which its own RMSD is subsequently measured.
            if fallback:
                fr = cmd.super(f"pred and chain {pch} and resi 1-120 and name CA and not resi {l3_pred}",
                               f"ref and chain {rch} and resi 1-120 and name CA and not resi {l3_ref}")
            else:
                fr = cmd.super(f"pred and chain {pch} and resi 1-120 and name CA and ss S and not resi {l3_pred}",
                               f"ref and chain {rch} and resi 1-120 and name CA and ss S and not resi {l3_ref}")
            res = cmd.align(f"pred and chain {pch} and resi {info[m]['resi']} and name CA",
                            f"ref and chain {rch} and resi {info['Reference']['resi']} and name CA",
                            transform=0, cycles=0)
            row.update({"L3_IMGT_RMSD": round(res[0],4), "L3_n_aligned": res[1],
                        "framework_rmsd": round(fr[0],4), "framework_n": fr[1],
                        "framework_fallback": fallback,
                        "n_ref_l3": info["Reference"]["n"], "n_pred_l3": info[m]["n"],
                        "Status": "Success"})
        except Exception as e:
            row.update({"Status": "ERR", "Notes": str(e)[:150]})
        out.append(row)
    print("done", t, file=sys.stderr)

json.dump(out, open("l3_imgt_v2.json","w"), indent=1)
print("TOTAL", len(out), file=sys.stderr)
