
import pymol, os, glob, json, re, sys
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd

resns_3to1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G","HIS":"H",
              "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def seq_identity(obj, ch):
    """N-terminal sequence signature check: LIGHT / HEAVY / UNKNOWN.
    Fixed 2026-07: the LIGHT motif check used seq[:8] against 7-character
    motif strings, so it could never match; also expanded the motif lists
    to cover additional observed N-terminal variants (kappa/lambda light,
    and non-canonical heavy)."""
    model = cmd.get_model(f"{obj} and chain {ch} and name CA")
    seq = "".join([resns_3to1.get(a.resn, "X") for a in model.atom])
    head = seq[:15]
    if any(m in head for m in ("TQSPS","TQSPL","TQSPG","TQPPS","TQSPA","TQSLA","TQSPD")):
        return "LIGHT"
    if seq[:7] in ("SDIVMTQ","DIQMTQS","EIVLTQS","QSALTQP","DIQMTQI","DIQMTQT",
                    "SDISVAP","SYVLTQP","QAVVTQP","QPVLTQP","QSVLTQP","NFMLTQP","QTVVTQP"):
        return "LIGHT"
    if head.startswith(("EVQL","QVQL","QVTL","QITL","EVTL","QLQL","EVKL","QMQL","EVMV","QVML",
                         "AVQL","AVKL","DVQL","EIQL","QIQL","QVQV","EVQV")):
        return "HEAVY"
    return "UNKNOWN"

def _get_h_chain(obj, valid_chains):
    """Shared heavy-chain selection, unchanged from the original heuristic."""
    h3_scores = {ch: cmd.count_atoms(f"{obj} and chain {ch} and resi 90-115 and name CA") for ch in valid_chains}
    max_score = max(h3_scores.values())
    tied = [ch for ch, s in h3_scores.items() if s == max_score]
    if len(tied) == 1:
        h_chain = tied[0]
    else:
        heavy_candidates = [ch for ch in tied if seq_identity(obj, ch) == "HEAVY"]
        h_chain = heavy_candidates[0] if heavy_candidates else tied[0]
    if seq_identity(obj, h_chain) == "LIGHT":
        heavy_alt = [ch for ch in valid_chains if ch != h_chain and seq_identity(obj, ch) == "HEAVY"]
        if heavy_alt:
            h_chain = heavy_alt[0]
    return h_chain

def _get_chains_HL_ref(obj):
    """Reference-structure chain ID. Fixed 2026-07: light-chain selection
    now picks the MAXIMUM-contact candidate (was: first candidate clearing
    a fixed threshold), correctly resolving multi-copy asymmetric units."""
    chains = cmd.get_chains(obj)
    if not chains: return None, None
    sizes = {ch: cmd.count_atoms(f"{obj} and chain {ch} and name CA") for ch in chains}
    valid_chains = [ch for ch, size in sizes.items() if size > 60]
    if not valid_chains: return None, None
    h_chain = _get_h_chain(obj, valid_chains)
    l_chain = None
    l_candidates = [ch for ch in valid_chains if ch != h_chain and sizes[ch] < 300]
    best_contact = -1
    for ch in l_candidates:
        contact_atoms = cmd.count_atoms(f"({obj} and chain {ch}) within 4.0 of ({obj} and chain {h_chain})")
        if contact_atoms > best_contact:
            best_contact = contact_atoms
            l_chain = ch
    if best_contact <= 10:
        l_chain = None
    return h_chain, l_chain

def _get_chains_HL_pred(obj):
    """Prediction-structure chain ID. Fixed 2026-07: light-chain selection
    checks sequence identity FIRST for all candidates; contact-distance
    (widened to 8.0A) is only a tiebreak among multiple sequence-confirmed
    LIGHT candidates, or a fallback (legacy 4.0A/>10 rule) when none reads
    LIGHT. Fixes IgFold's light chain being systematically undetected
    (its H-L interface registers only 4-9 CA contact atoms within 4.0A,
    below the old >10 threshold, for 100% of predictions tested)."""
    chains = cmd.get_chains(obj)
    if not chains: return None, None
    sizes = {ch: cmd.count_atoms(f"{obj} and chain {ch} and name CA") for ch in chains}
    valid_chains = [ch for ch, size in sizes.items() if size > 60]
    if not valid_chains: return None, None
    h_chain = _get_h_chain(obj, valid_chains)
    l_candidates = [ch for ch in valid_chains if ch != h_chain and sizes[ch] < 300]
    light_by_seq = [ch for ch in l_candidates if seq_identity(obj, ch) == "LIGHT"]
    l_chain = None
    if len(light_by_seq) == 1:
        l_chain = light_by_seq[0]
    elif len(light_by_seq) > 1:
        best_contact = -1
        for ch in light_by_seq:
            contact_atoms = cmd.count_atoms(f"({obj} and chain {ch}) within 8.0 of ({obj} and chain {h_chain})")
            if contact_atoms > best_contact:
                best_contact = contact_atoms
                l_chain = ch
    else:
        best_fv_score = -1
        for ch in l_candidates:
            contact_atoms = cmd.count_atoms(f"({obj} and chain {ch}) within 4.0 of ({obj} and chain {h_chain})")
            if contact_atoms > 10:
                fv_score = cmd.count_atoms(f"{obj} and chain {ch} and resi 1-120 and name CA")
                if fv_score > best_fv_score:
                    best_fv_score = fv_score
                    l_chain = ch
    return h_chain, l_chain

base = "/home/indiskaai_laptop_3/claude-science-data/Structure Prediction Paper/Outputs"
model_dirs = {
    "AB_Fold_3": os.path.join(base, "AB Fold 3 outputs"),
    "ABodyBuilder3": os.path.join(base, "ABodyBuilder3 outputs"),
    "Alpha_Fold_2": os.path.join(base, "Alpha Fold 2 outputs"),
    "Boltz_2": os.path.join(base, "Boltz 2 outputs"),
    "ESM_2": os.path.join(base, "ESM 2 outputs"),
    "Ibex": os.path.join(base, "Ibex outputs"),
    "IgFold": os.path.join(base, "IgFold outputs"),
    "Intellifold": os.path.join(base, "Intellifold outputs"),
}

def find_file_for_pdb(model_dir, pdb_id):
    all_files = glob.glob(os.path.join(model_dir, "**", "*.pdb"), recursive=True) + \
                glob.glob(os.path.join(model_dir, "**", "*.cif"), recursive=True)
    candidates = []
    for filepath in all_files:
        if "zone.identifier" in filepath.lower():
            continue
        filename = os.path.basename(filepath).lower()
        if "quarantine" in filepath.lower():
            continue
        if "rank" in filename and not any(best in filename for best in ["ranked_0", "rank_1", "rank_001"]):
            continue
        m = re.search(r'([1-9][a-zA-Z0-9]{3})', filepath, re.IGNORECASE)
        if m and m.group(1).upper() == pdb_id:
            candidates.append(filepath)
    candidates.sort()
    return candidates[0] if candidates else None

pdb_ids = json.loads(sys.argv[1])
results = {}

for model_name, model_dir in model_dirs.items():
    model_results = {}
    for pdb_id in pdb_ids:
        filepath = find_file_for_pdb(model_dir, pdb_id)
        if not filepath:
            model_results[pdb_id] = {"status": "NO_FILE"}
            continue
        cmd.reinitialize()
        entry = {"status": "Success"}
        try:
            cmd.load(filepath, "pred_clean")
            cmd.fetch(pdb_id, "ref_raw", type="pdb")
            cmd.create("ref_clean", "ref_raw and polymer")
            cmd.delete("ref_raw")
            for ch in cmd.get_chains("ref_clean"):
                n = cmd.count_atoms(f"ref_clean and chain {ch} and name CA")
                if n > 300:
                    cmd.remove(f"ref_clean and chain {ch}")
            cmd.remove("pred_clean and not polymer")

            ref_h, ref_l = _get_chains_HL_ref("ref_clean")
            pred_h, pred_l = _get_chains_HL_pred("pred_clean")
            entry["ref_h"] = ref_h; entry["ref_l"] = ref_l
            entry["pred_h"] = pred_h; entry["pred_l"] = pred_l

            # Masked Fv RMSD (combined chain selection, sanity check -- should be near-identical to original)
            ref_fv_sel = f"ref_clean and chain {ref_h}+{ref_l} and resi 1-120" if ref_l else f"ref_clean and chain {ref_h} and resi 1-120"
            pred_fv_sel = f"pred_clean and chain {pred_h}+{pred_l} and resi 1-120" if pred_l else f"pred_clean and chain {pred_h} and resi 1-120"
            try:
                fv_result_masked = cmd.super(pred_fv_sel, ref_fv_sel)
                if fv_result_masked[1] < 10:
                    entry["masked_fv"] = "ALIGN_FAIL"
                else:
                    entry["masked_fv"] = round(fv_result_masked[0], 4)
                    fv_result_zero = cmd.super(pred_fv_sel, ref_fv_sel, cycles=0, transform=0)
                    entry["zero_fv"] = round(fv_result_zero[0], 4)
            except Exception as e:
                entry["masked_fv"] = "ERR:" + str(e)[:60]

            # CDR-H3 (corrected chain identity)
            if not ref_h or not pred_h:
                entry["cdr_h3"] = "H_chain_missing"
            elif cmd.count_atoms(f"ref_clean and chain {ref_h} and resi 90-120 and name CA") < 5:
                entry["cdr_h3"] = "H3_unresolved_in_ref"
            else:
                try:
                    cmd.super(f"pred_clean and chain {pred_h} and resi 1-120 and ss s",
                              f"ref_clean and chain {ref_h} and resi 1-120 and ss s")
                    h3_result = cmd.align(
                        f"pred_clean and chain {pred_h} and resi 90-120 and name CA",
                        f"ref_clean  and chain {ref_h}  and resi 90-120 and name CA",
                        transform=0, cycles=0
                    )
                    entry["cdr_h3"] = round(h3_result[0], 4)
                except Exception as e:
                    entry["cdr_h3"] = "ERR:" + str(e)[:60]

            # CDR-L3 (corrected chain identity)
            if not ref_l or not pred_l:
                entry["cdr_l3"] = "L_chain_missing"
            elif cmd.count_atoms(f"ref_clean and chain {ref_l} and resi 89-105 and name CA") < 3:
                entry["cdr_l3"] = "L3_unresolved_in_ref"
            else:
                try:
                    cmd.super(f"pred_clean and chain {pred_l} and resi 1-120 and ss s",
                              f"ref_clean and chain {ref_l} and resi 1-120 and ss s")
                    l3_result = cmd.align(
                        f"pred_clean and chain {pred_l} and resi 89-105 and name CA",
                        f"ref_clean  and chain {ref_l}  and resi 89-105 and name CA",
                        transform=0, cycles=0
                    )
                    entry["cdr_l3"] = round(l3_result[0], 4)
                except Exception as e:
                    entry["cdr_l3"] = "ERR:" + str(e)[:60]

        except Exception as e:
            entry["status"] = "CRASH: " + str(e)[:100]
        model_results[pdb_id] = entry
    results[model_name] = model_results
    print(f"done: {model_name}", file=sys.stderr)

print(json.dumps(results))
