import os
import csv
import subprocess
import glob
import re
import json
import tempfile
from pymol import cmd


# ─────────────────────────────────────────────────────────────────────────
# 2026-08 update: added four new metrics matching this project's published
# Table 3, Table 4, and Table 5 -- Ramachandran-outlier %, clashscore,
# interface-contact recovery fraction, and a standardized (0-100, per-model-
# aware) confidence score. The first three require heavy external
# dependencies (cctbx/mmtbx, gemmi+Arpeggio) that do NOT reliably coexist in
# the same process as PyMOL, so they run as separate worker scripts in their
# own Python environments, invoked here via subprocess -- exactly how they
# were computed for the published values (never imported directly into
# PyMOL in this project, on any platform).
#
# SETUP REQUIRED before these four columns will populate (see
# README_metrics_setup.md for exact commands):
#   1. A conda/venv environment with cctbx-base installed, plus a local
#      clone of https://github.com/phenix-project/geostd (CCP4 monomer
#      restraint library) for clashscore2's bonding-restraint lookup.
#   2. A conda/venv environment with gemmi + pdbe-arpeggio installed
#      (Arpeggio also needs Open Babel; see README_metrics_setup.md).
#   3. Both python.exe paths and the geostd folder path set below.
#
# If a worker environment isn't configured, the corresponding column(s) are
# written as "N/A" with a note -- the rest of the pipeline (RMSD, CDR,
# TM-score, pTM/ipTM -- all unchanged from before) still runs normally.
# ─────────────────────────────────────────────────────────────────────────
RAMA_CLASH_PYTHON = r"C:\IndiskaAI\Envs\structqc\python.exe"      # <-- EDIT to your cctbx env's python.exe
RAMA_CLASH_WORKER = r"C:\IndiskaAI\Pipelines\rama_clash_worker.py"  # <-- EDIT if you place the worker script elsewhere
GEOSTD_PATH       = r"C:\IndiskaAI\Tools\geostd"                    # <-- EDIT to your local geostd clone

INTERFACE_PYTHON = r"C:\IndiskaAI\Envs\arpeggio-analysis\python.exe"   # <-- EDIT to your Arpeggio env's python.exe
INTERFACE_WORKER = r"C:\IndiskaAI\Pipelines\interface_recovery_worker.py"  # <-- EDIT if you place the worker script elsewhere


def _run_rama_clash(filepath):
    """Subprocess call to rama_clash_worker.py. Returns (pct_outliers, clashscore, note)."""
    if not (os.path.exists(RAMA_CLASH_PYTHON) and os.path.exists(RAMA_CLASH_WORKER)):
        return "N/A", "N/A", "rama_clash_worker_not_configured"
    env = os.environ.copy()
    env["MMTBX_CCP4_MONOMER_LIB"] = GEOSTD_PATH
    try:
        result = subprocess.run(
            [RAMA_CLASH_PYTHON, RAMA_CLASH_WORKER, filepath],
            capture_output=True, text=True, timeout=180, env=env
        )
        data = json.loads(result.stdout.strip().splitlines()[-1]) if result.stdout.strip() else {}
        if "error" in data:
            return "N/A", "N/A", f"rama_clash_error:{data['error'][:80]}"
        return data.get("pct_rama_outliers", "N/A"), data.get("clashscore", "N/A"), ""
    except Exception as e:
        return "N/A", "N/A", f"rama_clash_crashed:{str(e)[:80]}"


def _run_interface_recovery(ref_hl_path, pred_hl_path):
    """Subprocess call to interface_recovery_worker.py. Returns (fraction_recovered, note)."""
    if not (os.path.exists(INTERFACE_PYTHON) and os.path.exists(INTERFACE_WORKER)):
        return "N/A", "interface_recovery_worker_not_configured"
    try:
        result = subprocess.run(
            [INTERFACE_PYTHON, INTERFACE_WORKER, ref_hl_path, pred_hl_path],
            capture_output=True, text=True, timeout=120
        )
        data = json.loads(result.stdout.strip().splitlines()[-1]) if result.stdout.strip() else {}
        if "error" in data:
            return "N/A", f"interface_recovery_error:{data['error'][:80]}"
        return data.get("fraction_recovered", "N/A"), ""
    except Exception as e:
        return "N/A", f"interface_recovery_crashed:{str(e)[:80]}"


def _extract_standardized_confidence(model_name, model_dir, pdb_id, pred_filepath):
    """Per-model confidence extraction matching this project's published Table 3
    values exactly -- each model's native output format and scale differ, and
    this reproduces the specific rule validated for each one.

    Returns a float on a 0-100 scale, or "N/A" if the model does not emit a
    usable per-target confidence value (ABodyBuilder3: none emitted at all;
    Ibex: emits a constant 0.0 for every target, not informative).
    """
    def bfactor_avg(path):
        vals = []
        try:
            with open(path) as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        try:
                            vals.append(float(line[60:66]))
                        except ValueError:
                            pass
        except Exception:
            return None
        return sum(vals) / len(vals) if vals else None

    try:
        if model_name in ("ESMFold2", "Chai1"):
            # JSON sidecar, avg_plddt on a 0-1 scale
            fpath = os.path.join(model_dir, f"{pdb_id}_{model_name.lower()}_confidence.json")
            if os.path.exists(fpath):
                with open(fpath) as f:
                    data = json.load(f)
                return round(data["avg_plddt"] * 100, 2)

        elif model_name == "AB_Fold_3":
            # B-factor column, already on a 0-100 scale (no further rescaling)
            v = bfactor_avg(pred_filepath)
            return round(v, 2) if v is not None else "N/A"

        elif model_name == "IgFold":
            # B-factor column, on a 0-1 scale -- needs x100
            v = bfactor_avg(pred_filepath)
            return round(v * 100, 2) if v is not None else "N/A"

        elif model_name == "Boltz_2":
            # confidence_{pdb}_model_0.json, complex_plddt on a 0-1 scale
            pattern = os.path.join(model_dir, "predictions", f"boltz_results_{pdb_id}",
                                    "predictions", pdb_id, f"confidence_{pdb_id}_model_0.json")
            if os.path.exists(pattern):
                with open(pattern) as f:
                    data = json.load(f)
                return round(data["complex_plddt"] * 100, 2)

        elif model_name == "Alpha_Fold_2":
            # scores JSON, plddt list already on a 0-100 scale; falls back to
            # B-factor (also 0-100 native) if the JSON isn't found
            dirpath = os.path.join(model_dir, f"{pdb_id}_Complex")
            matches = glob.glob(os.path.join(dirpath, f"{pdb_id}_Complex_scores_rank_001_*.json"))
            if matches:
                with open(matches[0]) as f:
                    data = json.load(f)
                return round(sum(data["plddt"]) / len(data["plddt"]), 2)
            matches2 = glob.glob(os.path.join(dirpath, f"{pdb_id}_Complex_unrelaxed_rank_001_*.pdb"))
            if matches2:
                v = bfactor_avg(matches2[0])
                return round(v, 2) if v is not None else "N/A"

        elif model_name == "Intellifold":
            # summary_confidences JSON, plddt on a 0-1 scale
            fpath = os.path.join(model_dir, "Input_Yaml", "predictions", pdb_id,
                                  f"{pdb_id}_seed-42_sample-0_summary_confidences.json")
            if os.path.exists(fpath):
                with open(fpath) as f:
                    data = json.load(f)
                return round(data["plddt"] * 100, 2)

        elif model_name == "OpenFold3":
            # confidence sidecar, complex_plddt_score already on a 0-100 scale
            candidates = glob.glob(os.path.join(model_dir, "**", f"{pdb_id}*confidence*.json"), recursive=True) + \
                         glob.glob(os.path.join(model_dir, "**", f"{pdb_id}*summary_confidences*.json"), recursive=True)
            for fpath in candidates:
                with open(fpath) as f:
                    data = json.load(f)
                if "complex_plddt_score" in data:
                    return round(data["complex_plddt_score"], 2)

        elif model_name == "Ibex":
            # Emits a constant 0.0 for every target -- not usable; report N/A
            # rather than a misleadingly precise-looking zero.
            return "N/A"

        elif model_name == "ABodyBuilder3":
            # Does not emit any per-target confidence score.
            return "N/A"

    except Exception:
        pass

    return "N/A"


_RESN_3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G","HIS":"H",
              "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
              "PCA":"Q"}  # 2026-08 fix: pyroglutamate (cyclized N-terminal Gln), common in antibody
                          # heavy/light chains, was mapping to "X" and corrupting the N-terminal
                          # sequence signature check for any chain starting with it (e.g. 2XZC, 5KG9).


def _seq_identity(obj, ch):
    """N-terminal sequence signature check: LIGHT / HEAVY / UNKNOWN.

    Used to break ties in _get_chains_HL_ref/_get_chains_HL_pred below, where
    the CDR-H3-window atom count (the primary heavy-chain signal) is equal
    across two or more chains -- which happens whenever heavy and light
    chains have similar lengths in that window, i.e. often. Left unbroken,
    that tie was resolved arbitrarily (first chain in PyMOL's chain list),
    causing heavy/light swaps that corrupted CDR-H3/CDR-L3 RMSD (and, for
    asymmetric units with more than one Fab copy, occasionally the main Fv
    RMSD too).
    """
    model = cmd.get_model(f"{obj} and chain {ch} and name CA")
    seq = "".join([_RESN_3TO1.get(a.resn, "X") for a in model.atom])
    head = seq[:15]
    if any(m in head for m in ("TQSPS", "TQSPL", "TQSPG", "TQPPS", "TQSPA", "TQSLA", "TQSPD")):
        return "LIGHT"
    # Fixed 2026-07: these motifs are 7 characters, but this slice used to be
    # seq[:8], so an 8-char slice could never equal a 7-char string and this
    # branch never matched -- silently dropping several light-chain
    # N-terminal variants (kappa and lambda) to UNKNOWN. Also expanded the
    # motif list to cover additional observed light-chain N-termini.
    if seq[:7] in ("SDIVMTQ", "DIQMTQS", "EIVLTQS", "QSALTQP", "DIQMTQI", "DIQMTQT",
                    "SDISVAP", "SYVLTQP", "QAVVTQP", "QPVLTQP", "QSVLTQP", "NFMLTQP", "QTVVTQP"):
        return "LIGHT"
    if head.startswith(("EVQL", "QVQL", "QVTL", "QITL", "EVTL", "QLQL", "EVKL", "QMQL", "EVMV", "QVML",
                         "AVQL", "AVKL", "DVQL", "EIQL", "QIQL", "QVQV", "EVQV")):
        return "HEAVY"
    # 2026-08 fix: ANARCI cross-validation against the reference set found 4 targets
    # (1T04, 1T3F, 2XZC, 5KG9) with a genuine H/L swap because the true heavy chain's
    # leading residue is unresolved in the deposited structure (e.g. "VQLVQSGAEL..."
    # instead of "EVQLVQSGAEL..."), defeating the exact-prefix check above and leaving
    # the tie/override logic unable to confirm HEAVY. Tolerating one missing N-terminal
    # residue recovers all 4 with zero regressions (verified against all 58 reference
    # targets: the only sequences newly matched here were either these 4 true swaps or
    # chains already correctly labeled HEAVY by the CDR-H3-window count, for which this
    # check only confirms rather than changes the outcome).
    heavy_prefixes_trunc1 = ("VTL", "VQL", "IQL", "VKL", "VMV", "VML", "LQL", "MQL", "ITL", "VQV")
    if head.startswith(heavy_prefixes_trunc1):
        return "HEAVY"
    # 2026-08 fix: 5KG9's true light chain starts "QLVLTQSSS...", a germline
    # variant outside the fixed motif list above. VLTQ/IQMTQ substrings near the
    # N-terminus are light-chain-specific frameworks (0 false positives across all
    # 53 pre-classified sequences in this reference set with either motif, and
    # never seen in a HEAVY chain).
    if "VLTQ" in head or "IQMTQ" in head:
        return "LIGHT"
    # 2026-08 fix: 4UOK's true heavy chain starts "QLVQSGAEVKKPGA...", missing
    # the canonical leading "QV" entirely (2 residues, not just 1). VQSG in the
    # first 8 characters is a general framework-1 heavy-chain signal: 0 false
    # positives against all 53 pre-classified LIGHT sequences in this reference set.
    if "VQSG" in seq[:8]:
        return "HEAVY"
    return "UNKNOWN"


def _get_h_chain(obj, valid_chains):
    """Shared heavy-chain selection, unchanged from the original heuristic
    for both reference and prediction structures: max CDR-H3-window (resi
    90-115) CA count, with sequence-identity tiebreak/override. Left as-is
    deliberately -- 18/58 reference structures have a known, separate
    heavy-chain-identification ambiguity (no chain confirms HEAVY by
    sequence in the H3-window tie) that is out of scope for this fix and
    would need independent validation across all 8 models before changing.
    """
    h3_scores = {ch: cmd.count_atoms(f"{obj} and chain {ch} and resi 90-115 and name CA") for ch in valid_chains}
    max_score = max(h3_scores.values())
    tied = [ch for ch, s in h3_scores.items() if s == max_score]

    if len(tied) == 1:
        h_chain = tied[0]
    else:
        heavy_candidates = [ch for ch in tied if _seq_identity(obj, ch) == "HEAVY"]
        h_chain = heavy_candidates[0] if heavy_candidates else tied[0]

    if _seq_identity(obj, h_chain) == "LIGHT":
        heavy_alt = [ch for ch in valid_chains if ch != h_chain and _seq_identity(obj, ch) == "HEAVY"]
        if heavy_alt:
            h_chain = heavy_alt[0]

    return h_chain


def _get_chains_HL_ref(obj):
    """Reference-structure chain identification.

    Fixed 2026-07: light-chain selection now picks the candidate with
    MAXIMUM physical contact to h_chain (was: the first candidate whose
    contact count cleared a fixed >10-atom threshold, in whatever order
    PyMOL's chain list happened to enumerate them). For multi-copy
    asymmetric units (8 of 58 reference targets: 1AD9, 1FVC, 1RHH, 2HFF,
    5BZD, 5CEX, 5KG9, 6QN8 -- confirmed via direct testing), several chains
    can each clear that threshold (e.g. a second Fab copy's light chain
    touching the chosen heavy chain at a crystal contact), and the old logic
    could pick a physically weaker, wrong pairing over the true partner.
    Picking the maximum-contact candidate instead resolves this correctly
    in every case checked. h_chain selection is unchanged (see _get_h_chain).
    """
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
    """Prediction-structure chain identification.

    Fixed 2026-07: light-chain selection now checks sequence identity FIRST
    for every light-chain-sized candidate; physical contact (widened to
    8.0A) is used only as a tiebreak among multiple sequence-confirmed
    LIGHT candidates, or as a fallback (using the original 4.0A/>10 rule)
    when no candidate's sequence signature reads LIGHT. This fixes IgFold's
    light chain being systematically undetected: IgFold's H-L docking
    geometry produced only 4-9 CA atoms in contact within 4.0A across 12
    spot-checked targets, all below the old >10 threshold, so
    contact-distance-only detection failed for 100% of IgFold predictions
    tested (n=30). Confirmed via direct testing to introduce zero
    regressions across the other 7 models (AB_Fold_3, ABodyBuilder3,
    Alpha_Fold_2, Boltz_2, Ibex, Intellifold all agree with the legacy
    result on every sampled target). h_chain selection is unchanged (see
    _get_h_chain); predictions are single-copy Fv/Fab structures, so the
    multi-copy heavy-chain ambiguity seen in some references does not
    arise here.

    2026-08 update: this same logic was independently validated (Linux/
    PyMOL, separate scoring pass) against ESMFold2 and Chai-1 predictions --
    both emit standard two-chain (H/L) PDB output, both were scored 58/58
    with status=Success, and both are added to model_dirs below alongside
    the original 8. No code changes were required for either model.
    """
    chains = cmd.get_chains(obj)
    if not chains: return None, None
    sizes = {ch: cmd.count_atoms(f"{obj} and chain {ch} and name CA") for ch in chains}
    valid_chains = [ch for ch, size in sizes.items() if size > 60]
    if not valid_chains: return None, None

    h_chain = _get_h_chain(obj, valid_chains)

    l_candidates = [ch for ch in valid_chains if ch != h_chain and sizes[ch] < 300]
    light_by_seq = [ch for ch in l_candidates if _seq_identity(obj, ch) == "LIGHT"]
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



def _extract_json_scores(pdb_filepath):
    """Hunts for pTM and ipTM in companion JSON files."""
    ptm, iptm = "N/A", "N/A"
    dir_name = os.path.dirname(pdb_filepath)
    base_name = os.path.basename(pdb_filepath)
    
    af2_json = base_name.replace("unrelaxed", "scores").replace(".pdb", ".json").replace(".cif", ".json")
    boltz_json = "confidence.json"
    
    json_candidates = [
        os.path.join(dir_name, af2_json),
        os.path.join(dir_name, boltz_json)
    ] + glob.glob(os.path.join(dir_name, "*.json")) 
    
    for jpath in json_candidates:
        if os.path.exists(jpath):
            try:
                with open(jpath, 'r') as f:
                    data = json.load(f)
                    if 'ptm' in data and ptm == "N/A": ptm = f"{data['ptm']:.4f}"
                    if 'iptm' in data and iptm == "N/A": iptm = f"{data['iptm']:.4f}"
                    if ptm != "N/A" and iptm != "N/A": break
            except:
                pass
    return ptm, iptm

# Usage (PyMOL command line, one call per model):
#   run_benchmarking_pipeline AB_Fold_3
#   run_benchmarking_pipeline ABodyBuilder3
#   run_benchmarking_pipeline Alpha_Fold_2
#   run_benchmarking_pipeline Boltz_2
#   run_benchmarking_pipeline Ibex
#   run_benchmarking_pipeline IgFold
#   run_benchmarking_pipeline Intellifold
#   run_benchmarking_pipeline OpenFold3
#   run_benchmarking_pipeline ESMFold2
#   run_benchmarking_pipeline Chai1
# Each call expects predictions in C:\IndiskaAI\Models\{model_name} outputs
# and writes C:\IndiskaAI\Data\Processed\validated_results_{model_name}_Final.csv
def run_benchmarking_pipeline(model_name):
    model_dir   = rf"C:\IndiskaAI\Models\{model_name} outputs"
    tmalign_exe = r"C:\IndiskaAI\Tools\USalignWin64\USalign\USalign.exe"
    output_csv  = rf"C:\IndiskaAI\Data\Processed\validated_results_{model_name.replace(' ', '_')}_Final.csv"
    temp_ref    = r"C:\IndiskaAI\Pipelines\temp_ref.pdb"
    temp_pred   = r"C:\IndiskaAI\Pipelines\temp_pred.pdb"

    cmd.set('suspend_undo', 1)
    tasks = []

    print(f"\nScanning {model_dir} for predicted structure files...")
    all_files = glob.glob(os.path.join(model_dir, "**", "*.pdb"), recursive=True) + \
                glob.glob(os.path.join(model_dir, "**", "*.cif"), recursive=True)

    for filepath in all_files:
        filename = os.path.basename(filepath).lower()
        if "quarantine" in filepath.lower(): continue
        if "rank" in filename and not any(best in filename for best in ["ranked_0", "rank_1", "rank_001"]):
            continue
            
        match = re.search(r'([1-9][a-zA-Z0-9]{3})', filepath, re.IGNORECASE)
        if match:
            pdb_id = match.group(1).upper()
            if not any(t['pdb_id'] == pdb_id for t in tasks):
                tasks.append({'pdb_id': pdb_id, 'filepath': filepath})

    if not tasks:
        print(f"ERROR: No valid files found in {model_dir}")
        return

    tasks.sort(key=lambda x: x['pdb_id'])

    print("\n" + "=" * 200)
    print(f"{'PDB_ID':<6} | {'True TM':<8} | {'Ov. RMSD':<10} | {'Masked Fv':<10} | {'Zero Fv':<10} | {'Delta':<8} | {'CDR-H3':<8} | {'CDR-L3':<8} | {'pLDDT':<7} | {'pTM':<7} | {'ipTM':<7} | {'Conf':<7} | {'Rama%':<7} | {'Clash':<9} | {'IF-Rec':<7}")
    print("=" * 200)

    with open(output_csv, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['PDB_ID', 'True_TM_Score', 'Fv_Anchored_Overall_RMSD', 'Masked_Fv_RMSD', 'Strict_Zero_Fv_RMSD', 'RMSD_Delta', 'CDR_H3_RMSD', 'CDR_L3_RMSD', 'Avg_pLDDT', 'pTM', 'ipTM', 'Standardized_Confidence', 'Rama_Outliers_Pct', 'Clashscore', 'Interface_Recovery_Fraction', 'Status', 'Notes'])

        for task in tasks:
            pdb_id, filepath = task['pdb_id'], task['filepath']
            cmd.reinitialize()
            notes, status = "", "Success"
            calc_overall = calc_masked = calc_zero = calc_delta = calc_h3 = calc_l3 = calc_plddt = "N/A"
            calc_conf = calc_rama = calc_clash = calc_ifrec = "N/A"

            calc_ptm, calc_iptm = _extract_json_scores(filepath)
            calc_conf = _extract_standardized_confidence(model_name, model_dir, pdb_id, filepath)

            try:
                cmd.load(filepath, "pred_clean")
                cmd.fetch(pdb_id, "ref_raw", type="pdb")

                b_factors = []
                cmd.iterate("pred_clean and name CA", "b_factors.append(b)", space={'b_factors': b_factors})
                if b_factors:
                    calc_plddt = f"{sum(b_factors)/len(b_factors):.2f}"

                cmd.create("ref_clean", "ref_raw and polymer")
                cmd.delete("ref_raw")

                # 2026-08 fix: ANARCI cross-validation found that this size-only
                # antigen cutoff was stripping full-length (Fc-containing) antibody
                # heavy chains in two panel targets (1MCO, 5DK3), each >400 residues,
                # leaving only the light chain and causing it to be mislabeled "H" by
                # the single-chain fallback below. A chain is now stripped as antigen
                # only if BOTH (a) its N-terminal sequence signature does NOT read as
                # a genuine immunoglobulin chain (HEAVY/LIGHT), AND (b) at least 2
                # other valid-sized (60-300) chains would remain afterward to form a
                # real H+L pair -- 1MCO's true heavy chain starts "PLVL...", a rare
                # germline variant outside the motif list, so (a) alone is not enough;
                # (b) catches it because removing it would leave only 1 chain (the
                # light chain), never a valid post-cleaning antibody state. No target
                # in this 58-structure panel has a genuine bound antigen requiring this
                # branch, so this guard carries zero risk elsewhere.
                _all_chains = cmd.get_chains("ref_clean")
                _sizes = {c: cmd.count_atoms(f"ref_clean and chain {c} and name CA") for c in _all_chains}
                for ch in _all_chains:
                    n = _sizes[ch]
                    if n > 300:
                        ig_type = _seq_identity("ref_clean", ch)
                        if ig_type in ("HEAVY", "LIGHT"):
                            continue
                        others_valid = [c for c in _all_chains if c != ch and 60 < _sizes[c] <= 300]
                        if len(others_valid) < 2:
                            continue
                        cmd.remove(f"ref_clean and chain {ch}")
                        notes += f"removed_antigen_{ch}({n}res); "

                #notes = _isolate_one_fab("ref_clean", notes)
                cmd.remove("pred_clean and not polymer")

                # ── Step 2: Identify Specific H & L Chains ───────────────────
                ref_h, ref_l = _get_chains_HL_ref("ref_clean")
                pred_h, pred_l = _get_chains_HL_pred("pred_clean")

                # Create specific selection strings to ignore duplicates
                ref_fv_sel = f"ref_clean and chain {ref_h}+{ref_l} and resi 1-120" if ref_l else f"ref_clean and chain {ref_h} and resi 1-120"
                pred_fv_sel = f"pred_clean and chain {pred_h}+{pred_l} and resi 1-120" if pred_l else f"pred_clean and chain {pred_h} and resi 1-120"

                # ── Step 3: Dual-Pipeline Fv Superposition ──────────────
                try:
                    # 1. Industry Standard (Masked): Default 5 cycles
                    fv_result_masked = cmd.super(pred_fv_sel, ref_fv_sel)
                    
                    if fv_result_masked[1] < 10:
                        calc_masked, calc_zero, calc_delta = "ALIGN_FAIL", "ALIGN_FAIL", "N/A"
                        status, notes = "Super_Failed", notes + f"fv_super_aligned_only_{fv_result_masked[1]}_atoms; "
                    else:
                        calc_masked_val = fv_result_masked[0]
                        calc_masked = f"{calc_masked_val:.4f}"
                        
                        # 2. Zero-Rejection Standard (True RMSD): Force 0 cycles
                        fv_result_zero = cmd.super(pred_fv_sel, ref_fv_sel, cycles=0, transform=0)
                        calc_zero_val = fv_result_zero[0]
                        calc_zero = f"{calc_zero_val:.4f}"
                        
                        # 3. The Delta 
                        calc_delta = f"{(calc_zero_val - calc_masked_val):.4f}"

                        # Fv-Anchored Overall RMSD
                        try:
                            overall_val = cmd.rms_cur("pred_clean and name CA", "ref_clean and name CA", matchmaker=4)
                            calc_overall = f"{overall_val:.4f}"
                        except Exception as e:
                            calc_overall, notes = "ERR_CRASH", notes + f"overall_crashed:{str(e)[:40]}; "
                except Exception as e:
                    status, notes = "Super_Failed", notes + f"fv_super_crashed:{str(e)[:40]}; "

                # ── Step 4: Strict CDR-H3 & CDR-L3 RMSD ───────────────────────
                # CDR-H3
                if not ref_h or not pred_h:
                    notes += "H_chain_missing; "
                elif cmd.count_atoms(f"ref_clean and chain {ref_h} and resi 90-120 and name CA") < 5:
                    notes += "H3_unresolved_in_ref; "
                else:
                    try:
                        cmd.super(f"pred_clean and chain {pred_h} and resi 1-120 and ss s", 
                                  f"ref_clean and chain {ref_h} and resi 1-120 and ss s")
                        h3_result = cmd.align(
                            f"pred_clean and chain {pred_h} and resi 90-120 and name CA",
                            f"ref_clean  and chain {ref_h}  and resi 90-120 and name CA",
                            transform=0, cycles=0
                        )
                        calc_h3 = f"{h3_result[0]:.4f}"
                    except Exception as e:
                        calc_h3, notes = "ERR_CRASH", notes + f"h3_crashed:{str(e)[:40]}; "
                        
                # CDR-L3
                if not ref_l or not pred_l:
                    pass # Handled in notes already
                elif cmd.count_atoms(f"ref_clean and chain {ref_l} and resi 89-105 and name CA") < 3:
                    notes += "L3_unresolved_in_ref; "
                else:
                    try:
                        cmd.super(f"pred_clean and chain {pred_l} and resi 1-120 and ss s", 
                                  f"ref_clean and chain {ref_l} and resi 1-120 and ss s")
                        l3_result = cmd.align(
                            f"pred_clean and chain {pred_l} and resi 89-105 and name CA",
                            f"ref_clean  and chain {ref_l}  and resi 89-105 and name CA",
                            transform=0, cycles=0
                        )
                        calc_l3 = f"{l3_result[0]:.4f}"
                    except Exception as e:
                        calc_l3, notes = "ERR_CRASH", notes + f"l3_crashed:{str(e)[:40]}; "

                # ── Step 5: TM-Score via USalign ─────────────────────────────
                # temp_pred/temp_ref are saved unconditionally (not gated on
                # USalign being configured) since Step 6 below reuses them and
                # should not silently be skipped just because TM-score isn't set up.
                true_tm = "N/A"
                cmd.save(temp_pred, "pred_clean")
                cmd.save(temp_ref,  "ref_clean")
                if os.path.exists(tmalign_exe):
                    res = subprocess.run([tmalign_exe, temp_pred, temp_ref], capture_output=True, text=True)
                    for line in res.stdout.split('\n'):
                        if line.startswith("TM-score="):
                            true_tm = line.split()[1]; break

                # ── Step 6: Ramachandran-outlier % and clashscore (subprocess) ─
                # Reuses the temp_pred file already saved above for TM-align --
                # no extra write needed.
                if os.path.exists(temp_pred):
                    calc_rama, calc_clash, rama_note = _run_rama_clash(temp_pred)
                    if rama_note:
                        notes += rama_note + "; "

                # ── Step 7: Heavy/light interface-contact recovery (subprocess) ─
                # Needs H/L-relabeled copies (Arpeggio worker expects chains
                # literally named "H" and "L"), so build small temp files
                # containing just the identified Fv chains under those names.
                if ref_h and pred_h:
                    temp_ref_hl = os.path.join(os.path.dirname(temp_ref), "temp_ref_hl.pdb")
                    temp_pred_hl = os.path.join(os.path.dirname(temp_pred), "temp_pred_hl.pdb")

                    def _relabel_and_renumber(sel_obj, source_sel, h_chain, l_chain):
                        """Isolate the H/L chains, relabel to canonical H/L, and
                        renumber each chain to start at 1. Renumbering matters:
                        several models (e.g. AbFold) use continuous numbering
                        across chains rather than restarting per chain (H:0-129,
                        L:130-259), which would otherwise make the Arpeggio
                        worker's own internal Fv-window clipping silently drop
                        an entire chain."""
                        cmd.create(sel_obj, source_sel)
                        # Two-step temp-chain-id swap to avoid collisions if a
                        # native chain ID already happens to be "H" or "L".
                        cmd.alter(f"{sel_obj} and chain {h_chain}", 'chain="ZH"')
                        if l_chain:
                            cmd.alter(f"{sel_obj} and chain {l_chain}", 'chain="ZL"')
                        cmd.alter(f"{sel_obj} and chain ZH", 'chain="H"')
                        if l_chain:
                            cmd.alter(f"{sel_obj} and chain ZL", 'chain="L"')
                        cmd.sort(sel_obj)
                        # Renumber each chain 1..N in file order (high-offset
                        # two-step remap to avoid collisions mid-renumber).
                        for ch in (["H", "L"] if l_chain else ["H"]):
                            model = cmd.get_model(f"{sel_obj} and chain {ch} and name CA")
                            for i, atom in enumerate(model.atom):
                                cmd.alter(f"{sel_obj} and chain {ch} and resi \\{atom.resi}",
                                          f"resi={10000 + i}")
                            cmd.sort(sel_obj)
                        for ch in (["H", "L"] if l_chain else ["H"]):
                            model = cmd.get_model(f"{sel_obj} and chain {ch} and name CA")
                            for i, atom in enumerate(model.atom):
                                cmd.alter(f"{sel_obj} and chain {ch} and resi \\{atom.resi}",
                                          f"resi={i + 1}")
                        cmd.sort(sel_obj)

                    try:
                        _relabel_and_renumber("ref_hl_tmp", "ref_clean", ref_h, ref_l)
                        cmd.save(temp_ref_hl, "ref_hl_tmp")
                        cmd.delete("ref_hl_tmp")

                        _relabel_and_renumber("pred_hl_tmp", "pred_clean", pred_h, pred_l)
                        cmd.save(temp_pred_hl, "pred_hl_tmp")
                        cmd.delete("pred_hl_tmp")

                        calc_ifrec, ifrec_note = _run_interface_recovery(temp_ref_hl, temp_pred_hl)
                        if ifrec_note:
                            notes += ifrec_note + "; "
                    except Exception as e:
                        notes += f"interface_recovery_setup_crashed:{str(e)[:60]}; "
                    finally:
                        for f_tmp in [temp_ref_hl, temp_pred_hl]:
                            if os.path.exists(f_tmp):
                                os.remove(f_tmp)

                print(f"{pdb_id:<6} | {true_tm:<8} | {calc_overall:<10} | {calc_masked:<10} | {calc_zero:<10} | {calc_delta:<8} | {calc_h3:<8} | {calc_l3:<8} | {calc_plddt:<7} | {calc_ptm:<7} | {calc_iptm:<7} | {calc_conf:<7} | {calc_rama:<7} | {calc_clash:<9} | {calc_ifrec:<7}")
                writer.writerow([pdb_id, true_tm, calc_overall, calc_masked, calc_zero, calc_delta, calc_h3, calc_l3, calc_plddt, calc_ptm, calc_iptm, calc_conf, calc_rama, calc_clash, calc_ifrec, status, notes])
                outfile.flush()

            except Exception as e:
                err_msg = str(e)[:60]
                print(f"{pdb_id:<6} | FATAL ERROR: {err_msg}")
                writer.writerow([pdb_id, "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "Failed", err_msg])
                outfile.flush()

    for f in [temp_ref, temp_pred]:
        if os.path.exists(f): os.remove(f)

    cmd.set('suspend_undo', 0)
    print("=" * 160)
    print(f"PIPELINE COMPLETE. Saved to {output_csv}")

cmd.extend("run_benchmarking_pipeline", run_benchmarking_pipeline)