#!/usr/bin/env python
"""
rama_clash_worker.py -- standalone Ramachandran-outlier and clashscore computation.

Runs in a SEPARATE Python environment from PyMOL (see README_metrics_setup.md),
using cctbx/mmtbx -- the same MolProbity-equivalent implementation used to
compute this project's published Table 3/4 values (validated: reproduces
AbFold/1R24's clashscore of 418.9 exactly, the value cited in the manuscript).
Invoked via subprocess from benchmark_master_pipeline.py, not imported
directly into the PyMOL process (cctbx and PyMOL's own dependency stack do
not reliably coexist in one interpreter).

Requires:
  1. conda install -c conda-forge cctbx-base   (into THIS script's own env,
     not PyMOL's)
  2. The CCP4/geostd monomer restraint library, needed by the modern
     probe2-based clashscore2 (no external "Probe" binary required, unlike
     the legacy mmtbx.validation.clashscore class):
       git clone --depth 1 https://github.com/phenix-project/geostd.git
     Then set the environment variable MMTBX_CCP4_MONOMER_LIB to point at
     that cloned directory before running this script (see
     README_metrics_setup.md for the exact command).

Usage:
    python rama_clash_worker.py <pdb_or_cif_path>

Prints one line of JSON to stdout:
    {"pct_rama_outliers": 0.68, "clashscore": 15.3, "n_residues": 213}
or, on failure:
    {"error": "..."}
"""
import sys
import json
import os
import io
import contextlib


def main():
    pdb_path = sys.argv[1]
    # mmtbx prints informational hydrogen-placement/bonding chatter directly to
    # stdout from inside the library (not via logging), which would otherwise
    # get mixed in before our JSON result line. Redirect it to a buffer so the
    # ONLY thing this process ever writes to real stdout is the final JSON.
    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured):
            if not os.environ.get("MMTBX_CCP4_MONOMER_LIB") and not os.environ.get("CLIBD_MON"):
                raise RuntimeError(
                    "MMTBX_CCP4_MONOMER_LIB (or CLIBD_MON) is not set -- point it at a "
                    "local clone of https://github.com/phenix-project/geostd before running "
                    "this script. See README_metrics_setup.md."
                )

            import iotbx.phil
            from iotbx.data_manager import DataManager
            from mmtbx.validation.ramalyze import ramalyze
            from mmtbx.validation.clashscore2 import clashscore2
            import mmtbx.probe.Helpers as probeHelpers

            dm = DataManager()
            dm.process_model_file(pdb_path)
            model = dm.get_model()

            rama = ramalyze(model.get_hierarchy(), out=None, quiet=True)
            n_total = rama.n_total
            n_outliers = rama.n_outliers
            pct_outliers = 100.0 * n_outliers / n_total if n_total else 0.0

            master_phil = iotbx.phil.parse(probeHelpers.probe_phil_parameters)
            probe_params = master_phil.fetch().extract().probe
            clash = clashscore2(probe_params, dm, fast=True, condensed_probe=True)
            score = clash.get_clashscore()

        print(json.dumps({
            "pct_rama_outliers": round(pct_outliers, 4),
            "clashscore": round(score, 4),
            "n_residues": n_total,
        }))
    except Exception as e:
        print(json.dumps({"error": str(e)[:300]}))
        sys.exit(1)


if __name__ == "__main__":
    main()
