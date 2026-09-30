"""Compare TNP's surface patch metrics across implementations.

Calculates the PSH, PPC and PNC metrics for the clinical-stage nanobody models in
`paper/paper_data/vhh_clinical_set_models`, and compares them with the values
published in `paper/paper_data/insilico_descriptors`.  Optionally, also compares
them with the metrics from another TNP source tree, such as a checkout from
before `psa` was replaced with FreeSASA:

    git worktree add /tmp/tnp-psa <commit>
    python scripts/compare_sasa.py --baseline /tmp/tnp-psa

This is a development aid and is not part of the package.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parent.parent
MODELS = ROOT / "paper/paper_data/vhh_clinical_set_models"
PUBLISHED = ROOT / "paper/paper_data/insilico_descriptors/VHH_TSD_all_properties_FINAL.csv"
METRICS = {
    "PSH": "Patch_Hydrophob_CDR",
    "PPC": "Patch_Pos_Charge_CDR",
    "PNC": "Patch_Neg_Charge_CDR",
}
PUBLISHED_COLUMNS = {
    "PSH": "Patches CDR Surface Hydrophobicity",
    "PPC": "Patches CDR Positive Charge",
    "PNC": "Patches CDR Negative Charge",
}


def compute():
    """Calculate the metrics for each model with the importable TNP."""
    from theraprofnano.Hydrophobicity_and_Charge_Profiler.Hydrophobicity_and_Charge_Assigner import (  # noqa: E501
        CreateAnnotation,
    )

    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        for model in sorted(MODELS.glob("*.pdb")):
            # CreateAnnotation writes an annotated copy beside its input.
            path = shutil.copy(model, tmp)
            stats = CreateAnnotation(0, 7.4, path, "IG", "imgt", verbose=False)[0]
            results[model.stem] = {m: stats[key] for m, key in METRICS.items()}
    return results


def baseline(source_tree):
    """Calculate the metrics with the TNP in another source tree."""
    env = dict(os.environ, PYTHONPATH=str(Path(source_tree).resolve() / "src"))
    output = subprocess.run(
        [sys.executable, __file__, "--compute-only"],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return json.loads(output.splitlines()[-1])


def published():
    """Return the published metrics, keyed by model name via sequence."""
    parser = PDBParser(QUIET=True)
    by_sequence = {}
    for model in MODELS.glob("*.pdb"):
        residues = parser.get_structure(model.stem, model).get_residues()
        by_sequence["".join(seq1(r.get_resname()) for r in residues)] = model.stem

    results = {}
    with open(PUBLISHED) as f:
        for row in csv.DictReader(f):
            model = by_sequence.get(row["Sequence"])
            if model:
                results[model] = {
                    m: float(row[column]) for m, column in PUBLISHED_COLUMNS.items()
                }
            else:
                print(f"No model found for {row['SeqID']}", file=sys.stderr)
    return results


def report(name, reference, new):
    from theraprofnano.cli import assign_flag

    models = sorted(set(reference) & set(new), key=lambda m: int(m[3:].split("_")[0]))
    print(f"\n## FreeSASA vs {name} ({len(models)} models)\n")
    print("| Metric | Pearson r | Mean abs. diff. | Max abs. diff. | Flag changes |")
    print("| --- | --- | --- | --- | --- |")
    changes = []
    for metric in METRICS:
        a = np.array([reference[m][metric] for m in models])
        b = np.array([new[m][metric] for m in models])
        r = np.corrcoef(a, b)[0, 1] if a.std() and b.std() else float("nan")
        flipped = [
            (m, metric, reference[m][metric], new[m][metric])
            for m in models
            if assign_flag(metric, reference[m][metric])
            != assign_flag(metric, new[m][metric])
        ]
        changes += flipped
        diff = np.abs(a - b)
        print(
            f"| {metric} | {r:.4f} | {diff.mean():.4f} | {diff.max():.4f} "
            f"| {len(flipped)} |"
        )
    for m, metric, old, new_value in changes:
        print(
            f"\n- {m} {metric}: {old:.4f} ({assign_flag(metric, old)}) → "
            f"{new_value:.4f} ({assign_flag(metric, new_value)})"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--baseline", help="another TNP source tree to compare")
    parser.add_argument("--compute-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.compute_only:
        print(json.dumps(compute()))
        return

    new = compute()
    report("published values", published(), new)
    if args.baseline:
        report(f"`{args.baseline}`", baseline(args.baseline), new)


if __name__ == "__main__":
    main()
