"""Plot geometry counts and raw peptide/receptor confidence for all 30 states."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior-audit", type=Path, required=True)
    parser.add_argument("--new-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prior = json.loads(args.prior_audit.read_text())
    current = json.loads(args.new_audit.read_text())
    samples = prior["samples"] + current["samples"]
    summaries = {s["state"]: s for s in prior["summary"] + current["summary"]}
    previous = {s["state"] for s in prior["summary"]}
    names = sorted(summaries)
    assert len(names) == 30 and len(samples) == 450
    assert all(s["samples"] == 15 for s in summaries.values())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (counts, confidence) = plt.subplots(1, 2, figsize=(14, 15.5), sharey=True,
                                            gridspec_kw={"width_ratios": [1.0, 1.1], "wspace": .16})
    green, gray, blue = "#128277", "#d7dee3", "#234e70"
    rng = np.random.default_rng(0)
    for i, name in enumerate(names):
        summary = summaries[name]
        counts.barh(i, 15, color=gray, height=.66)
        counts.barh(i, summary["complex_geometry_accepted"], color=green, height=.66)
        counts.text(15.4, i, f'{summary["complex_geometry_accepted"]}/15', va="center", fontsize=10)
        local = sorted((s for s in samples if s["state"] == name), key=lambda s: s["sample"])
        xs = [s["receptor_peptide_iptm"] for s in local]
        confidence.plot([min(xs), max(xs)], [i, i], color="#c5cdd3", lw=1, zorder=1)
        for sample in local:
            confidence.scatter(sample["receptor_peptide_iptm"], i + rng.uniform(-.20, .20),
                               s=25, color=green if sample["accepted"] else "#a5afb8",
                               alpha=.8, edgecolors="white", linewidths=.3, zorder=2)
        if i % 2 == 0:
            for ax in (counts, confidence):
                ax.axhspan(i - .5, i + .5, color="#f6f8f9", zorder=0)
    labels = [name.replace("_S1-", "  S1=").replace("_S3-", ", S3=").replace("_S-", "  S=")
              + (" *" if name in previous else "") for name in names]
    counts.set_yticks(range(len(names)), labels)
    counts.set_ylim(len(names) - .5, -.5)
    counts.set_xlim(0, 17.3)
    counts.set_xticks([0, 5, 10, 15])
    counts.set_xlabel("Samples passing the complex geometry screen")
    counts.set_title("Geometry: all 15 samples retained", loc="left", fontsize=14, pad=16, color=blue)
    confidence.set_xlim(0, 1)
    confidence.set_xticks(np.arange(0, 1.01, .2))
    confidence.set_xlabel("Peptide / TAS2R14 pair ipTM (all raw samples)")
    confidence.set_title("Interface confidence: one dot per sample", loc="left", fontsize=14, pad=16, color=blue)
    for ax in (counts, confidence):
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color="#e7ecef", lw=.7)
        ax.set_axisbelow(True)
    total = sum(s["complex_geometry_accepted"] for s in summaries.values())
    fig.suptitle("HighFold3 | TAS2R14–G protein + flaxseed cyclic peptides", x=.045, y=.982,
                 ha="left", fontsize=20, fontweight="bold", color="#16334a")
    fig.text(.045, .952, f"30 chemical states  ·  450 baseline predictions  ·  {total} pass the geometry screen",
             fontsize=12, color="#4e6170")
    legend = [Line2D([], [], marker="o", lw=0, color=green, label="Pass geometry screen"),
              Line2D([], [], marker="o", lw=0, color="#a5afb8", label="Fail one or more geometry screens")]
    fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(.04, .041), frameon=False, ncol=2)
    fig.text(.045, .027, "* Earlier five-state baseline. Other rows: the 25-state extension. Three seeds × five samples; 10 recycles.",
             fontsize=10, color="#536573")
    fig.text(.045, .013, "Geometry and ipTM do not measure affinity or validate binding poses. The separate 10-sample CLE renumbering trial is excluded.",
             fontsize=10, color="#536573")
    fig.subplots_adjust(left=.225, right=.97, top=.91, bottom=.092)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.with_suffix(".png"), dpi=180, facecolor="white")
    fig.savefig(args.output.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
