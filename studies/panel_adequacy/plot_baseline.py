#!/usr/bin/env python3
"""Plot frozen results; the truth-oracle bound is a post-hoc interpretation aid."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=Path, default=Path("results/panel_adequacy_baseline"))
    args = p.parse_args()
    summary = json.loads((args.results / "summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), sharex=True)
    styles = [
        ("quality", "Genotype quality", "#8e6c8a", "-"),
        (summary["selected_native"], "Unexplained-read fraction (selected)", "#007a87", "-"),
        ("logistic", "Combined logistic baseline", "#ce671b", "-"),
        ("unexpl_reads", "Unexplained-read count", "#3761b0", "--"),
    ]
    for ax, condition, title in zip(axes, ["loo", "full"],
                                   ["Leave-one-out panel", "Matched full-panel control"]):
        m = summary["test"][condition]["metrics_by_divergence_threshold"]["0.005"]
        for method, label, color, style in styles:
            points = [m[f"matched_{r}:{method}"] for r in [.8, .9, .95, .99, 1.0]]
            ax.plot([100*x["coverage"] for x in points], [x["accepted_wrong"] for x in points],
                    label=label, color=color, linestyle=style, marker="o", markersize=4, linewidth=1.8)
        n = int(m["cutoff:no_filter"]["eligible"])
        errors = int(m["cutoff:no_filter"]["accepted_wrong"])
        ks = list(range(int(.8*n), n+1))
        ax.plot([100*k/n for k in ks], [max(0, errors-(n-k)) for k in ks],
                color="#222222", linestyle=":", linewidth=1.5, label="Truth-oracle lower bound")
        ax.axvline(100*766/806, color="#bbbbbb", linewidth=.8, zorder=0)
        ax.set(title=title, xlabel="Scorable calls retained (%)", ylabel="Retained discrepant calls")
        ax.set_xlim(79.5, 100.5)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#e8e8e8", linewidth=.7)
        ax.spines[["top", "right"]].set_visible(False)
        if condition == "loo":
            for method, dy in [(summary["selected_native"], -4), ("logistic", 5), ("unexpl_reads", -6)]:
                val = m[f"matched_0.95:{method}"]["accepted_wrong"]
                ax.annotate(str(int(val)), (100*766/806, val), xytext=(7,dy),
                            textcoords="offset points", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5,.01), ncol=3, frameon=False, fontsize=9)
    fig.suptitle("Ordinary confidence filters leave residual sequence discrepancies", fontsize=13, x=.5, y=.99)
    fig.text(.5,.905,"Discrepancy = >0.5% sequence difference; 806 held-out cases; historical calls; matched retention",
             ha="center", fontsize=9, color="#555555")
    fig.subplots_adjust(bottom=.27, top=.83, wspace=.32, left=.07, right=.98)
    fig.savefig(args.results / "retention_curves.png", dpi=180)
    fig.savefig(args.results / "retention_curves.svg")


if __name__ == "__main__":
    main()
