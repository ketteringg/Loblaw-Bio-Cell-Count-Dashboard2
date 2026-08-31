"""
generate_outputs.py

Produces the required Part 2-4 deliverable files by calling the same
analysis.py functions the interactive dashboard is built on. This keeps
these files guaranteed to match what the dashboard shows for the
equivalent fixed cohort, rather than being a separately maintained copy
of the same logic.

The dashboard (app.py) is deliberately general purpose (any cohort, any
comparison) and doesn't have fixed tabs for Parts 2, 3, or 4
specifically, since those are just special cases of what the dashboard
can already do. This script exists so the required, graded outputs
still exist as concrete files in the repo, without needing to click
through the dashboard UI to reconstruct the exact required filter
combination.

Run after load_data.py:
    python generate_outputs.py

Writes to the repository root (file names match the assignment part
each one answers):
    part2_frequency_table.csv           Part 2
    part2_sample_total_counts.csv       Part 2 (per-sample total cell
                                         count step, shown standalone)
    part3_stats_results.csv             Part 3
    part3_boxplot_responders.png        Part 3
    part4_baseline_melanoma_samples.csv Part 4.1
    part4_summary.txt                   Part 4.2
"""
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display available in Codespaces/CI
import matplotlib.pyplot as plt
import numpy as np

from analysis import (
    POPULATIONS,
    get_frequency_table,
    get_sample_totals,
    get_responder_comparison,
    run_stats_test,
    get_baseline_melanoma_samples,
    get_baseline_summary,
)

ROOT = Path(__file__).parent
DB_PATH = ROOT / "cell_counts.db"


def main() -> None:
    if not DB_PATH.exists():
        raise SystemExit(
            f"{DB_PATH} not found. Run `python load_data.py` first to build it."
        )

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        # --- Part 2: frequency table ---
        freq = get_frequency_table(conn)
        freq.to_csv(ROOT / "part2_frequency_table.csv", index=False)
        print(f"Wrote part2_frequency_table.csv ({len(freq)} rows)")

        # Part 2's specific intermediate step, spelled out as its own
        # file: "For each sample, calculate the total number of cells by
        # summing the counts across all five populations." One row per
        # sample rather than one row per (sample, population), since
        # that total already appears repeated 5 times inside
        # part2_frequency_table.csv and this makes the standalone step
        # easy to check directly.
        totals = get_sample_totals(conn)
        totals.to_csv(ROOT / "part2_sample_total_counts.csv", index=False)
        print(f"Wrote part2_sample_total_counts.csv ({len(totals)} rows)")

        # --- Part 3: stats + boxplot ---
        # get_responder_comparison() already restricts to melanoma,
        # miraclib treated, PBMC samples, matching the assignment's
        # wording exactly ("melanoma patients receiving miraclib...
        # Please only include PBMC samples").
        comparison = get_responder_comparison(conn)
        stats = run_stats_test(comparison)
        stats.to_csv(ROOT / "part3_stats_results.csv", index=False)
        print(f"Wrote part3_stats_results.csv ({len(stats)} rows)")

        # Plain matplotlib, not seaborn: the boxplot/stripplot combination
        # doesn't need an extra dependency for something matplotlib
        # already does natively, and keeping requirements.txt minimal
        # means fewer ways `make setup` can fail in a fresh Codespaces
        # environment.
        fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(3 * len(POPULATIONS), 5), sharey=True)
        rng = np.random.default_rng(0)  # fixed seed, so output is reproducible run to run
        for ax, pop in zip(axes, POPULATIONS):
            pop_df = comparison[comparison["population"] == pop]
            responders = pop_df[pop_df["response"] == "yes"]["percentage"]
            non_responders = pop_df[pop_df["response"] == "no"]["percentage"]
            ax.boxplot([responders, non_responders], tick_labels=["Responder", "Non-responder"])
            # Individual points, jittered slightly on x so they don't all
            # overlap in a single vertical line: the matplotlib
            # equivalent of what sns.stripplot was doing.
            for i, values in enumerate([responders, non_responders], start=1):
                jitter = rng.uniform(-0.08, 0.08, size=len(values))
                ax.scatter(i + jitter, values, color="black", alpha=0.25, s=6, zorder=3)
            p = stats.loc[stats["population"] == pop, "p_value"].values[0]
            ax.set_title(f"{pop}\np={p:.4f}", fontsize=10)
            ax.tick_params(axis="x", rotation=20)
        axes[0].set_ylabel("% of total cells")
        fig.suptitle("Responder vs non-responder comparison (melanoma, miraclib, PBMC)")
        fig.tight_layout()
        fig.savefig(ROOT / "part3_boxplot_responders.png", dpi=150)
        plt.close(fig)
        print("Wrote part3_boxplot_responders.png")

        sig = stats[stats["significant_bonferroni"]]["population"].tolist()
        print(f"  Conclusion: {', '.join(sig)} show a significant difference "
              f"(Bonferroni corrected p < 0.05).")

        # --- Part 4: baseline cohort ---
        baseline = get_baseline_melanoma_samples(conn)
        baseline.to_csv(ROOT / "part4_baseline_melanoma_samples.csv", index=False)
        print(f"Wrote part4_baseline_melanoma_samples.csv ({len(baseline)} rows)")

        summary = get_baseline_summary(conn)
        with open(ROOT / "part4_summary.txt", "w") as f:
            f.write("Part 4: Melanoma PBMC baseline (time=0), miraclib treated\n")
            f.write("=" * 60 + "\n")
            f.write(f"Total samples: {len(baseline)}\n\n")
            f.write("Samples per project:\n")
            f.write(summary["samples_per_project"].to_string(index=False) + "\n\n")
            f.write("Subjects by response:\n")
            f.write(summary["subjects_by_response"].to_string(index=False) + "\n\n")
            f.write("Subjects by sex:\n")
            f.write(summary["subjects_by_sex"].to_string(index=False) + "\n")
        print("Wrote part4_summary.txt")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
