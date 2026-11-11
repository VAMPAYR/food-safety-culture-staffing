"""Draw sensitivity figures from saved aggregate outputs; never fit a model.

Example, from the reproducibility directory:
    python manuscript/make_priority123_figures.py

Both directory arguments may be replaced for a public reproduction run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
OUTCOMES = ("Q6", "Q4")
METHODS = ("site", "CV3")
STEM = "figure_A1_jurisdiction_variance_shares"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def truth(value: str) -> bool:
    return value.strip().upper() == "TRUE"


def figure_data(input_dir: Path) -> tuple[list[dict], dict]:
    paths = {
        name: input_dir / name
        for name in (
            "joint_covariance_summary.csv",
            "joint_covariance_geometry.csv",
            "jurisdiction_covariance_contributions.csv",
        )
    }
    summary = read_csv(paths["joint_covariance_summary.csv"])
    geometry = read_csv(paths["joint_covariance_geometry.csv"])
    contributions = read_csv(paths["jurisdiction_covariance_contributions.csv"])
    if len(summary) != 4 or len(geometry) != 8 or len(contributions) != 32:
        raise ValueError("Expected four covariances, eight directions, and 32 jurisdiction rows.")

    rows = []
    group_checks = []
    for outcome in OUTCOMES:
        outcome_n = {int(r["n"]) for r in summary if r["outcome"] == outcome}
        if len(outcome_n) != 1:
            raise ValueError(f"Covariance methods have different worker counts for {outcome}.")
        for method in METHODS:
            group = [r for r in contributions if (r["outcome"], r["covariance_method"]) == (outcome, method)]
            covariance = [r for r in summary if (r["outcome"], r["covariance_method"]) == (outcome, method)]
            direction = [r for r in geometry if (r["outcome"], r["covariance_method"], r["eigen_index"]) == (outcome, method, "1")]
            if len(covariance) != 1 or len(direction) != 1:
                raise ValueError(f"Missing or duplicated covariance/direction for {outcome}/{method}.")
            if sorted(int(r["site"]) for r in group) != list(range(1, 9)):
                raise ValueError(f"Jurisdictions are missing or duplicated for {outcome}/{method}.")
            if not truth(covariance[0]["valid"]) or not truth(direction[0]["valid"]) or not truth(direction[0]["orientation_unique"]):
                raise ValueError(f"An accepted unique covariance direction is required for {outcome}/{method}.")
            for row in group:
                if not truth(row["valid"]) or not truth(row["orientation_unique"]):
                    raise ValueError(f"An invalid contribution cannot be displayed: {outcome}/{method}/{row['site']}.")
            values = [float(r["small_eigen_direction_share"]) for r in group]
            if any(not math.isfinite(v) or v < 0 or v > 1 for v in values):
                raise ValueError(f"Invalid variance share for {outcome}/{method}.")
            if abs(sum(values) - 1) > 1e-8:
                raise ValueError(f"Variance shares do not sum to one for {outcome}/{method}.")
            if max(values) * 100 > 60:
                raise ValueError("A share exceeds the prespecified 0-60% figure axis; revise the display explicitly.")
            group_checks.append({"outcome": outcome, "covariance_method": method, "share_sum": sum(values), "n_jurisdictions": len(group)})
            for row in sorted(group, key=lambda r: int(r["site"])):
                share = float(row["small_eigen_direction_share"])
                rows.append({
                    "outcome": outcome,
                    "covariance_method": method,
                    "jurisdiction": int(row["site"]),
                    "variance_share": share,
                    "variance_share_percent": 100 * share,
                    "n_workers": int(covariance[0]["n"]),
                    "direction_busy_coverage": float(direction[0]["direction_B"]),
                    "direction_absence_coverage": float(direction[0]["direction_A"]),
                })
    provenance = {
        "source_sha256": {name: sha256(path) for name, path in paths.items()},
        "group_checks": group_checks,
        "n_plotted_rows": len(rows),
        "source": "Saved priority 1 aggregate output; no model estimation.",
        "interpretation": "Each share is the jurisdiction's contribution to estimated covariance along that outcome and method's own smallest-eigenvalue direction. It is not a share of the Wald statistic, a calibrated influence statistic, or an effective degree of freedom.",
        "axis_percent": [0, 60],
        "matplotlib_version": matplotlib.__version__,
    }
    return rows, provenance


def draw(rows: list[dict], output_dir: Path) -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })
    colors = {"site": "#245A81", "CV3": "#CA8C34"}
    offsets = {"site": -0.16, "CV3": 0.16}
    labels = {"site": "Site covariance", "CV3": "Cluster jackknife (CV3)"}
    outcomes_text = {"Q6": "Perceived sick-worker exclusion", "Q4": "Perceived handwashing"}
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.7), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.85, bottom=0.20, wspace=0.17)
    for ax, outcome in zip(axes, OUTCOMES):
        for method in METHODS:
            group = [r for r in rows if r["outcome"] == outcome and r["covariance_method"] == method]
            y = [r["jurisdiction"] + offsets[method] for r in group]
            x = [r["variance_share_percent"] for r in group]
            ax.barh(y, x, height=0.27, color=colors[method], edgecolor="none", zorder=3)
            for xx, yy in zip(x, y):
                label = "<0.1" if 0 < xx < 0.1 else f"{xx:.1f}"
                ax.annotate(label, (xx, yy), xytext=(3, 0), textcoords="offset points", va="center", ha="left", fontsize=8.5, color="#111111", clip_on=False)
        ax.set_xlim(0, 60)
        ax.set_ylim(8.6, 0.4)
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_yticks(list(range(1, 9)))
        n_workers = next(r["n_workers"] for r in rows if r["outcome"] == outcome)
        title = f"{outcomes_text[outcome]}\n{outcome}, n = {n_workers}"
        ax.set_title(title, loc="left", pad=10, fontweight="bold", linespacing=1.25)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#E4E4E4", linewidth=0.6)
        ax.grid(axis="y", visible=False)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#555555")
        ax.spines["bottom"].set_linewidth(0.6)
        ax.tick_params(axis="y", length=0, pad=6)
        ax.tick_params(axis="x", length=3, width=0.6)
    axes[0].set_ylabel("Jurisdiction", labelpad=9)
    fig.text(0.54, 0.107, "Share of variance in the narrowest direction (%)", ha="center", va="center", fontsize=10)
    fig.legend(handles=[Patch(facecolor=colors[m], edgecolor="none", label=labels[m]) for m in METHODS], loc="lower center", bbox_to_anchor=(0.54, 0.015), ncol=2, frameon=False, fontsize=9.5, handlelength=1.2, columnspacing=1.8)
    metadata = {"Title": "Jurisdiction contributions to staffing covariance", "Creator": "make_priority123_figures.py"}
    fig.savefig(output_dir / f"{STEM}.png", dpi=300, facecolor="white")
    fig.savefig(output_dir / f"{STEM}.pdf", facecolor="white", metadata=metadata)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=ROOT / "outputs" / "public_priority1_20261111_v1")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "priority123_publication_figures")
    parser.add_argument("--overwrite", action="store_true", help="Replace this script's figure, aggregate data, and validation outputs.")
    args = parser.parse_args()
    expected = [args.output_dir / f"{STEM}.{suffix}" for suffix in ("png", "pdf")]
    expected += [args.output_dir / f"{STEM}_data.csv", args.output_dir / "figure_A1_validation.json"]
    if not args.overwrite and any(path.exists() for path in expected):
        raise FileExistsError("Figure outputs already exist; use a new directory or the explicit --overwrite option.")
    rows, provenance = figure_data(args.input_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with expected[2].open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    draw(rows, args.output_dir)
    provenance["output_sha256"] = {path.name: sha256(path) for path in expected[:3]}
    provenance["script_sha256"] = sha256(Path(__file__))
    expected[3].write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} verified aggregate rows and Figure A1 to {args.output_dir}")


if __name__ == "__main__":
    main()
