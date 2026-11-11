"""Summarize generated donor support without exporting recipient records."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


SCENARIOS = ("baseline", "both_lower", "both_higher", "busy_lower", "absence_lower")
MEANS = (
    "unique_categories",
    "candidate_mean",
    "tilted_mean",
    "conditional_mean_shift",
    "agreement_pp_shift",
)


def summarize(input_root: Path, chains: int, iterations: int) -> list[dict]:
    groups = defaultdict(list)
    for scenario in SCENARIOS:
        source = input_root / "scenarios" / scenario / f"diagnostics_m{chains}_i{iterations}" / "staffing_updates.csv"
        with source.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"scenario", "seed", "target", "chain", "iteration", "single_category", *MEANS}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"Required donor-support fields are absent for {scenario}.")
            for row in reader:
                if row["scenario"] != scenario or row["target"] not in ("Q11", "Q12"):
                    raise ValueError(f"Unexpected scenario or staffing target for {scenario}.")
                if row["single_category"] not in ("TRUE", "FALSE"):
                    raise ValueError("Invalid single-category indicator.")
                key = (scenario, int(row["seed"]), row["target"])
                # Retain only computational counts and aggregate ingredients.
                groups[key].append({
                    "chain": int(row["chain"]),
                    "iteration": int(row["iteration"]),
                    "single_category": int(row["single_category"] == "TRUE"),
                    **{field: float(row[field]) for field in MEANS},
                })
    expected = {(scenario, seed, target) for scenario in SCENARIOS
                for seed in (20261111, 20261112) for target in ("Q11", "Q12")}
    if set(groups) != expected:
        raise ValueError("The five scenarios, two seeds, and two staffing targets are required.")
    output = []
    for key in sorted(groups):
        scenario, seed, target = key
        rows = groups[key]
        if ({row["chain"] for row in rows} != set(range(1, chains + 1)) or
                {row["iteration"] for row in rows} != set(range(1, iterations + 1))):
            raise ValueError(f"Incomplete chain or iteration coverage for {key}.")
        expected_per_update = 3 if target == "Q11" else 2
        if len(rows) != chains * iterations * expected_per_update:
            raise ValueError(f"Unexpected number of staffing updates for {key}.")
        multiplicities = Counter((row["chain"], row["iteration"]) for row in rows)
        if (len(multiplicities) != chains * iterations or
                any(count != expected_per_update for count in multiplicities.values())):
            raise ValueError(f"Incomplete staffing-cell coverage within a chain and iteration for {key}.")
        if any(not math.isfinite(row[field]) for row in rows for field in MEANS):
            raise ValueError(f"Nonfinite donor-support value for {key}.")
        output.append({
            "scenario": scenario,
            "seed": seed,
            "target": target,
            "n_updates": len(rows),
            "n_chains": chains,
            "n_iterations": iterations,
            "single_category_fraction": math.fsum(row["single_category"] for row in rows) / len(rows),
            **{f"mean_{field}": math.fsum(row[field] for row in rows) / len(rows) for field in MEANS},
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-file", type=Path, required=True)
    parser.add_argument("--m", type=int, choices=(50, 100, 200), default=50)
    parser.add_argument("--iterations", type=int, choices=(50, 100, 200), default=50)
    args = parser.parse_args()
    if args.output_file.exists():
        raise FileExistsError("Use a new output file; existing summaries are preserved.")
    rows = summarize(args.input_root, args.m, args.iterations)
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    with args.output_file.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} scenario/seed/target donor-support summaries.")


if __name__ == "__main__":
    main()
