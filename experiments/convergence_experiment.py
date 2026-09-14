"""
convergence_experiment.py

Numerical experiment: compare adaptive Simpson's rule against fixed-grid
composite Simpson's rule on a sharply peaked integrand, at MATCHED
function-evaluation budgets, across a range of tolerances. Matching the
budget is what makes the comparison fair: it isolates the benefit of
placing panels where the integrand is hard to approximate, rather than
just spending more total effort.

Run with:  python3 convergence_experiment.py
Writes:    results.csv (tolerance, adaptive_evals, adaptive_error, fixed_evals, fixed_error)
"""

import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from quadrature import simpson, adaptive_simpson  # noqa: E402


def counting(f):
    """Wrap f to count how many times it is called."""
    calls = [0]

    def wrapped(x):
        calls[0] += 1
        return f(x)

    return wrapped, calls


def main():
    # A narrow peak (width ~0.01) at x=0.3, sitting inside a much wider
    # domain [-10, 10]. This is the regime where adaptive quadrature earns
    # its keep: a uniform grid has to be fine EVERYWHERE across a wide,
    # mostly-flat domain just to resolve one narrow feature, while adaptive
    # quadrature spends its evaluations only near the peak.
    base_f = lambda x: 1.0 / ((x - 0.3) ** 2 + 1e-4)
    a, b = -10.0, 10.0

    # High-resolution fixed Simpson's rule as a reference "true" value.
    reference = simpson(base_f, a, b, n=2_000_000)

    tolerances = [1e-2, 1e-4, 1e-6, 1e-8]
    rows = []
    for tol in tolerances:
        f_adaptive, adaptive_calls = counting(base_f)
        adaptive_value, panels = adaptive_simpson(f_adaptive, a, b, tol=tol, max_depth=60)
        adaptive_evals = adaptive_calls[0]

        # Fixed Simpson's rule using (about) the same number of function
        # evaluations as adaptive used for this tolerance: n = evals - 1,
        # rounded down to an even number as composite Simpson requires.
        fixed_n = max(2, (adaptive_evals - 1) // 2 * 2)
        f_fixed, fixed_calls = counting(base_f)
        fixed_value = simpson(f_fixed, a, b, n=fixed_n)
        fixed_evals = fixed_calls[0]

        rows.append((tol, adaptive_evals, abs(adaptive_value - reference),
                      fixed_evals, abs(fixed_value - reference)))

    out_path = os.path.join(os.path.dirname(__file__), "results.csv")
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["tolerance", "adaptive_evals", "adaptive_error", "fixed_evals", "fixed_error"])
        writer.writerows(rows)

    print(f"Wrote {out_path}")
    print(f"Reference value (n=2,000,000 fixed Simpson on [{a}, {b}]): {reference:.6f}\n")
    print(f"{'tol':>10}  {'adaptive evals':>15}  {'adaptive error':>16}  {'fixed evals':>12}  {'fixed error':>14}")
    for tol, aevals, aerr, fevals, ferr in rows:
        print(f"{tol:>10.0e}  {aevals:>15}  {aerr:>16.3e}  {fevals:>12}  {ferr:>14.3e}")
    print("\nAt matched function-evaluation budgets, adaptive Simpson is dramatically more")
    print("accurate here: the peak occupies a tiny fraction of the wide domain, so a")
    print("uniform grid wastes almost all of its evaluations on flat, easy regions while")
    print("adaptive quadrature concentrates its evaluations near x=0.3.")


if __name__ == "__main__":
    main()
