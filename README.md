# Adaptive Quadrature

Error-controlled numerical integration: adaptive Simpson's rule versus a fixed-grid Simpson's rule, compared at matched function-evaluation budgets.

A fixed grid has to be fine everywhere to resolve a sharp, localized feature, even where the integrand is flat and easy. Adaptive quadrature instead recurses only where a local error estimate says it needs to, concentrating function evaluations near the hard part of the integrand. The experiment here integrates a narrow peak sitting inside a wide, mostly-flat domain — the regime where that difference actually matters — and shows adaptive Simpson's rule reaching several orders of magnitude lower error than a fixed grid given the same number of evaluations.

**Topics:** numerical-analysis, numerical-integration, quadrature, adaptive-methods, math

## Structure

- `README.md` — this file
- `src/` — composite Simpson's rule and adaptive Simpson's rule (recursive, tolerance-driven)
- `experiments/` — matched-budget comparison of adaptive vs. fixed quadrature on a narrow peak in a wide domain, across four tolerances
- `app/` — placeholder for an interactive browser demo (visualize where adaptive quadrature places its panels)
- `pdf/` — placeholder for the write-up
