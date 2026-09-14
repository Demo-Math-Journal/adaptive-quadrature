"""
quadrature.py

Numerical integration via a fixed-panel composite Simpson's rule and an
adaptive Simpson's rule that recursively subdivides an interval until a
local error estimate falls below a tolerance. Adaptive quadrature
concentrates panels where the integrand is changing quickly and uses
far fewer function evaluations than a fixed grid for the same accuracy
on integrands with localized features.
"""

from __future__ import annotations
from typing import Callable, List, Tuple


def simpson(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    """Composite Simpson's rule on n subintervals (n must be even)."""
    if n % 2 != 0:
        raise ValueError("n must be even for composite Simpson's rule")
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        x = a + i * h
        total += (4 if i % 2 else 2) * f(x)
    return total * h / 3.0


def _simpson_panel(f: Callable[[float], float], a: float, b: float) -> Tuple[float, float]:
    """Return (Simpson estimate, midpoint) for one panel [a, b]."""
    m = 0.5 * (a + b)
    s = (b - a) / 6.0 * (f(a) + 4 * f(m) + f(b))
    return s, m


def adaptive_simpson(
    f: Callable[[float], float],
    a: float,
    b: float,
    tol: float = 1e-8,
    max_depth: int = 50,
) -> Tuple[float, List[Tuple[float, float]]]:
    """Adaptive Simpson's rule. Returns (integral estimate, panel list),
    where panel list is the list of (left, right) endpoints of every leaf
    panel actually used, useful for visualizing where the method refined."""
    panels: List[Tuple[float, float]] = []

    def recurse(a: float, b: float, whole: float, depth: int) -> float:
        m = 0.5 * (a + b)
        left, _ = _simpson_panel(f, a, m)
        right, _ = _simpson_panel(f, m, b)
        if depth <= 0 or abs(left + right - whole) < 15 * tol:
            panels.append((a, b))
            return left + right + (left + right - whole) / 15.0
        return recurse(a, m, left, depth - 1) + recurse(m, b, right, depth - 1)

    whole, _ = _simpson_panel(f, a, b)
    total = recurse(a, b, whole, max_depth)
    return total, panels


if __name__ == "__main__":
    import math

    # Integrand with a sharp peak: adaptive quadrature should concentrate
    # panels near x = 0.3 while the fixed rule spreads them uniformly.
    f = lambda x: 1.0 / ((x - 0.3) ** 2 + 1e-4)

    fixed = simpson(f, -10.0, 10.0, n=2000)
    adaptive, panels = adaptive_simpson(f, -10.0, 10.0, tol=1e-6, max_depth=60)

    print(f"Fixed Simpson (n=1000):     {fixed:.6f}")
    print(f"Adaptive Simpson (tol=1e-6): {adaptive:.6f} using {len(panels)} panels")
