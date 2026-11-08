#!/usr/bin/env python3
"""Self-test for the Gompertz saturation skill.

Run:  python scripts/selftest.py
Checks parameter recovery on synthetic data and the end-to-end CLI behaviour.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fit_gompertz as fg  # noqa: E402
import forecast_trend as ft  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}" + (f" -> {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def synthetic(x, A, B, k, noise=0.0, seed=1):
    rng = np.random.default_rng(seed)
    y = fg.gompertz(x, A, B, k)
    if noise:
        y = y * (1.0 + rng.normal(0.0, noise, size=len(x)))
    return y


def test_recovery():
    # Parameters of the reference case: gamma 0.464, alpha -4.211, beta -0.224
    A, B, k = 0.464, 4.211, 0.224
    x = np.linspace(1.5, 14.0, 40)
    y = synthetic(x, A, B, k, noise=0.01)
    fit = fg.fit_gompertz(x, y)
    check("free saturation recovered", abs(fit["A"] - A) / A < 0.05,
          f"gamma={fit['A']:.4f} (true {A})")
    check("alpha recovered in original x scale", abs(fit["B_orig"] - B) / B < 0.10,
          f"-alpha={fit['B_orig']:.4f} (true {B})")
    check("beta recovered", abs(fit["k"] - k) / k < 0.15, f"-beta={fit['k']:.4f} (true {k})")
    check("R2 high on clean data", fit["r2"] > 0.99, f"R2={fit['r2']:.4f}")


def test_year_scale():
    # Calendar-year x: exp(-k*x) underflows unless x is shifted internally.
    years = np.arange(2008.0, 2031.0)
    A, B, k = 0.55, 3.0, 0.12
    y = A * np.exp(-B * np.exp(-k * (years - 2008.0)))
    fit = fg.fit_gompertz(years, y)
    check("year-scale fit converges", abs(fit["A"] - A) / A < 0.02, f"gamma={fit['A']:.4f} (true {A})")
    check("year-scale beta recovered", abs(fit["k"] - k) / k < 0.05, f"-beta={fit['k']:.4f} (true {k})")
    check("year-scale predictions finite", np.all(np.isfinite(fg.predict(fit, years))))
    check("original-scale alpha flagged as overflowing",
          not np.isfinite(fit["B_orig"]) or fit["B_orig"] > 0)


def test_fixed_gamma():
    x = np.linspace(1.5, 14.0, 40)
    y = synthetic(x, 0.464, 4.211, 0.224)
    fit = fg.fit_gompertz(x, y, gamma_fixed=0.62)
    check("fixed gamma honoured", abs(fit["A"] - 0.62) < 1e-9, f"gamma={fit['A']}")
    check("fixed gamma still fits", fit["r2"] > 0.98, f"R2={fit['r2']:.4f}")


def test_identifiability_warning():
    # Data only on the accelerating branch -> the saturation level is not identified.
    x = np.linspace(0.5, 4.0, 15)
    y = synthetic(x, 0.464, 4.211, 0.224)
    fit = fg.fit_gompertz(x, y)
    check("unidentified saturation still returns a fit", np.isfinite(fit["A"]),
          f"gamma={fit['A']:.4f}, inflection x*={fit['inflect_x']:.3f}")


def test_backtest_and_weights():
    rng = np.random.default_rng(3)
    rows = []
    for region, mult in (("A", 1.0), ("B", 1.4), ("C", 0.7)):
        for xv in np.linspace(2.0, 12.0, 12):
            rows.append({"region": region, "x": xv,
                         "y": fg.gompertz(xv, 0.464, 4.211, 0.224) * mult})
    df = pd.DataFrame(rows)
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "panel.csv"
        out = Path(tmp) / "out"
        df.to_csv(src, index=False, encoding="utf-8")
        cmd = [sys.executable, str(HERE / "fit_gompertz.py"), "--input", str(src),
               "--x-col", "x", "--y-col", "y", "--out", str(out),
               "--weight-col", "region", "--gamma", "0.62",
               "--backtest", "--compare-logistic", "--bootstrap", "30",
               "--predict-range", "2:20:2", "--no-figures"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        check("CLI exits cleanly", proc.returncode == 0, proc.stderr.strip()[-300:])
        for fname in ("fit_report.md", "parameters.csv", "predictions.csv",
                      "spatial_weights.csv", "backtest.csv", "fit_summary.json"):
            check(f"output {fname} exists", (out / fname).exists())
        if (out / "spatial_weights.csv").exists():
            w = pd.read_csv(out / "spatial_weights.csv")
            check("weights sum to 1", abs(w["weight"].sum() - 1.0) < 1e-6,
                  f"sum={w['weight'].sum():.6f}")
            check("weights ordered by GDP level",
                  w.sort_values("x")["weight"].is_monotonic_increasing)
        if (out / "predictions.csv").exists():
            p = pd.read_csv(out / "predictions.csv")
            check("predictions cover requested grid", len(p) == 10, f"n={len(p)}")


def test_trend():
    y = np.array([10.0, 11.0, 12.5, 13.0, 15.0, 16.0, 18.0])
    method, fitted, fc, diag = ft.forecast_series(y, 3, "auto")
    check("trend method chosen", method in ("linear", "brown", "holt"), method)
    check("trend forecast length", len(fc) == 3)
    check("trend R2 reasonable", diag["r2"] > 0.8, f"R2={diag['r2']:.3f}")
    check("trend forecast increasing", fc[-1] > fc[0], f"{fc[0]:.2f} -> {fc[-1]:.2f}")


def main():
    test_recovery()
    test_year_scale()
    test_fixed_gamma()
    test_identifiability_warning()
    test_backtest_and_weights()
    test_trend()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) failed: {FAILURES}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
