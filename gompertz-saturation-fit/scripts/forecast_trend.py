#!/usr/bin/env python3
"""Univariate linear-trend extrapolation (Brown / Holt / OLS) for exogenous drivers.

Use it to project the predictors of a Gompertz model (population, GDP per
capita) before running fit_gompertz.py. Zhang et al. (2020) used the same
idea with SPSS "Brown linear trend" to project Chinese population and GDP per
capita to 2050, then fed the projection into a Gompertz ownership model.

Models
    brown : double exponential smoothing with a linear trend
    holt  : Holt's additive trend, optionally damped
    linear: ordinary least squares trend (baseline)

Selection is by one-step-ahead RMSE on the observed series; --method auto
picks the winner. No statsmodels dependency.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def brown_fit(y, alpha, beta, horizon):
    """Brown double exponential smoothing with a linear trend.

    Returns (fitted_one_step, level, trend, forecast_values).
    """
    n = len(y)
    if n < 3:
        raise ValueError("need at least 3 points")
    level = float(y[0])
    trend = float(y[1] - y[0])
    fitted = np.empty(n)
    for t in range(1, n):
        prev_level = level
        level = alpha * y[t] + (1.0 - alpha) * (level + trend)
        trend = beta * (level - prev_level) + (1.0 - beta) * trend
        fitted[t] = prev_level + trend
    fitted[0] = y[0]
    fc = np.array([level + (m + 1) * trend for m in range(horizon)])
    return fitted, level, trend, fc


def holt_fit(y, alpha, beta, horizon, phi=1.0):
    """Holt additive trend with optional damping."""
    n = len(y)
    level = float(y[0])
    trend = float(y[1] - y[0])
    fitted = np.empty(n)
    for t in range(1, n):
        prev_level = level
        level = alpha * y[t] + (1.0 - alpha) * (level + phi * trend)
        trend = beta * (level - prev_level) + (1.0 - beta) * phi * trend
        fitted[t] = prev_level + trend
    fitted[0] = y[0]
    fc = np.array([level + sum(phi ** (i + 1) for i in range(m + 1)) * trend for m in range(horizon)])
    return fitted, level, trend, fc


def linear_fit(y, horizon):
    n = len(y)
    t = np.arange(n, dtype=float)
    coef = np.polyfit(t, y, 1)
    fitted = np.polyval(coef, t)
    fc = np.polyval(coef, np.arange(n, n + horizon, dtype=float))
    return fitted, fc


def _best_smooth(y, horizon, kind):
    """Grid-search the smoothing parameters; return the best (fitted, forecast)."""
    best = None
    for alpha in np.arange(0.05, 1.0, 0.05):
        for beta in np.arange(0.05, 1.0, 0.05):
            try:
                if kind == "brown":
                    fitted, _, _, fc = brown_fit(y, alpha, beta, horizon)
                else:
                    fitted, _, _, fc = holt_fit(y, alpha, beta, horizon)
            except ValueError:
                continue
            err = float(np.sum((y[1:] - fitted[1:]) ** 2))
            if best is None or err < best[0]:
                best = (err, fitted, fc, alpha, beta)
    if best is None:
        raise ValueError(f"{kind} smoothing failed")
    return best[1], best[2], best[3], best[4]


def score(model, y, horizon):
    """One-step-ahead RMSE/R2 for the fitted values."""
    if model == "linear":
        fitted, _ = linear_fit(y, horizon)
    else:
        fitted, _, _, _ = _best_smooth(y, horizon, model)
    resid = y - fitted
    sse = float(np.sum(resid[1:] ** 2))
    rmse = float(np.sqrt(sse / max(len(y) - 1, 1)))
    mape = float(np.mean(np.abs(resid[1:] / y[1:])) * 100.0)
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum(resid ** 2)) / sst if sst > 0 else float("nan")
    return fitted, rmse, mape, r2


def forecast_series(y, horizon, method="auto"):
    """Return (method_used, fitted, forecast, diagnostics_dict)."""
    cands = ["linear", "brown", "holt"] if method == "auto" else [method]
    results = []
    for name in cands:
        try:
            fitted, rmse, mape, r2 = score(name, y, horizon)
        except Exception:
            continue
        if name == "linear":
            _, fc = linear_fit(y, horizon)
        elif name == "brown":
            _, fc, _, _ = _best_smooth(y, horizon, "brown")
        else:
            _, fc, _, _ = _best_smooth(y, horizon, "holt")
        results.append({"method": name, "fitted": fitted, "forecast": fc, "rmse": rmse, "mape": mape, "r2": r2})
    if not results:
        raise ValueError("no trend model could be fitted")
    best = min(results, key=lambda r: r["rmse"])
    diag = {"chosen": best["method"], "rmse": best["rmse"], "mape": best["mape"], "r2": best["r2"],
            "candidates": [{"method": r["method"], "rmse": r["rmse"], "mape": r["mape"], "r2": r["r2"]}
                           for r in results]}
    return best["method"], best["fitted"], best["forecast"], diag


def load_table(path):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p)
    if p.suffix.lower() == ".csv":
        for enc in ("utf-8-sig", "utf-8", "gb18030", "cp1252", "latin1"):
            try:
                return pd.read_csv(p, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise SystemExit("failed to decode CSV")
    raise SystemExit(f"unsupported input type: {p.suffix}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Brown/Holt linear-trend extrapolation of a yearly series.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--year-col", required=True)
    ap.add_argument("--value-col", action="append", required=True,
                    help="series to extrapolate; repeatable")
    ap.add_argument("--horizon", type=int, default=10, help="number of steps to extrapolate")
    ap.add_argument("--method", choices=["auto", "brown", "holt", "linear"], default="auto")
    ap.add_argument("--out", default="trend_out")
    args = ap.parse_args(argv)

    df = load_table(args.input)
    for col in (args.year_col,) + tuple(args.value_col):
        if col not in df.columns:
            raise SystemExit(f"column '{col}' not found. Available: {list(df.columns)}")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {}
    report = ["# Linear-trend extrapolation", ""]

    for vcol in args.value_col:
        sub = df[[args.year_col, vcol]].dropna().sort_values(args.year_col)
        years = sub[args.year_col].astype(float).to_numpy()
        y = sub[vcol].astype(float).to_numpy()
        if len(y) < 4:
            report += [f"## {vcol}", "", "Skipped: fewer than 4 usable points.", ""]
            continue
        method, fitted, fc, diag = forecast_series(y, args.horizon, args.method)
        fyears = years[-1] + np.arange(1, args.horizon + 1, dtype=float)
        rows = []
        for yr, val, fv in zip(years, y, fitted):
            rows.append({"year": yr, "observed": val, "fitted": fv, "forecast": None})
        for yr, fv in zip(fyears, fc):
            rows.append({"year": yr, "observed": None, "fitted": None, "forecast": fv})
        out_csv = outdir / f"trend_{vcol}.csv"
        pd.DataFrame(rows).to_csv(out_csv, index=False, encoding="utf-8-sig")
        summary[vcol] = diag
        report += [
            f"## {vcol}",
            "",
            f"- method: **{method}**  |  in-sample RMSE = {diag['rmse']:.4g}  |  MAPE = {diag['mape']:.3f}%"
            f"  |  R2 = {diag['r2']:.4f}",
            f"- horizon: {args.horizon} steps, last observed year {years[-1]:.0f} -> {fyears[-1]:.0f}",
            f"- candidates: " + ", ".join(
                f"{c['method']} (RMSE {c['rmse']:.4g}, R2 {c['r2']:.4f})" for c in diag["candidates"]),
            "",
            "| year | observed | fitted | forecast |",
            "|---|---|---|---|",
        ]
        for r in rows:
            obs = "" if r["observed"] is None else f"{r['observed']:.6g}"
            fitv = "" if r["fitted"] is None else f"{r['fitted']:.6g}"
            fv = "" if r["forecast"] is None else f"{r['forecast']:.6g}"
            report.append(f"| {r['year']:.0f} | {obs} | {fitv} | {fv} |")
        report.append("")
        print(f"OK  {vcol}: method={method} rmse={diag['rmse']:.5g} r2={diag['r2']:.4f} "
              f"-> {out_csv}")

    (outdir / "trend_report.md").write_text("\n".join(report), encoding="utf-8")
    (outdir / "trend_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"    outputs -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
