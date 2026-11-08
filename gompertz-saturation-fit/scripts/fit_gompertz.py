#!/usr/bin/env python3
"""Gompertz saturation model: fit, compare saturation scenarios, extrapolate,
build spatial allocation weights, backtest.

Model (paper notation, Zhang et al. 2020, China Environmental Science):
    y = gamma * exp(alpha * exp(beta * x)),   alpha < 0, beta < 0
Internally the equivalent positive-parameter form is used:
    y = A * exp(-B * exp(-k * x)),   A = gamma,  B = -alpha,  k = -beta

Two numerical safeguards matter:

* x is shifted by x_ref = min(x) before fitting. Without the shift, a
  calendar-year x (values around 2010) makes exp(-k*x) underflow, the curve
  turns numerically flat and the optimiser cannot separate B from k.
* The saturation level A is profiled over a log grid, with B and k obtained
  in closed form at each candidate A (transformed OLS). This avoids the
  local minima that trap a naive multi-start fit when the data sit on the
  accelerating branch of the S-curve.

Only numpy / scipy / pandas / matplotlib are required.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_PLT = True
except Exception:  # pragma: no cover
    HAS_PLT = False


# --------------------------------------------------------------------------- model


def gompertz(x, A, B, k):
    """Gompertz curve. Works with scalars and arrays."""
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", under="ignore"):
        return A * np.exp(-B * np.exp(-k * x))


def logistic(x, A, k, x0):
    """Logistic baseline used only for model comparison."""
    x = np.asarray(x, dtype=float)
    z = -k * (x - x0)
    with np.errstate(over="ignore"):
        return A / (1.0 + np.exp(np.clip(z, -700, 700)))


def predict(fit, x):
    """Evaluate a fitted Gompertz at x, honouring the internal x shift."""
    x = np.asarray(x, dtype=float)
    return gompertz(x - fit.get("x_ref", 0.0), fit["A"], fit["B"], fit["k"])


def jac_natural(x, params, A_fixed=None):
    """Analytic Jacobian of y wrt the free parameters (natural scale)."""
    x = np.asarray(x, dtype=float)
    if A_fixed is None:
        A, B, k = params
    else:
        A, (B, k) = A_fixed, params
    e = np.exp(-k * x)
    y = A * np.exp(-B * e)
    dA = y / A
    dB = -y * e
    dk = y * B * e * x
    if A_fixed is None:
        return np.column_stack([dA, dB, dk])
    return np.column_stack([dB, dk])


# --------------------------------------------------------------------------- estimation


def ols_Bk(x, y, A):
    """Closed-form B, k for a fixed saturation A via transformed OLS.

    ln(ln(A/y)) = ln(B) - k*x  ->  linear in x.
    Returns None when the candidate A is invalid for this data.
    """
    ratio = A / y
    if np.any(ratio <= 1.0):
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.log(np.log(ratio))
    if not np.all(np.isfinite(t)):
        return None
    design = np.vstack([np.ones_like(x), x]).T
    try:
        coef, *_ = np.linalg.lstsq(design, t, rcond=None)
    except np.linalg.LinAlgError:
        return None
    c, d = float(coef[0]), float(coef[1])
    if not np.isfinite(c) or not np.isfinite(d) or d >= 0:
        return None
    B, k = float(np.exp(c)), float(-d)
    if not (np.isfinite(B) and np.isfinite(k)) or k <= 0 or B <= 0:
        return None
    return B, k


def profile_fit(x, y, grid_points=250, mult_hi=1000.0):
    """Profile the saturation A over a log grid; pick the best by log-space SSE."""
    ymax = float(np.max(y))
    A_grid = ymax * np.exp(np.linspace(np.log(1.0005), np.log(mult_hi), grid_points))
    best = None
    for A in A_grid:
        res = ols_Bk(x, y, A)
        if res is None:
            continue
        B, k = res
        yhat = gompertz(x, A, B, k)
        if np.any(yhat <= 0) or not np.all(np.isfinite(yhat)):
            continue
        sse = float(np.sum((np.log(y) - np.log(yhat)) ** 2))
        if best is None or sse < best["sse_log"]:
            best = {"A": float(A), "B": B, "k": k, "sse_log": sse}
    if best is None:
        raise RuntimeError("profile search failed: no valid saturation level found")
    return best


def refine(x, y, p0, A_fixed=None, scale="y"):
    """Bounded least-squares refinement; returns params, covariance, residuals."""
    ymax = float(np.max(y))
    if A_fixed is None:
        bounds = ([ymax * (1.0 + 1e-9), 1e-12, 1e-12], [ymax * 1000.0, 1e8, 1e8])
        start = [p0[0], p0[1], p0[2]]
    else:
        bounds = ([1e-12, 1e-12], [1e8, 1e8])
        start = [p0[1], p0[2]]

    lo = np.array(bounds[0])
    hi = np.array(bounds[1])
    start = [min(max(v, lo[i] * 1.0001), hi[i] * 0.9999) for i, v in enumerate(start)]

    def resid(p):
        yhat = gompertz(x, A_fixed, p[0], p[1]) if A_fixed is not None else gompertz(x, *p)
        yhat = np.clip(yhat, 1e-300, None)
        if scale == "log":
            return np.log(y) - np.log(yhat)
        return y - yhat

    res = least_squares(resid, start, bounds=bounds, method="trf", max_nfev=40000)
    if A_fixed is None:
        A, B, k = (float(v) for v in res.x)
        free = (A, B, k)
    else:
        B, k = (float(v) for v in res.x)
        A = float(A_fixed)
        free = (B, k)

    J = jac_natural(x, free, A_fixed=A_fixed)
    n, p = len(y), J.shape[1]
    dof = max(n - p, 1)
    yhat = gompertz(x, A, B, k)
    r = y - yhat
    s2 = float(np.sum(r ** 2)) / dof
    try:
        cov = np.linalg.pinv(J.T @ J) * s2
    except np.linalg.LinAlgError:
        cov = np.full((p, p), np.nan)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    tcrit = float(stats.t.ppf(0.975, dof))
    return {
        "A": A,
        "B": B,
        "k": k,
        "se": se,
        "cov": cov,
        "tcrit": tcrit,
        "dof": dof,
        "resid": r,
        "yhat": yhat,
        "success": bool(res.success),
    }


def diagnostics(x, y, fit):
    """Fit statistics; x is the (already shifted) design vector."""
    yhat = fit["yhat"]
    r = y - yhat
    ss_res = float(np.sum(r ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    n = len(y)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(np.mean(r ** 2)))
    mape = float(np.mean(np.abs(r / y)) * 100.0)
    p = 2 if fit["gamma_fixed"] else 3
    aic = n * np.log(ss_res / n) + 2 * p if ss_res > 0 else float("inf")
    bic = n * np.log(ss_res / n) + p * np.log(n) if ss_res > 0 else float("inf")
    adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p) if n > p else float("nan")
    A, B, k = fit["A"], fit["B"], fit["k"]
    inflect_x = float(np.log(B / k) / k) if B > 0 and k > 0 else float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        y_inflect = float(A * np.exp(-B / np.e))
    x90 = float(np.log(B / (-np.log(0.9))) / k) if B > 0 and k > 0 else float("nan")
    return {
        "n": n,
        "r2": r2,
        "adj_r2": adj_r2,
        "rmse": rmse,
        "mae": float(np.mean(np.abs(r))),
        "mape": mape,
        "aic": float(aic),
        "bic": float(bic),
        "ss_res": ss_res,
        "inflect_x": inflect_x,
        "y_inflect": y_inflect,
        "x90": x90,
    }


def fit_gompertz(x, y, gamma_fixed=None, scale="y", grid_points=250):
    """Fit one scenario on the shifted design. gamma_fixed=None -> free saturation."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_ref = float(np.min(x))
    xc = x - x_ref

    if gamma_fixed is None:
        prof = profile_fit(xc, y, grid_points=grid_points)
        out = refine(xc, y, [prof["A"], prof["B"], prof["k"]], A_fixed=None, scale=scale)
        out["gamma_fixed"] = False
    else:
        A = float(gamma_fixed)
        seed = ols_Bk(xc, y, A)
        if seed is None:
            # Saturation at or below the observed maximum: the transform is invalid,
            # so start from a generic guess and let least squares decide.
            out = refine(xc, y, [A, 1.0, 0.05], A_fixed=A, scale=scale)
            out["gamma_fixed"] = True
            out["gamma_below_data"] = True
        else:
            out = refine(xc, y, [A, seed[0], seed[1]], A_fixed=A, scale=scale)
            out["gamma_fixed"] = True
            out["gamma_below_data"] = False

    out["x_ref"] = x_ref
    # Original-space alpha (-alpha) is reported when it is numerically representable.
    expo = fit_exponent(out["k"], x_ref)
    out["B_orig"] = out["B"] * expo if np.isfinite(expo) else float("inf")
    out.update(diagnostics(xc, y, out))
    out["inflect_x"] += x_ref  # back to the original x scale
    out["x90"] += x_ref
    return out


def fit_exponent(k, x_ref):
    """exp(k * x_ref), guarded against overflow."""
    z = k * x_ref
    if z > 700:
        return float("inf")
    return float(np.exp(z))


# --------------------------------------------------------------------------- extras


def prediction_band(x_pred, fit, draws=2000, seed=20201105):
    """Percentile band from the parameter covariance (Monte Carlo)."""
    cov = fit["cov"]
    if not np.all(np.isfinite(cov)):
        return np.full(len(x_pred), np.nan), np.full(len(x_pred), np.nan)
    x_ref = fit.get("x_ref", 0.0)
    xs = np.asarray(x_pred, dtype=float) - x_ref
    if fit["gamma_fixed"]:
        mean = np.array([fit["B"], fit["k"]])
        A = fit["A"]
    else:
        mean = np.array([fit["A"], fit["B"], fit["k"]])
        A = None
    rng = np.random.default_rng(seed)
    try:
        samples = rng.multivariate_normal(mean, cov, size=draws, method="cholesky")
    except np.linalg.LinAlgError:
        samples = rng.multivariate_normal(mean, cov, size=draws, method="svd")
    ys = []
    for row in samples:
        if A is None:
            a, b, kk = row
        else:
            a, (b, kk) = A, (row[0], row[1])
        if b <= 0 or kk <= 0 or a <= 0:
            continue
        ys.append(gompertz(xs, a, b, kk))
    if not ys:
        return np.full(len(x_pred), np.nan), np.full(len(x_pred), np.nan)
    arr = np.vstack(ys)
    return np.nanpercentile(arr, 2.5, axis=0), np.nanpercentile(arr, 97.5, axis=0)


def bootstrap_ci(x, y, base_fit, n_boot=200, seed=7):
    """Residual bootstrap; returns percentile CIs for A, B, k."""
    rng = np.random.default_rng(seed)
    x_ref = base_fit.get("x_ref", 0.0)
    xc = np.asarray(x, dtype=float) - x_ref
    yhat = base_fit["yhat"]
    resid = np.asarray(y, dtype=float) - yhat
    A, B, k = base_fit["A"], base_fit["B"], base_fit["k"]
    draws = []
    for _ in range(n_boot):
        yb = np.clip(yhat + rng.choice(resid, size=len(resid), replace=True), 1e-12, None)
        try:
            if base_fit["gamma_fixed"]:
                out = refine(xc, yb, [A, B, k], A_fixed=A)
            else:
                out = refine(xc, yb, [A, B, k], A_fixed=None)
        except Exception:
            continue
        draws.append([out["A"], out["B"], out["k"]])
    if len(draws) < 20:
        return None
    arr = np.vstack(draws)
    return {
        "A": (float(np.percentile(arr[:, 0], 2.5)), float(np.percentile(arr[:, 0], 97.5))),
        "B": (float(np.percentile(arr[:, 1], 2.5)), float(np.percentile(arr[:, 1], 97.5))),
        "k": (float(np.percentile(arr[:, 2], 2.5)), float(np.percentile(arr[:, 2], 97.5))),
        "n_ok": len(draws),
    }


def backtest(x, y, folds=4, gamma_fixed=None, scale="y"):
    """Chronological hold-out: train on the past, predict the next block."""
    order = np.argsort(x)
    x, y = np.asarray(x)[order], np.asarray(y)[order]
    n = len(x)
    if n < 6:
        return []
    blocks = np.array_split(np.arange(n), folds + 1)
    rows = []
    for i in range(1, folds + 1):
        test_idx = blocks[i]
        train_idx = np.concatenate(blocks[:i])
        if len(test_idx) < 1 or len(train_idx) < 4:
            continue
        try:
            fit = fit_gompertz(x[train_idx], y[train_idx], gamma_fixed=gamma_fixed, scale=scale)
        except Exception as exc:
            rows.append({"fold": i, "error": str(exc)})
            continue
        pred = predict(fit, x[test_idx])
        rows.append(
            {
                "fold": i,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "x_test_min": float(x[test_idx].min()),
                "x_test_max": float(x[test_idx].max()),
                "rmse": float(np.sqrt(np.mean((y[test_idx] - pred) ** 2))),
                "mape": float(np.mean(np.abs((y[test_idx] - pred) / y[test_idx])) * 100.0),
            }
        )
    return rows


def fit_logistic(x, y):
    """Logistic baseline for AIC comparison."""
    ymax = float(np.max(y))
    x = np.asarray(x, dtype=float)
    span = max(float(np.ptp(x)), 1e-6)
    A0, k0, x0 = ymax * 1.5, 0.5 / span, float(np.mean(x))
    lo = [ymax * 1.0000001, 1e-9, float(np.min(x)) - 10 * span]
    hi = [ymax * 1e4, 1e4 / span, float(np.max(x)) + 10 * span]
    start = [min(max(v, lo[i] * 1.0001), hi[i] * 0.9999) for i, v in enumerate([A0, k0, x0])]

    def resid(p):
        return y - logistic(x, *p)

    res = least_squares(resid, start, bounds=(lo, hi), method="trf", max_nfev=40000)
    pred = logistic(x, *res.x)
    r = y - pred
    ss_res = float(np.sum(r ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    n = len(y)
    return {
        "A": float(res.x[0]),
        "k": float(res.x[1]),
        "x0": float(res.x[2]),
        "r2": 1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "rmse": float(np.sqrt(np.mean(r ** 2))),
        "mape": float(np.mean(np.abs(r / y)) * 100.0),
        "aic": n * np.log(ss_res / n) + 6 if ss_res > 0 else float("inf"),
        "yhat": pred,
    }


# --------------------------------------------------------------------------- io


def load_table(path):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        try:
            return pd.read_excel(p)
        except Exception as exc:
            raise SystemExit(f"failed to read Excel file: {exc}")
    if p.suffix.lower() == ".csv":
        for enc in ("utf-8-sig", "utf-8", "gb18030", "cp1252", "latin1"):
            try:
                return pd.read_csv(p, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise SystemExit("failed to decode CSV; pass --encoding")
    raise SystemExit(f"unsupported input type: {p.suffix}")


def parse_x_values(args):
    vals = []
    if args.predict_x:
        vals.extend([float(v) for v in str(args.predict_x).replace(";", ",").split(",") if v.strip()])
    if args.predict_range:
        parts = [float(v) for v in str(args.predict_range).split(":")]
        if len(parts) != 3:
            raise SystemExit("--predict-range expects start:stop:step")
        start, stop, step = parts
        if step == 0:
            raise SystemExit("--predict-range step must not be zero")
        count = int(np.floor((stop - start) / step + 1e-9)) + 1
        vals.extend([start + i * step for i in range(count)])
    if args.predict_from is not None:
        start = float(args.predict_from)
        stop = float(args.predict_to) if args.predict_to is not None else start
        step = float(args.predict_step)
        count = int(np.floor((stop - start) / step + 1e-9)) + 1
        vals.extend([start + i * step for i in range(count)])
    return sorted(set(round(v, 10) for v in vals))


# --------------------------------------------------------------------------- reporting


def fmt(v, digits=4):
    if v is None:
        return ""
    if isinstance(v, float) and not np.isfinite(v):
        return "n/a"
    return f"{v:.{digits}f}"


BIG_ALPHA = 1e12  # above this the original-scale alpha is numerically meaningless


def fmt_alpha(fit, digits=4):
    """Report -alpha; fall back to the shifted scale when the original one blows up."""
    b_orig = fit.get("B_orig", fit["B"])
    if not np.isfinite(b_orig) or abs(b_orig) > BIG_ALPHA:
        return f"{-fit['B']:.{digits}f} (shifted)"
    return f"{-b_orig:.{digits}f}"


def build_report(ctx):
    best = ctx["best"]
    lines = []
    a = lines.append
    a("# Gompertz saturation model report")
    a("")
    a(f"- Input: `{ctx['input']}`")
    a(f"- x column: `{ctx['x_col']}`  |  y column: `{ctx['y_col']}`")
    a(f"- Observations used: {best['n']}   x range: [{fmt(ctx['x_min'],3)}, {fmt(ctx['x_max'],3)}]"
      f"   y range: [{fmt(ctx['y_min'],4)}, {fmt(ctx['y_max'],4)}]")
    a("- Model: `y = gamma * exp(alpha * exp(beta * x))`,  alpha < 0, beta < 0")
    a(f"- Selected scenario: **{best['label']}** (by AIC)")
    a(f"- x shift used internally: x_ref = {fmt(best['x_ref'],3)} (parameters below are given for the original x)")
    a("")
    a("## Scenario comparison")
    a("")
    a("| scenario | gamma | alpha | beta | R2 | adj R2 | RMSE | MAPE % | AIC | inflection x | y at inflection |")
    a("|---|---|---|---|---|---|---|---|---|---|---|")
    for row in ctx["rows"]:
        a(
            f"| {row['label']} | {fmt(row['A'],4)} | {fmt_alpha(row,4)} | {fmt(-row['k'],4)} | "
            f"{fmt(row['r2'],4)} | {fmt(row['adj_r2'],4)} | {fmt(row['rmse'],5)} | {fmt(row['mape'],2)} | "
            f"{fmt(row['aic'],2)} | {fmt(row['inflect_x'],3)} | {fmt(row['y_inflect'],4)} |"
        )
    a("")
    a("## Selected model parameters (95% CI)")
    a("")
    a("| parameter | symbol | estimate | se | 95% CI |")
    a("|---|---|---|---|---|")
    if best["gamma_fixed"]:
        a(f"| A | gamma (saturation, fixed) | {fmt(best['A'],4)} | - | fixed input |")
    else:
        lo, hi = best["ci"]["A"]
        a(f"| A | gamma (saturation) | {fmt(best['A'],4)} | {fmt(best['se'][0],4)} | [{fmt(lo,4)}, {fmt(hi,4)}] |")
    off = 0 if best["gamma_fixed"] else 1
    se_idx = {"B": 0, "k": 1}
    for key, label in (("B", "-alpha"), ("k", "-beta")):
        lo, hi = best["ci"][key]
        a(f"| {key} | {label} | {fmt(best[key],4)} | {fmt(best['se'][se_idx[key] + off],4)} | "
          f"[{fmt(lo,4)}, {fmt(hi,4)}] |")
    a("")
    if np.isfinite(best["B_orig"]) and abs(best["B_orig"]) <= BIG_ALPHA:
        a(f"Model in the original x scale: `y = {best['A']:.6g} * exp({-best['B_orig']:.6g} * "
          f"exp({-best['k']:.6g} * x))`")
    else:
        a(f"Model (shifted, x' = x - {best['x_ref']:.6g}): `y = {best['A']:.6g} * exp({-best['B']:.6g} * "
          f"exp({-best['k']:.6g} * x'))`")
        a("")
        a("The original-scale alpha is not representable for this x range (calendar-year style x); "
          "quote the shifted form and state the shift.")
    a("")
    if best.get("bootstrap"):
        a("Residual bootstrap CIs (shifted scale):")
        for key, label in (("A", "gamma"), ("B", "-alpha"), ("k", "-beta")):
            lo, hi = best["bootstrap"][key]
            a(f"- {label}: [{fmt(lo,4)}, {fmt(hi,4)}]  (n_ok={best['bootstrap']['n_ok']})")
        a("")
    a("## Diagnostics")
    a("")
    a(f"- R2 = {fmt(best['r2'],4)}, adj R2 = {fmt(best['adj_r2'],4)}, RMSE = {fmt(best['rmse'],5)},"
      f" MAE = {fmt(best['mae'],5)}, MAPE = {fmt(best['mape'],2)}%")
    a(f"- Inflection point x* = {fmt(best['inflect_x'],3)} (y at inflection = {fmt(best['y_inflect'],4)})")
    a(f"- x where the curve reaches 90% of saturation: {fmt(best['x90'],3)}")
    if ctx.get("logistic"):
        lg = ctx["logistic"]
        a(f"- Logistic baseline: R2 = {fmt(lg['r2'],4)}, RMSE = {fmt(lg['rmse'],5)}, AIC = {fmt(lg['aic'],2)}"
          f" -> {'Gompertz' if best['aic'] < lg['aic'] else 'Logistic'} is preferred by AIC")
    a("")
    if ctx.get("bt"):
        a("## Chronological hold-out validation")
        a("")
        a("| fold | n train | n test | x test range | RMSE | MAPE % |")
        a("|---|---|---|---|---|---|")
        for r in ctx["bt"]:
            if "error" in r:
                a(f"| {r['fold']} | - | - | - | - | failed: {r['error']} |")
            else:
                a(f"| {r['fold']} | {r['n_train']} | {r['n_test']} | "
                  f"[{fmt(r['x_test_min'],3)}, {fmt(r['x_test_max'],3)}] | {fmt(r['rmse'],5)} | {fmt(r['mape'],2)} |")
        a("")
    if ctx.get("warnings"):
        a("## Warnings")
        a("")
        for w in ctx["warnings"]:
            a(f"- {w}")
        a("")
    if ctx.get("weights"):
        a("## Spatial allocation weights")
        a("")
        a(f"Weight = fitted value / sum of fitted values over {ctx['weight_rows']} regions.")
        a("")
        a("| region | x used | fitted y | weight | allocated |")
        a("|---|---|---|---|---|")
        for row in ctx["weights"]:
            alloc = f"{row['allocated']:.4g}" if row.get("allocated") is not None else ""
            a(f"| {row['region']} | {fmt(row['x'],4)} | {fmt(row['predicted'],5)} | {fmt(row['weight'],5)} | {alloc} |")
        a("")
    a("## Files")
    a("")
    for name, desc in ctx["files"]:
        a(f"- `{name}`: {desc}")
    a("")
    return "\n".join(lines)


def make_figure(ctx, path):
    if not HAS_PLT:
        return False
    best = ctx["best"]
    x, y = ctx["x"], ctx["y"]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.scatter(x, y, s=22, color="#2f6fb0", label="observed", zorder=3)
    x_lo = min(float(np.min(x)), float(ctx["x_min"]))
    x_hi = max(float(np.max(x)), float(ctx["x_max"]))
    pad = 0.05 * (x_hi - x_lo) if x_hi > x_lo else 1.0
    xs = np.linspace(x_lo, x_hi + pad, 400)
    ax.plot(xs, predict(best, xs), color="#c0392b", lw=2, label=f"Gompertz fit ({best['label']})", zorder=2)
    lo, hi = prediction_band(xs, best)
    if np.all(np.isfinite(lo)):
        ax.fill_between(xs, lo, hi, color="#c0392b", alpha=0.15, label="95% prediction band")
    ax.axhline(best["A"], color="#7f8c8d", ls="--", lw=1, label=f"saturation gamma={best['A']:.4g}")
    if np.isfinite(best["inflect_x"]):
        ax.axvline(best["inflect_x"], color="#27ae60", ls=":", lw=1.2)
        ax.plot([best["inflect_x"]], [best["y_inflect"]], "o", color="#27ae60", ms=6)
    if ctx.get("pred_x"):
        px = np.array(ctx["pred_x"])
        ax.plot(px, predict(best, px), "^", color="#8e44ad", ms=6, label="forecast points")
    ax.set_xlabel(ctx["x_col"])
    ax.set_ylabel(ctx["y_col"])
    ax.set_title("Gompertz saturation fit")
    ax.grid(alpha=0.25, ls=":")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


# --------------------------------------------------------------------------- main


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fit a Gompertz saturation model and export diagnostics.")
    ap.add_argument("--input", required=True, help="CSV or Excel file with the series")
    ap.add_argument("--x-col", required=True, help="independent variable column (income, year, ...)")
    ap.add_argument("--y-col", required=True, help="dependent variable column (must be > 0)")
    ap.add_argument("--stock-col", help="stock/count column; combined with --pop-col to build y = stock / population")
    ap.add_argument("--pop-col", help="population column, required with --stock-col")
    ap.add_argument("--out", default="gompertz_out", help="output directory")
    ap.add_argument("--gamma", action="append", default=[], type=float,
                    help="fixed saturation scenario; repeatable. Default: none (free saturation only)")
    ap.add_argument("--include-free", action="store_true", default=True,
                    help="always include the free-saturation scenario (default on)")
    ap.add_argument("--no-free", dest="include_free", action="store_false")
    ap.add_argument("--scale", choices=["y", "log"], default="y",
                    help="residual scale used for refinement (default y)")
    ap.add_argument("--grid-points", type=int, default=250, help="profile grid size for the saturation level")
    ap.add_argument("--bootstrap", type=int, default=0, help="bootstrap repetitions for robust CIs")
    ap.add_argument("--backtest", action="store_true", help="run chronological hold-out validation")
    ap.add_argument("--backtest-folds", type=int, default=4)
    ap.add_argument("--compare-logistic", action="store_true", help="also fit a logistic baseline")
    ap.add_argument("--predict-x", help="comma separated x values to predict")
    ap.add_argument("--predict-range", help="start:stop:step of x values to predict")
    ap.add_argument("--predict-from", type=float, help="start of a prediction grid")
    ap.add_argument("--predict-to", type=float, help="end of a prediction grid")
    ap.add_argument("--predict-step", type=float, default=1.0, help="step of the prediction grid")
    ap.add_argument("--weight-col", help="region/category column used to build normalized allocation weights")
    ap.add_argument("--weight-year", help="value of the year column to use for the weights")
    ap.add_argument("--allocate", type=float, help="national total to distribute over the weights")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    df = load_table(args.input)
    for col in (args.x_col, args.y_col):
        if col not in df.columns:
            raise SystemExit(f"column '{col}' not found. Available: {list(df.columns)}")

    work = pd.DataFrame({"x": df[args.x_col], "y": df[args.y_col]})
    if args.stock_col:
        if not args.pop_col:
            raise SystemExit("--stock-col requires --pop-col")
        for col in (args.stock_col, args.pop_col):
            if col not in df.columns:
                raise SystemExit(f"column '{col}' not found")
        work["y"] = pd.to_numeric(df[args.stock_col], errors="coerce") / pd.to_numeric(df[args.pop_col], errors="coerce")
    if args.weight_col:
        if args.weight_col not in df.columns:
            raise SystemExit(f"column '{args.weight_col}' not found")
        work["region"] = df[args.weight_col].astype(str)
    if args.weight_year:
        for cand in ("year", "年份", "Year", "年"):
            if cand in df.columns:
                work["year"] = df[cand]
                break
        else:
            raise SystemExit("--weight-year needs a year-like column (year/年份/Year/年)")

    before = len(work)
    work["x"] = pd.to_numeric(work["x"], errors="coerce")
    work["y"] = pd.to_numeric(work["y"], errors="coerce")
    work = work.dropna(subset=["x", "y"])
    work = work[np.isfinite(work["x"]) & np.isfinite(work["y"])]
    work = work[work["y"] > 0]
    dropped = before - len(work)
    if len(work) < 5:
        raise SystemExit(f"only {len(work)} usable rows; need at least 5 positive observations")

    x = work["x"].astype(float).to_numpy()
    y = work["y"].astype(float).to_numpy()
    if np.ptp(x) == 0:
        raise SystemExit("x has no variation")

    scenarios = []
    if args.include_free:
        scenarios.append((None, "free gamma"))
    for g in args.gamma:
        scenarios.append((float(g), f"gamma={g:g}"))
    if not scenarios:
        raise SystemExit("nothing to fit")

    rows, fits = [], {}
    warnings = []
    for g, label in scenarios:
        try:
            fit = fit_gompertz(x, y, gamma_fixed=g, scale=args.scale, grid_points=args.grid_points)
        except Exception as exc:
            warnings.append(f"scenario {label} failed: {exc}")
            continue
        fit["label"] = label
        fits[label] = fit
        row = {k: fit[k] for k in ("A", "B", "B_orig", "k", "r2", "adj_r2", "rmse", "mae", "mape", "aic", "bic",
                                   "inflect_x", "y_inflect", "x90", "n", "gamma_fixed", "x_ref")}
        row["label"] = label
        rows.append(row)

    if not rows:
        raise SystemExit("all scenarios failed:\n" + "\n".join(warnings))

    best_row = min(rows, key=lambda r: r["aic"] if np.isfinite(r["aic"]) else np.inf)
    best = fits[best_row["label"]]

    tcrit = best["tcrit"]
    se = best["se"]
    if best["gamma_fixed"]:
        ci = {"A": (best["A"], best["A"]),
              "B": (best["B"] - tcrit * se[0], best["B"] + tcrit * se[0]),
              "k": (best["k"] - tcrit * se[1], best["k"] + tcrit * se[1])}
    else:
        ci = {
            "A": (best["A"] - tcrit * se[0], best["A"] + tcrit * se[0]),
            "B": (best["B"] - tcrit * se[1], best["B"] + tcrit * se[1]),
            "k": (best["k"] - tcrit * se[2], best["k"] + tcrit * se[2]),
        }
    best["ci"] = ci
    if args.bootstrap > 0:
        best["bootstrap"] = bootstrap_ci(x, y, best, n_boot=args.bootstrap)

    # Identifiability checks
    if np.isfinite(best["inflect_x"]) and best["inflect_x"] > x.max():
        warnings.append(
            "The inflection point lies beyond the observed x range: the saturation level is "
            "extrapolated, not identified by the data. Treat gamma as a scenario assumption.")
    if not best["gamma_fixed"] and np.isfinite(ci["A"][1]) and ci["A"][1] / max(ci["A"][0], 1e-12) > 10:
        warnings.append("The 95% CI of gamma spans an order of magnitude: the upper branch of the "
                        "S-curve is unconstrained by the data.")
    if best["k"] < 1e-6:
        warnings.append("beta (k) is near zero: the curve behaves like a straight line in log space.")
    if best.get("gamma_below_data"):
        warnings.append("The fixed gamma is at or below the observed maximum of y: the fit is not "
                        "meaningful for this data.")
    if best["r2"] < 0.9:
        warnings.append(f"Poor fit (R2 = {best['r2']:.3f}). Check the x variable, units and the "
                        "assumption of a single saturation level.")

    logistic_fit = fit_logistic(x, y) if args.compare_logistic else None

    bt = backtest(x, y, folds=args.backtest_folds,
                  gamma_fixed=best["A"] if best["gamma_fixed"] else None,
                  scale=args.scale) if args.backtest else []

    pred_x = parse_x_values(args)
    pred_rows = []
    if pred_x:
        px = np.array(pred_x, dtype=float)
        py = predict(best, px)
        lo, hi = prediction_band(px, best)
        for i, xv in enumerate(px):
            pred_rows.append({"x": xv, "y_pred": float(py[i]),
                              "y_lower95": float(lo[i]) if np.isfinite(lo[i]) else None,
                              "y_upper95": float(hi[i]) if np.isfinite(hi[i]) else None})

    weight_rows = []
    if args.weight_col:
        wdf = work.copy()
        if args.weight_year and "year" in wdf.columns:
            wdf = wdf[wdf["year"].astype(str) == str(args.weight_year)]
            if wdf.empty:
                raise SystemExit(f"no rows match --weight-year {args.weight_year}")
        idx = wdf.groupby("region")["x"].idxmax()
        sel = wdf.loc[idx]
        vals = predict(best, sel["x"].to_numpy())
        total = float(np.sum(vals))
        for reg, xv, yv in zip(sel["region"].tolist(), sel["x"].tolist(), vals.tolist()):
            w = float(yv / total) if total > 0 else 0.0
            weight_rows.append({"region": reg, "x": float(xv), "predicted": float(yv), "weight": w,
                                "allocated": (args.allocate * w) if args.allocate else None})

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    files = []
    pd.DataFrame(rows).to_csv(outdir / "parameters.csv", index=False, encoding="utf-8-sig")
    files.append(("parameters.csv", "per-scenario parameters and fit statistics"))

    if pred_rows:
        pd.DataFrame(pred_rows).to_csv(outdir / "predictions.csv", index=False, encoding="utf-8-sig")
        files.append(("predictions.csv", "forecast values with 95% prediction band"))

    if weight_rows:
        pd.DataFrame(weight_rows).to_csv(outdir / "spatial_weights.csv", index=False, encoding="utf-8-sig")
        files.append(("spatial_weights.csv", "normalized per-region allocation weights"))

    if bt:
        pd.DataFrame(bt).to_csv(outdir / "backtest.csv", index=False, encoding="utf-8-sig")
        files.append(("backtest.csv", "chronological hold-out diagnostics"))

    if not args.no_figures and HAS_PLT:
        if make_figure({"best": best, "x": x, "y": y, "x_min": float(x.min()), "x_max": float(x.max()),
                        "x_col": args.x_col, "y_col": args.y_col, "pred_x": pred_x},
                       outdir / "gompertz_fit.png"):
            files.append(("gompertz_fit.png", "observed data, fitted curve and prediction band"))

    ctx = {
        "input": args.input,
        "x_col": args.x_col,
        "y_col": args.y_col,
        "x": x, "y": y,
        "x_min": float(x.min()), "x_max": float(x.max()),
        "y_min": float(y.min()), "y_max": float(y.max()),
        "rows": rows,
        "best": best,
        "logistic": logistic_fit,
        "bt": bt,
        "warnings": warnings,
        "weights": weight_rows,
        "weight_rows": len(weight_rows),
        "files": files,
        "pred_x": pred_x,
    }
    files.append(("fit_report.md", "methods, diagnostics, warnings"))
    files.append(("fit_summary.json", "machine-readable run summary"))
    (outdir / "fit_report.md").write_text(build_report(ctx), encoding="utf-8")

    def num(v):
        return float(v) if isinstance(v, (int, float, np.floating)) and np.isfinite(v) else None

    summary = {
        "input": args.input,
        "x_col": args.x_col,
        "y_col": args.y_col,
        "n": int(len(y)),
        "dropped_rows": int(dropped),
        "selected_scenario": best["label"],
        "best": {k: num(best[k]) for k in ("A", "B", "k", "r2", "adj_r2", "rmse", "mape", "aic",
                                           "inflect_x", "y_inflect", "x_ref")},
        "alpha_original_scale": num(-best["B_orig"]) if (np.isfinite(best["B_orig"]) and
                                                        abs(best["B_orig"]) <= BIG_ALPHA) else None,
        "ci95": {k: [num(v[0]), num(v[1])] for k, v in ci.items()},
        "scenarios": [{k: num(r[k]) for k in ("label", "A", "k", "r2", "rmse", "mape", "aic")}
                      for r in rows],
        "warnings": warnings,
    }
    (outdir / "fit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"OK  n={len(y)}  dropped={dropped}  scenario='{best['label']}'")
    alpha_txt = (f"{-best['B_orig']:.5f}" if np.isfinite(best["B_orig"]) and abs(best["B_orig"]) <= BIG_ALPHA
                 else f"{-best['B']:.5f} (shifted, x' = x - {best['x_ref']:.4g})")
    print(f"    gamma={best['A']:.5f}  alpha={alpha_txt}  beta={-best['k']:.5f}")
    print(f"    R2={best['r2']:.4f}  RMSE={best['rmse']:.5f}  MAPE={best['mape']:.2f}%  AIC={best['aic']:.2f}")
    print(f"    inflection x*={best['inflect_x']:.3f}  y*={best['y_inflect']:.4f}")
    for w in warnings:
        print(f"    WARN {w}")
    print(f"    outputs -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
