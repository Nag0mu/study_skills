#!/usr/bin/env python3
"""Complete an unbalanced city-year vehicle ownership panel.

The script preserves observed values, fits logistic or Gompertz saturation
curves, falls back to pooled curves for sparse cities, optionally reconciles
the completed panel to known yearly totals, and writes auditable outputs.
"""

from __future__ import annotations

import argparse
import codecs
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ENCODING_CANDIDATES = (
    "utf-8",
    "cp1252",
    "iso-8859-1",
    "gb18030",
    "big5",
    "shift_jis",
)
MISSING_TOKENS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    ".",
    "..",
    "...",
    "-",
    "--",
    "—",
    "–",
    "n.a.",
}
YEAR_PATTERN = re.compile(r"(?<!\d)((?:19|20)\d{2})(?!\d)")


@dataclass
class CurveFit:
    model: str
    params: np.ndarray
    rmse: float
    aicc: float
    nobs: int
    residual_sigma: float
    covariance: np.ndarray | None

    def predict(self, years: np.ndarray) -> np.ndarray:
        return predict_curve(self.model, self.params, years)

    def interval(self, years: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return prediction_interval(self, years)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Complete an unbalanced city-year vehicle ownership panel."
    )
    parser.add_argument("--input", required=True, help="Input CSV/XLSX file.")
    parser.add_argument("--out", required=True, help="Output directory.")
    parser.add_argument(
        "--format",
        choices=["auto", "long", "wide"],
        default="auto",
        help="Input layout. Wide means city rows and year columns.",
    )
    parser.add_argument("--sheet", default=None, help="Excel sheet name.")
    parser.add_argument(
        "--encoding",
        default="auto",
        help="Input text encoding, or 'auto'. Supports UTF-8 BOM, UTF-16, UTF-32, CP1252, Latin-1, and other Python codecs.",
    )
    parser.add_argument(
        "--output-encoding",
        default=None,
        help="CSV output encoding. Defaults to a Unicode-safe derivative of the detected input encoding.",
    )
    parser.add_argument("--city-col", default="city", help="City identifier column.")
    parser.add_argument("--year-col", default="year", help="Year column.")
    parser.add_argument(
        "--value-col",
        default="vehicles",
        help="Vehicle ownership/stock column in long format.",
    )
    parser.add_argument(
        "--population-col",
        default=None,
        help="Population column. Used for rates and count reconciliation.",
    )
    parser.add_argument(
        "--gdp-col",
        default=None,
        help="Optional GDP-per-capita column for sparse-city neighbour matching.",
    )
    parser.add_argument(
        "--total-col",
        default=None,
        help="Yearly regional/national total repeated in the input file.",
    )
    parser.add_argument(
        "--total-file",
        default=None,
        help="Optional CSV/XLSX file with yearly totals.",
    )
    parser.add_argument(
        "--total-sheet", default=None, help="Excel sheet name for --total-file."
    )
    parser.add_argument(
        "--total-year-col", default="year", help="Year column in --total-file."
    )
    parser.add_argument(
        "--total-value-col",
        default="total",
        help="Total-value column in --total-file.",
    )
    parser.add_argument(
        "--value-unit",
        choices=["counts", "per_1000"],
        default="counts",
        help="Unit of --value-col. Counts are total vehicles; per_1000 is already normalized.",
    )
    parser.add_argument(
        "--basis",
        choices=["auto", "raw", "rate"],
        default="auto",
        help="Fit raw levels or vehicles per 1000 population.",
    )
    parser.add_argument(
        "--model",
        choices=["auto", "logistic", "gompertz"],
        default="auto",
        help="Curve family. Auto compares AICc when both converge.",
    )
    parser.add_argument(
        "--min-obs",
        type=int,
        default=5,
        help="Minimum observed years for a city-specific curve.",
    )
    parser.add_argument(
        "--min-span",
        type=int,
        default=3,
        help="Minimum observed year span for a city-specific curve.",
    )
    parser.add_argument(
        "--aggregate-duplicates",
        choices=["sum", "mean"],
        default="sum",
        help="How to aggregate duplicate city-year observations.",
    )
    parser.add_argument("--start-year", type=int, default=None)
    parser.add_argument("--end-year", type=int, default=None)
    parser.add_argument(
        "--extrapolate",
        action="store_true",
        help="Fill years outside the observed range when --start-year/--end-year extend it.",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip matplotlib figures.",
    )
    parser.add_argument(
        "--no-excel",
        action="store_true",
        help="Skip the Unicode-safe XLSX workbook output.",
    )
    return parser.parse_args()


def clean_string(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def to_numeric_series(series: pd.Series) -> pd.Series:
    text = clean_string(series).str.lower()
    text = text.mask(text.isin(MISSING_TOKENS))
    return pd.to_numeric(text, errors="coerce")


def parse_year(value: Any) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    match = YEAR_PATTERN.search(text)
    if match:
        return float(match.group(1))
    try:
        number = float(text)
    except ValueError:
        return np.nan
    if 1900 <= number <= 2200:
        return number
    return np.nan


def normalize_encoding_name(encoding: str) -> str:
    return codecs.lookup(encoding).name


def western_latin_score(text: str) -> int:
    western = set(
        "ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÑÒÓÔÕÖØÙÚÛÜÝŸ"
        "àáâãäåæçèéêëìíîïñòóôõöøùúûüýÿ"
    )
    return sum(character in western for character in text)


def detect_text_encoding(
    path: Path, requested: str
) -> tuple[str, dict[str, Any]]:
    data = path.read_bytes()
    if not data:
        raise ValueError(f"Input file is empty: {path}")

    if requested.lower() != "auto":
        encoding = normalize_encoding_name(requested)
        try:
            data.decode(encoding, errors="strict")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"File cannot be decoded as {encoding}: {exc}"
            ) from exc
        return encoding, {
            "encoding": encoding,
            "method": "explicit",
            "confidence": 1.0,
        }

    if data.startswith(codecs.BOM_UTF32_LE) or data.startswith(codecs.BOM_UTF32_BE):
        encoding = normalize_encoding_name("utf-32")
    elif data.startswith(codecs.BOM_UTF16_LE) or data.startswith(codecs.BOM_UTF16_BE):
        encoding = normalize_encoding_name("utf-16")
    elif data.startswith(codecs.BOM_UTF8):
        encoding = normalize_encoding_name("utf-8-sig")
    else:
        encoding = ""

    if encoding:
        try:
            data.decode(encoding, errors="strict")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"Encoding BOM indicated {encoding}, but strict decoding failed: {exc}"
            ) from exc
        return encoding, {
            "encoding": encoding,
            "method": "bom",
            "confidence": 1.0,
        }

    try:
        data.decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        pass
    else:
        return "utf-8", {
            "encoding": "utf-8",
            "method": "strict-utf8",
            "confidence": 1.0,
        }

    try:
        from charset_normalizer import from_bytes

        best = from_bytes(data).best()
        if best is not None and best.encoding:
            candidate = normalize_encoding_name(best.encoding)
            decoded = str(best)
            if "\ufffd" not in decoded:
                try:
                    data.decode(candidate, errors="strict")
                except UnicodeDecodeError:
                    pass
                else:
                    method = "charset-normalizer"
                    if candidate in {"cp1250", "iso8859-2"}:
                        try:
                            western_text = data.decode("cp1252", errors="strict")
                        except UnicodeDecodeError:
                            western_text = ""
                        if western_latin_score(western_text) > western_latin_score(
                            decoded
                        ):
                            candidate = "cp1252"
                            decoded = western_text
                            method = "charset-normalizer+western-latin"
                    chaos = getattr(best, "percent_chaos", None)
                    confidence = (
                        max(0.0, 1.0 - float(chaos) / 100.0)
                        if chaos is not None
                        else 0.75
                    )
                    return candidate, {
                        "encoding": candidate,
                        "method": method,
                        "confidence": confidence,
                    }
    except ImportError:
        pass

    for candidate in ENCODING_CANDIDATES:
        try:
            decoded = data.decode(candidate, errors="strict")
        except UnicodeDecodeError:
            continue
        if "\x00" in decoded:
            continue
        return normalize_encoding_name(candidate), {
            "encoding": normalize_encoding_name(candidate),
            "method": "strict-fallback",
            "confidence": 0.6,
        }

    raise ValueError(
        "Could not detect the text encoding safely. Re-run with an explicit codec, "
        "for example --encoding cp1252, --encoding iso-8859-1, or --encoding gb18030."
    )


def read_table(
    path: Path, sheet: str | None, encoding: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".txt"}:
        detected, info = detect_text_encoding(path, encoding)
        frame = pd.read_csv(path, encoding=detected)
    elif suffix == ".tsv":
        detected, info = detect_text_encoding(path, encoding)
        frame = pd.read_csv(path, sep="\t", encoding=detected)
    elif suffix in {".xlsx", ".xls", ".xlsm"}:
        frame = pd.read_excel(path, sheet_name=sheet or 0)
        info = {
            "encoding": "unicode-xlsx",
            "method": "excel-unicode",
            "confidence": 1.0,
        }
    else:
        raise ValueError(f"Unsupported input type: {path.suffix}")
    frame.attrs["encoding_info"] = info
    return frame, info


def infer_format(df: pd.DataFrame, year_col: str, requested: str) -> str:
    if requested != "auto":
        return requested
    if year_col in df.columns:
        return "long"
    year_columns = [col for col in df.columns if YEAR_PATTERN.search(str(col))]
    if year_columns:
        return "wide"
    raise ValueError(
        f"Cannot infer panel format: '{year_col}' is absent and no year columns were found."
    )


def wide_to_long(
    df: pd.DataFrame, city_col: str, value_unit: str
) -> pd.DataFrame:
    year_columns = [
        col for col in df.columns if YEAR_PATTERN.search(str(col))
    ]
    if city_col not in df.columns:
        raise ValueError(f"City column '{city_col}' was not found.")
    if not year_columns:
        raise ValueError("No four-digit year columns were found for wide input.")
    work = df[[city_col] + year_columns].copy()
    work = work.rename(columns={city_col: "city"})
    long = work.melt(
        id_vars=["city"],
        value_vars=year_columns,
        var_name="year",
        value_name="vehicles",
    )
    long["year"] = long["year"].map(parse_year)
    long["vehicles"] = to_numeric_series(long["vehicles"])
    long.attrs["value_unit"] = value_unit
    return long


def load_panel(
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, str, dict[str, Any]]:
    input_path = Path(args.input).expanduser().resolve()
    raw, encoding_info = read_table(input_path, args.sheet, args.encoding)
    panel_format = infer_format(raw, args.year_col, args.format)

    if panel_format == "wide":
        panel = wide_to_long(raw, args.city_col, args.value_unit)
        if args.population_col and args.population_col in raw.columns:
            population = raw[[args.city_col, args.population_col]].copy()
            population = population.rename(
                columns={args.city_col: "city", args.population_col: "population"}
            )
            population["population"] = to_numeric_series(population["population"])
            population = population.dropna(subset=["population"]).drop_duplicates("city")
            panel = panel.merge(population, on="city", how="left")
        if args.gdp_col and args.gdp_col in raw.columns:
            gdp = raw[[args.city_col, args.gdp_col]].copy()
            gdp = gdp.rename(columns={args.city_col: "city", args.gdp_col: "gdp"})
            gdp["gdp"] = to_numeric_series(gdp["gdp"])
            gdp = gdp.dropna(subset=["gdp"]).drop_duplicates("city")
            panel = panel.merge(gdp, on="city", how="left")
        return panel, panel_format, encoding_info

    required = [args.city_col, args.year_col, args.value_col]
    missing_cols = [col for col in required if col not in raw.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    selected = {
        args.city_col: "city",
        args.year_col: "year",
        args.value_col: "vehicles",
    }
    optional = {
        args.population_col: "population",
        args.gdp_col: "gdp",
        args.total_col: "total",
    }
    for source, target in optional.items():
        if source and source in raw.columns:
            selected[source] = target
    panel = raw[list(selected)].rename(columns=selected).copy()
    panel["city"] = clean_string(panel["city"])
    panel["year"] = panel["year"].map(parse_year)
    panel["vehicles"] = to_numeric_series(panel["vehicles"])
    for column in ["population", "gdp", "total"]:
        if column in panel.columns:
            panel[column] = to_numeric_series(panel[column])
    return panel, panel_format, encoding_info


def load_external_totals(
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = Path(args.total_file).expanduser().resolve()
    frame, encoding_info = read_table(path, args.total_sheet, args.encoding)
    required = [args.total_year_col, args.total_value_col]
    missing_cols = [col for col in required if col not in frame.columns]
    if missing_cols:
        raise ValueError(f"Missing total-file columns: {missing_cols}")
    totals = frame[required].rename(
        columns={args.total_year_col: "year", args.total_value_col: "total"}
    )
    totals["year"] = totals["year"].map(parse_year)
    totals["total"] = to_numeric_series(totals["total"])
    result = totals.dropna(subset=["year", "total"]).groupby(
        "year", as_index=False
    )["total"].median()
    return result, encoding_info


def extract_totals(
    panel: pd.DataFrame, args: argparse.Namespace
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    total_frames: list[pd.DataFrame] = []
    encoding_infos: list[dict[str, Any]] = []
    if "total" in panel.columns:
        same_file = panel[["year", "total"]].dropna()
        if not same_file.empty:
            total_frames.append(
                same_file.groupby("year", as_index=False)["total"].median()
            )
    if args.total_file:
        external, encoding_info = load_external_totals(args)
        total_frames.append(external)
        encoding_infos.append(encoding_info)
    if not total_frames:
        return pd.DataFrame(columns=["year", "total"]), encoding_infos
    combined = pd.concat(total_frames, ignore_index=True)
    return (
        combined.groupby("year", as_index=False)["total"].median(),
        encoding_infos,
    )


def aggregate_duplicates(panel: pd.DataFrame, method: str) -> pd.DataFrame:
    subset = panel.dropna(subset=["year"])
    duplicated = subset.duplicated(["city", "year"], keep=False)
    if not duplicated.any():
        return subset
    numeric_cols = [
        col
        for col in ["vehicles", "population", "gdp", "total"]
        if col in subset.columns
    ]
    group_cols = ["city", "year"]
    if method == "sum":
        aggregator = {
            "vehicles": "sum",
            "population": "mean",
            "gdp": "mean",
            "total": "median",
        }
    else:
        aggregator = {col: "mean" for col in numeric_cols}
    aggregator = {
        col: operation
        for col, operation in aggregator.items()
        if col in subset.columns
    }
    return subset.groupby(group_cols, as_index=False).agg(aggregator)


def determine_year_range(
    panel: pd.DataFrame, args: argparse.Namespace
) -> tuple[int, int]:
    observed_years = panel["year"].dropna()
    if observed_years.empty:
        raise ValueError("No valid years were found.")
    observed_min = int(observed_years.min())
    observed_max = int(observed_years.max())
    start = args.start_year or observed_min
    end = args.end_year or observed_max
    if end < start:
        raise ValueError("--end-year must be greater than or equal to --start-year.")
    if not args.extrapolate:
        start = max(start, observed_min)
        end = min(end, observed_max)
    return start, end


def complete_grid(panel: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    cities = sorted(panel["city"].dropna().unique())
    years = list(range(start, end + 1))
    grid = pd.MultiIndex.from_product([cities, years], names=["city", "year"]).to_frame(
        index=False
    )
    merged = grid.merge(panel, on=["city", "year"], how="left")
    return merged.sort_values(["city", "year"]).reset_index(drop=True)


def interpolate_population(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    if "population" not in df.columns:
        return df, 0
    result = df.copy()
    imputed = 0
    for city, indices in result.groupby("city").groups.items():
        city_rows = result.loc[indices].sort_values("year")
        series = city_rows.set_index("year")["population"]
        observed = series.dropna()
        if observed.empty:
            continue
        filled = series.interpolate(method="index", limit_direction="both")
        newly = series.isna() & filled.notna()
        imputed += int(newly.sum())
        result.loc[city_rows.index, "population"] = filled.to_numpy()
    return result, imputed


def determine_basis(df: pd.DataFrame, args: argparse.Namespace) -> str:
    if args.basis != "auto":
        return args.basis
    if (
        "population" in df.columns
        and df["population"].notna().any()
        and args.value_unit == "counts"
    ):
        return "rate"
    return "raw"


def add_model_basis(
    df: pd.DataFrame, basis: str, value_unit: str
) -> tuple[pd.DataFrame, list[str]]:
    result = df.copy()
    warnings: list[str] = []
    if value_unit == "per_1000":
        result["model_value"] = result["vehicles"]
        return result, warnings
    if basis == "raw":
        result["model_value"] = result["vehicles"]
        return result, warnings
    valid_population = result["population"].notna() & (result["population"] > 0)
    result["model_value"] = np.where(
        valid_population,
        1000.0 * result["vehicles"] / result["population"],
        np.nan,
    )
    lost = int(
        result["vehicles"].notna().sum()
        - result["vehicles"].notna().mul(valid_population).sum()
    )
    if lost:
        warnings.append(
            f"{lost} observed vehicle values lacked a positive population and were excluded from rate fitting."
        )
    return result, warnings


def logistic_curve(years: np.ndarray, params: np.ndarray) -> np.ndarray:
    saturation, growth, inflection = params
    exponent = np.clip(-growth * (years - inflection), -700, 700)
    return saturation / (1.0 + np.exp(exponent))


def gompertz_curve(years: np.ndarray, params: np.ndarray) -> np.ndarray:
    saturation, growth, inflection = params
    exponent = np.clip(-growth * (years - inflection), -700, 700)
    return saturation * np.exp(-np.exp(exponent))


def predict_curve(model: str, params: np.ndarray, years: np.ndarray) -> np.ndarray:
    if model == "logistic":
        return logistic_curve(years, params)
    if model == "gompertz":
        return gompertz_curve(years, params)
    raise ValueError(f"Unknown curve model: {model}")


def curve_residuals(
    params: np.ndarray, years: np.ndarray, values: np.ndarray, model: str
) -> np.ndarray:
    return predict_curve(model, params, years) - values


def parameter_bounds(
    years: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    observed_max = float(np.nanmax(values))
    observed_min = float(np.nanmin(values))
    span = float(np.nanmax(years) - np.nanmin(years))
    scale = max(observed_max, observed_min, 1e-6)
    lower_k = max(scale * 1.001, 1e-9)
    upper_k = max(scale * 10.0 + 1.0, lower_k * 2.0)
    pad = max(20.0, span * 2.0)
    lower = np.array([lower_k, 1e-5, float(np.min(years)) - pad])
    upper = np.array([upper_k, 5.0, float(np.max(years)) + pad])
    return lower, upper


def initial_starts(
    years: np.ndarray, values: np.ndarray
) -> list[np.ndarray]:
    max_value = float(np.nanmax(values))
    median_year = float(np.nanmedian(years))
    starts = []
    for factor in (1.05, 1.5, 3.0):
        for growth in (0.05, 0.2, 0.5):
            starts.append(
                np.array([max(max_value * factor, 1e-8), growth, median_year])
            )
    return starts


def fit_curve(
    years: np.ndarray,
    values: np.ndarray,
    model: str,
) -> CurveFit | None:
    valid = np.isfinite(years) & np.isfinite(values) & (values >= 0)
    years = np.asarray(years[valid], dtype=float)
    values = np.asarray(values[valid], dtype=float)
    if len(values) < 4 or np.nanmax(values) <= 0:
        return None
    lower, upper = parameter_bounds(years, values)
    best = None
    for start in initial_starts(years, values):
        start = np.clip(start, lower + 1e-9, upper - 1e-9)
        try:
            result = least_squares(
                curve_residuals,
                start,
                bounds=(lower, upper),
                args=(years, values, model),
                max_nfev=5000,
                ftol=1e-10,
                xtol=1e-10,
                gtol=1e-10,
            )
        except (ValueError, FloatingPointError):
            continue
        if not result.success:
            continue
        residuals = curve_residuals(result.x, years, values, model)
        rss = float(np.sum(residuals**2))
        if best is None or rss < best[0]:
            best = (rss, result)
    if best is None:
        return None
    rss, result = best
    nobs = len(values)
    sigma = math.sqrt(rss / max(nobs - 3, 1))
    aicc = nobs * math.log(max(rss / nobs, 1e-20)) + 2 * 3
    if nobs - 3 - 1 > 0:
        aicc += (2 * 3 * (3 + 1)) / (nobs - 3 - 1)
    covariance = None
    if result.jac is not None and nobs > 3:
        try:
            covariance = np.linalg.pinv(result.jac.T @ result.jac) * sigma**2
        except np.linalg.LinAlgError:
            covariance = None
    return CurveFit(
        model=model,
        params=result.x,
        rmse=math.sqrt(rss / nobs),
        aicc=float(aicc),
        nobs=nobs,
        residual_sigma=float(sigma),
        covariance=covariance,
    )


def prediction_gradient(
    model: str, params: np.ndarray, years: np.ndarray
) -> np.ndarray:
    saturation, growth, inflection = params
    epsilon = 1e-6
    base = predict_curve(model, params, years)
    gradient = np.empty((len(years), 3), dtype=float)
    for index in range(3):
        shifted = params.copy()
        shifted[index] += epsilon * max(abs(params[index]), 1.0)
        gradient[:, index] = (
            predict_curve(model, shifted, years) - base
        ) / (epsilon * max(abs(params[index]), 1.0))
    return gradient


def prediction_interval(
    fit: CurveFit, years: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    years = np.asarray(years, dtype=float)
    prediction = fit.predict(years)
    variance = np.full(len(years), fit.residual_sigma**2)
    if fit.covariance is not None:
        gradient = prediction_gradient(fit.model, fit.params, years)
        variance += np.einsum("ij,jk,ik->i", gradient, fit.covariance, gradient)
    half_width = 1.96 * np.sqrt(np.clip(variance, 0, None))
    return np.clip(prediction - half_width, 0, None), prediction + half_width


def choose_fit(
    years: np.ndarray,
    values: np.ndarray,
    requested_model: str,
) -> CurveFit | None:
    candidates: list[CurveFit] = []
    models = (
        ["logistic", "gompertz"]
        if requested_model == "auto"
        else [requested_model]
    )
    for model in models:
        fit = fit_curve(years, values, model)
        if fit is not None:
            candidates.append(fit)
    if not candidates:
        return None
    return min(candidates, key=lambda item: item.aicc)


def city_scales(
    city: str,
    years: np.ndarray,
    values: np.ndarray,
    pooled_fit: CurveFit,
    city_features: dict[str, np.ndarray],
    scale_lookup: dict[str, float],
) -> tuple[float, float]:
    base = pooled_fit.predict(years)
    valid = np.isfinite(values) & np.isfinite(base) & (base > 1e-12)
    if valid.any():
        ratios = values[valid] / base[valid]
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        if len(ratios):
            return float(np.median(ratios)), 0.75
    observed_scales = [
        (name, scale)
        for name, scale in scale_lookup.items()
        if name != city and np.isfinite(scale) and scale > 0
    ]
    if not observed_scales:
        return 1.0, 0.1
    target_features = city_features.get(city)
    if target_features is None or not np.isfinite(target_features).any():
        return float(np.median([item[1] for item in observed_scales])), 0.35
    names = [item[0] for item in observed_scales]
    matrix = np.vstack([city_features[name] for name in names])
    target = target_features
    valid_columns = np.isfinite(matrix).all(axis=0) & np.isfinite(target)
    if not valid_columns.any():
        return float(np.median([item[1] for item in observed_scales])), 0.35
    matrix = matrix[:, valid_columns]
    target = target[valid_columns]
    std = matrix.std(axis=0)
    std[std == 0] = 1.0
    distances = np.sqrt(np.sum(((matrix - target) / std) ** 2, axis=1))
    nearest = np.argsort(distances)[: min(5, len(distances))]
    weights = 1.0 / (distances[nearest] + 0.25)
    values = np.array([observed_scales[index][1] for index in nearest])
    return float(np.average(values, weights=weights)), 0.45


def fit_all_cities(
    df: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    city_rows: list[dict[str, Any]] = []
    individual_fits: dict[str, CurveFit] = {}
    normalized: list[dict[str, Any]] = []
    city_feature_rows: list[dict[str, Any]] = []

    for city, group in df.groupby("city", sort=True):
        observed = group.dropna(subset=["model_value"])
        years = observed["year"].to_numpy(dtype=float)
        values = observed["model_value"].to_numpy(dtype=float)
        feature_values = []
        if "population" in observed.columns and observed["population"].notna().any():
            feature_values.append(
                math.log1p(max(float(observed["population"].median()), 0))
            )
        if "gdp" in observed.columns and observed["gdp"].notna().any():
            feature_values.append(math.log1p(max(float(observed["gdp"].median()), 0)))
        city_feature_rows.append(
            {
                "city": city,
                "features": np.asarray(feature_values, dtype=float)
                if feature_values
                else None,
            }
        )
        span = float(np.nanmax(years) - np.nanmin(years)) if len(years) else 0
        fit = None
        if len(observed) >= args.min_obs and span >= args.min_span:
            fit = choose_fit(years, values, args.model)
        if fit is not None:
            individual_fits[city] = fit
            city_rows.append(
                {
                    "city": city,
                    "n_obs": len(observed),
                    "first_year": int(years.min()),
                    "last_year": int(years.max()),
                    "fit_type": "city_specific",
                    "selected_model": fit.model,
                    "K": fit.params[0],
                    "growth": fit.params[1],
                    "inflection": fit.params[2],
                    "rmse": fit.rmse,
                    "aicc": fit.aicc,
                    "scale": 1.0,
                    "confidence": 0.9,
                }
            )
        else:
            city_rows.append(
                {
                    "city": city,
                    "n_obs": len(observed),
                    "first_year": int(years.min()) if len(years) else None,
                    "last_year": int(years.max()) if len(years) else None,
                    "fit_type": "pooled",
                    "selected_model": None,
                    "K": np.nan,
                    "growth": np.nan,
                    "inflection": np.nan,
                    "rmse": np.nan,
                    "aicc": np.nan,
                    "scale": np.nan,
                    "confidence": 0.45 if len(observed) else 0.2,
                }
            )
        if len(values) >= 2 and np.nanmax(values) > 0:
            scale = float(np.median(values[values > 0]))
            if np.isfinite(scale) and scale > 0:
                normalized.append(
                    {
                        "city": city,
                        "year": years,
                        "value": values / scale,
                    }
                )

    pooled_fit = None
    if normalized:
        pooled_years = np.concatenate([item["year"] for item in normalized])
        pooled_values = np.concatenate([item["value"] for item in normalized])
        pooled_fit = choose_fit(pooled_years, pooled_values, args.model)
    if pooled_fit is None:
        raise ValueError("Could not fit a pooled fallback curve from the available data.")

    city_features = {
        item["city"]: item["features"]
        for item in city_feature_rows
        if item["features"] is not None
    }
    row_lookup = {row["city"]: row for row in city_rows}
    scale_lookup: dict[str, float] = {
        city: 1.0 for city in individual_fits
    }
    for city, group in df.groupby("city", sort=True):
        observed = group.dropna(subset=["model_value"])
        years = observed["year"].to_numpy(dtype=float)
        values = observed["model_value"].to_numpy(dtype=float)
        if city in individual_fits or len(observed) == 0:
            continue
        scale, confidence = city_scales(
            city,
            years,
            values,
            pooled_fit,
            city_features,
            scale_lookup,
        )
        scale_lookup[city] = scale
        row = row_lookup[city]
        row["scale"] = scale
        row["confidence"] = confidence
    for city, group in df.groupby("city", sort=True):
        if city in scale_lookup:
            continue
        observed = group.dropna(subset=["model_value"])
        years = observed["year"].to_numpy(dtype=float)
        values = observed["model_value"].to_numpy(dtype=float)
        scale, confidence = city_scales(
            city,
            years,
            values,
            pooled_fit,
            city_features,
            scale_lookup,
        )
        scale_lookup[city] = scale
        row = row_lookup[city]
        row["scale"] = scale
        row["confidence"] = confidence

    fits_by_city = {
        city: {
            "kind": "individual",
            "fit": fit,
            "scale": 1.0,
            "confidence": 0.9,
        }
        for city, fit in individual_fits.items()
    }
    for row in city_rows:
        city = row["city"]
        if city in fits_by_city:
            continue
        fits_by_city[city] = {
            "kind": "pooled" if row["n_obs"] > 0 else "pooled_all_missing",
            "fit": pooled_fit,
            "scale": scale_lookup.get(city, 1.0),
            "confidence": row["confidence"],
        }
        row["selected_model"] = pooled_fit.model
        row["K"] = pooled_fit.params[0] * scale_lookup.get(city, 1.0)
        row["growth"] = pooled_fit.params[1]
        row["inflection"] = pooled_fit.params[2]
        row["rmse"] = pooled_fit.rmse * scale_lookup.get(city, 1.0)
        row["aicc"] = pooled_fit.aicc

    return fits_by_city, {
        "pooled_fit": pooled_fit,
        "city_diagnostics": city_rows,
    }


def fill_model_values(
    df: pd.DataFrame,
    fits_by_city: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    result = df.copy()
    completed_model = np.full(len(result), np.nan)
    lower = np.full(len(result), np.nan)
    upper = np.full(len(result), np.nan)
    method = np.full(len(result), "unfilled", dtype=object)
    confidence = np.full(len(result), np.nan)
    for city, indices in result.groupby("city").groups.items():
        block = result.loc[indices].sort_values("year")
        info = fits_by_city[city]
        fit: CurveFit = info["fit"]
        years = block["year"].to_numpy(dtype=float)
        predicted = fit.predict(years) * info["scale"]
        low, high = fit.interval(years)
        low = low * info["scale"]
        high = high * info["scale"]
        observed = block["model_value"].notna().to_numpy()
        completed_model[block.index.to_numpy()] = np.where(
            observed, block["model_value"].to_numpy(), predicted
        )
        lower[block.index.to_numpy()] = np.where(observed, np.nan, low)
        upper[block.index.to_numpy()] = np.where(observed, np.nan, high)
        if info["kind"] == "individual":
            method[block.index.to_numpy()] = np.where(
                observed, "observed", f"city_{fit.model}"
            )
        elif info["kind"] == "pooled":
            method[block.index.to_numpy()] = np.where(
                observed, "observed", f"pooled_{fit.model}"
            )
        else:
            method[block.index.to_numpy()] = np.where(
                observed, "observed", f"pooled_all_missing_{fit.model}"
            )
        confidence[block.index.to_numpy()] = np.where(
            observed, 1.0, info["confidence"]
        )
    result["completed_model"] = completed_model
    result["lower_95"] = lower
    result["upper_95"] = upper
    result["method"] = method
    result["confidence"] = confidence
    return result


def model_values_to_original(
    df: pd.DataFrame, basis: str, value_unit: str
) -> pd.Series:
    if value_unit == "per_1000" or basis == "raw":
        return df["completed_model"]
    population = df["population"].where(df["population"] > 0)
    return df["completed_model"] * population / 1000.0


def model_column_to_original(
    df: pd.DataFrame, column: str, basis: str, value_unit: str
) -> pd.Series:
    if value_unit == "per_1000" or basis == "raw":
        return df[column]
    population = df["population"].where(df["population"] > 0)
    return df[column] * population / 1000.0


def calibrate_to_totals(
    df: pd.DataFrame,
    totals: pd.DataFrame,
    basis: str,
    value_unit: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    result = df.copy()
    preliminary = model_values_to_original(result, basis, value_unit)
    result["preliminary_value"] = preliminary
    result["lower_95"] = model_column_to_original(
        result, "lower_95", basis, value_unit
    )
    result["upper_95"] = model_column_to_original(
        result, "upper_95", basis, value_unit
    )
    result["completed_value"] = result["vehicles"].where(
        result["vehicles"].notna(), preliminary
    )
    result["total_calibrated"] = False
    diagnostics: list[dict[str, Any]] = []
    if totals.empty:
        return result, diagnostics

    total_map = dict(zip(totals["year"], totals["total"]))
    for year, indices in result.groupby("year").groups.items():
        target = total_map.get(float(year), total_map.get(int(year), np.nan))
        if not np.isfinite(target):
            continue
        block = result.loc[indices]
        observed_mask = block["vehicles"].notna()
        missing_mask = ~observed_mask
        observed_sum = float(block.loc[observed_mask, "vehicles"].sum())
        missing_target = float(target) - observed_sum
        missing_indices = block.index[missing_mask]
        diagnostic = {
            "year": int(year),
            "target_total": float(target),
            "observed_sum": observed_sum,
            "n_missing": int(missing_mask.sum()),
            "missing_target": missing_target,
            "status": "no_missing_rows",
        }
        if missing_indices.empty:
            difference = float(target) - observed_sum
            diagnostic["status"] = (
                "observed_total_matches"
                if abs(difference) <= max(1e-6, abs(target) * 1e-8)
                else "observed_total_mismatch"
            )
            diagnostics.append(diagnostic)
            continue
        tolerance = max(1e-6, abs(float(target)) * 1e-8)
        if missing_target < -tolerance:
            diagnostic["status"] = "inconsistent_observed_exceeds_total"
            diagnostics.append(diagnostic)
            continue
        prelim = block.loc[missing_indices, "preliminary_value"].astype(float)
        weights = prelim.where(prelim > 0)
        if "population" in block.columns:
            weights = weights.fillna(block.loc[missing_indices, "population"])
        weights = weights.where(weights > 0).fillna(1.0)
        weight_sum = float(weights.sum())
        if weight_sum <= 0:
            weights = pd.Series(1.0, index=missing_indices)
            weight_sum = float(len(missing_indices))
        allocation = weights / weight_sum * missing_target
        result.loc[missing_indices, "completed_value"] = allocation
        result.loc[missing_indices, "total_calibrated"] = True
        result.loc[missing_indices, "method"] = (
            result.loc[missing_indices, "method"].astype(str) + "+total"
        )
        diagnostic["status"] = "calibrated"
        diagnostic["allocated_sum"] = float(allocation.sum())
        diagnostics.append(diagnostic)
    return result, diagnostics


def refresh_rates(
    df: pd.DataFrame, basis: str, value_unit: str
) -> pd.DataFrame:
    result = df.copy()
    if value_unit == "per_1000":
        result["completed_per_1000"] = result["completed_value"]
    elif (
        basis == "rate"
        and value_unit == "counts"
        and "population" in result.columns
    ):
        valid = result["population"].notna() & (result["population"] > 0)
        result["completed_per_1000"] = np.where(
            valid,
            1000.0 * result["completed_value"] / result["population"],
            np.nan,
        )
    else:
        result["completed_per_1000"] = np.nan
    return result


def select_output_encoding(
    encoding_info: dict[str, Any], requested: str | None
) -> str:
    if requested:
        try:
            return normalize_encoding_name(requested)
        except LookupError as exc:
            raise ValueError(f"Unknown output encoding: {requested}") from exc
    source = str(encoding_info.get("encoding", "utf-8"))
    if source in {"unicode-xlsx", "ascii", "utf-8", "utf-8-sig"}:
        return "utf-8-sig"
    if source.startswith(("utf-16", "utf-32")):
        return "utf-8-sig"
    return source


def write_missingness(panel: pd.DataFrame, out: Path, encoding: str) -> pd.DataFrame:
    missing = (
        panel.assign(missing=panel["vehicles"].isna())
        .groupby("year", as_index=False)
        .agg(
            n_cities=("city", "size"),
            n_observed=("missing", lambda values: int((~values).sum())),
            n_missing=("missing", "sum"),
        )
    )
    missing.to_csv(out / "missingness_by_year.csv", index=False, encoding=encoding)
    return missing


def write_figures(
    completed: pd.DataFrame,
    city_diagnostics: pd.DataFrame,
    calibration: pd.DataFrame,
    out: Path,
) -> list[str]:
    figure_paths: list[str] = []
    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import font_manager
        import matplotlib.pyplot as plt
    except ImportError:
        return figure_paths

    preferred_fonts = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Noto Sans CJK JP",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    for family in preferred_fonts:
        try:
            font_manager.findfont(
                font_manager.FontProperties(family=family),
                fallback_to_default=False,
            )
        except ValueError:
            continue
        matplotlib.rcParams["font.family"] = family
        break
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    selected = (
        city_diagnostics.sort_values(
            ["fit_type", "n_obs"], ascending=[False, True]
        )["city"]
        .drop_duplicates()
        .head(9)
        .tolist()
    )
    if selected:
        fig, axes = plt.subplots(3, 3, figsize=(12, 9), squeeze=False)
        for axis, city in zip(axes.flat, selected):
            block = completed[completed["city"] == city].sort_values("year")
            observed = block[block["vehicles"].notna()]
            axis.plot(
                block["year"],
                block["completed_value"],
                color="#2f6f73",
                linewidth=1.8,
            )
            axis.scatter(
                observed["year"],
                observed["vehicles"],
                color="#b34d32",
                s=18,
                zorder=3,
            )
            axis.set_title(str(city), fontsize=9)
            axis.grid(alpha=0.2)
        for axis in list(axes.flat)[len(selected) :]:
            axis.axis("off")
        fig.suptitle("Observed values and completed trajectories", fontsize=12)
        fig.tight_layout()
        path = fig_dir / "trajectory_check.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        figure_paths.append(str(path))

    if not calibration.empty:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        ax.plot(
            calibration["year"],
            calibration["target_total"],
            label="Reported total",
            color="#2f6f73",
            linewidth=2,
        )
        ax.plot(
            calibration["year"],
            calibration["completed_sum"],
            label="Completed sum",
            color="#b34d32",
            linewidth=2,
        )
        ax.grid(alpha=0.2)
        ax.legend()
        ax.set_title("Total reconciliation")
        path = fig_dir / "total_reconciliation.png"
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        figure_paths.append(str(path))
    return figure_paths


def write_excel_workbook(
    path: Path,
    long_output: pd.DataFrame,
    wide: pd.DataFrame,
    city_diagnostics: pd.DataFrame,
    missingness: pd.DataFrame,
    calibration: pd.DataFrame,
) -> bool:
    try:
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            long_output.to_excel(writer, sheet_name="completed_long", index=False)
            wide.to_excel(writer, sheet_name="completed_wide", index=True)
            city_diagnostics.to_excel(writer, sheet_name="city_parameters", index=False)
            missingness.to_excel(writer, sheet_name="missingness", index=False)
            if not calibration.empty:
                calibration.to_excel(
                    writer, sheet_name="total_calibration", index=False
                )
            for worksheet in writer.book.worksheets:
                worksheet.freeze_panes = "A2"
                worksheet.auto_filter.ref = worksheet.dimensions
    except ImportError:
        return False
    return True


def build_report(
    args: argparse.Namespace,
    panel_format: str,
    basis: str,
    encoding_info: dict[str, Any],
    output_encoding: str,
    workbook_created: bool,
    panel: pd.DataFrame,
    completed: pd.DataFrame,
    city_diagnostics: pd.DataFrame,
    pooled_fit: CurveFit,
    missingness: pd.DataFrame,
    calibration_rows: list[dict[str, Any]],
    warnings: list[str],
) -> str:
    total_cities = int(panel["city"].nunique())
    observed_cells = int(panel["vehicles"].notna().sum())
    missing_cells = int(panel["vehicles"].isna().sum())
    city_counts = city_diagnostics["fit_type"].value_counts().to_dict()
    model_counts = city_diagnostics["selected_model"].fillna("none").value_counts()
    calibrated_cells = int(
        completed.get("total_calibrated", pd.Series(dtype=bool)).fillna(False).sum()
    )
    lines = [
        "# 城市汽车保有量面板补全报告",
        "",
        "## 输入概览",
        "",
        f"- 输入文件：`{Path(args.input).resolve()}`",
        f"- 输入布局：{panel_format}",
        f"- 输入文本编码：{encoding_info.get('encoding')} "
        f"（{encoding_info.get('method')}）",
        f"- CSV 输出编码：{output_encoding}",
        f"- Unicode Excel 工作簿：{'已生成' if workbook_created else '未生成'}",
        f"- 城市数：{total_cities}",
        f"- 年份范围：{int(panel['year'].min())}–{int(panel['year'].max())}",
        f"- 观测像元：{observed_cells}",
        f"- 缺失像元：{missing_cells}",
        f"- 建模口径：{basis}",
        f"- 原始值单位：{args.value_unit}",
        f"- 请求模型：{args.model}",
        "",
        "## 模型结果",
        "",
        f"- 单城市曲线：{city_counts.get('city_specific', 0)} 个城市",
        f"- 借用同类城市曲线：{city_counts.get('pooled', 0)} 个城市",
        f"- 所选模型：{dict(model_counts)}",
        f"- 汇总备用曲线：{pooled_fit.model}",
        f"- 汇总备用曲线 RMSE：{pooled_fit.rmse:.6g}",
        f"- 按总量校准的像元：{calibrated_cells}",
        "",
    ]
    dense_missing = city_diagnostics[
        city_diagnostics["n_obs"] == 0
    ]["city"].astype(str).tolist()
    sparse_missing = city_diagnostics[
        (city_diagnostics["n_obs"] > 0)
        & (city_diagnostics["n_obs"] < args.min_obs)
    ]["city"].astype(str).tolist()
    if dense_missing:
        preview = ", ".join(dense_missing[:20])
        lines.extend(
            [
                "## 低可识别城市",
                "",
                f"以下 {len(dense_missing)} 个城市没有任何实测值，补全结果属于低置信度外推：",
                "",
                preview + (" ..." if len(dense_missing) > 20 else ""),
                "",
            ]
        )
    if sparse_missing:
        preview = ", ".join(sparse_missing[:20])
        lines.extend(
            [
                "## 稀疏观测城市",
                "",
                f"以下 {len(sparse_missing)} 个城市观测点少于 {args.min_obs} 个，主要借用汇总曲线：",
                "",
                preview + (" ..." if len(sparse_missing) > 20 else ""),
                "",
            ]
        )
    if calibration_rows:
        calibration = pd.DataFrame(calibration_rows)
        lines.extend(
            [
                "## 总量校准",
                "",
                "```text",
                calibration.to_string(index=False),
                "```",
                "",
            ]
        )
    if warnings:
        lines.extend(["## 警告", ""])
        lines.extend(f"- {warning}" for warning in warnings)
        lines.append("")
    lines.extend(
        [
            "## 解释与限制",
            "",
            "- `observed` 表示原始实测值，脚本不会修改该值。",
            "- `city_logistic`、`city_gompertz` 表示由该城市自身数据拟合。",
            "- `pooled_*` 表示观测太少，借用其他城市的曲线形状。",
            "- `pooled_all_missing_*` 表示该城市没有实测值，只能依据人口、GDP 近邻或总体水平补全。",
            "- `+total` 表示该值随后按年度总量校准。",
            "- 95% 区间为近似区间，未完整传播模型选择和总量校准的不确定性。",
            "- 补全值不应被标注为官方观测值；正式研究应报告补全比例和敏感性分析。",
            "- CSV 的编码影响文本读写；字体只影响 Excel/PDF/图片的实际显示。"
            " 挪威语 æ、ø、å 属于 Unicode 字符，只要编码和字体支持即可正确保存。",
            "",
        ]
    )
    return "\n".join(lines)


def run(args: argparse.Namespace) -> dict[str, Any]:
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    panel, panel_format, input_encoding = load_panel(args)
    output_encoding = select_output_encoding(
        input_encoding, args.output_encoding
    )
    panel = aggregate_duplicates(panel, args.aggregate_duplicates)
    if panel.empty:
        raise ValueError("No valid city-year rows remain after cleaning.")
    if panel["vehicles"].dropna().lt(0).any():
        raise ValueError("Vehicle values must be non-negative.")
    totals, total_encoding_infos = extract_totals(panel, args)
    start, end = determine_year_range(panel, args)
    panel = complete_grid(panel, start, end)
    panel, population_filled = interpolate_population(panel)
    basis = determine_basis(panel, args)
    panel, warnings = add_model_basis(panel, basis, args.value_unit)
    if population_filled:
        warnings.append(
            f"{population_filled} missing population cells were linearly interpolated within city."
        )
    fits_by_city, fit_info = fit_all_cities(panel, args)
    completed = fill_model_values(panel, fits_by_city)
    completed, calibration_rows = calibrate_to_totals(
        completed, totals, basis, args.value_unit
    )
    completed = refresh_rates(completed, basis, args.value_unit)
    completed["is_filled"] = completed["vehicles"].isna()
    notes = np.where(
        completed["method"].astype(str).str.contains("total"),
        "Total-calibrated value",
        "",
    )
    completed["notes"] = notes
    missingness = write_missingness(completed, out, output_encoding)

    long_columns = [
        "city",
        "year",
        "vehicles",
        "completed_value",
        "completed_per_1000",
        "population",
        "gdp",
        "total",
        "is_filled",
        "method",
        "confidence",
        "lower_95",
        "upper_95",
        "total_calibrated",
        "notes",
    ]
    long_output = completed[
        [column for column in long_columns if column in completed.columns]
    ].rename(columns={"vehicles": "observed_value"})
    long_output.to_csv(
        out / "completed_panel_long.csv",
        index=False,
        encoding=output_encoding,
    )
    wide = completed.pivot(index="city", columns="year", values="completed_value")
    wide.to_csv(out / "completed_panel_wide.csv", encoding=output_encoding)

    city_diagnostics = pd.DataFrame(fit_info["city_diagnostics"]).sort_values(
        "city"
    )
    city_diagnostics.to_csv(
        out / "city_parameters.csv", index=False, encoding=output_encoding
    )
    calibration = pd.DataFrame(calibration_rows)
    if not calibration.empty:
        completed_sum = (
            completed.groupby("year", as_index=False)["completed_value"]
            .sum()
            .rename(columns={"completed_value": "completed_sum"})
        )
        calibration = calibration.merge(completed_sum, on="year", how="left")
        calibration.to_csv(
            out / "total_calibration.csv", index=False, encoding=output_encoding
        )
    workbook_path = out / "completion_workbook.xlsx"
    workbook_created = False
    if not args.no_excel:
        workbook_created = write_excel_workbook(
            workbook_path,
            long_output,
            wide,
            city_diagnostics,
            missingness,
            calibration,
        )
        if not workbook_created:
            warnings.append(
                "XLSX output was skipped because openpyxl is not installed."
            )
    figure_paths = []
    if not args.no_figures:
        figure_paths = write_figures(
            completed, city_diagnostics, calibration, out
        )

    report = build_report(
        args,
        panel_format,
        basis,
        input_encoding,
        output_encoding,
        workbook_created,
        panel,
        completed,
        city_diagnostics,
        fit_info["pooled_fit"],
        missingness,
        calibration_rows,
        warnings,
    )
    (out / "completion_report.md").write_text(report, encoding="utf-8")
    summary = {
        "input": str(Path(args.input).resolve()),
        "output": str(out),
        "panel_format": panel_format,
        "input_encoding": input_encoding,
        "total_file_encodings": total_encoding_infos,
        "output_encoding": output_encoding,
        "excel_workbook": str(workbook_path) if workbook_created else None,
        "basis": basis,
        "requested_model": args.model,
        "year_start": start,
        "year_end": end,
        "n_cities": int(panel["city"].nunique()),
        "n_cells": int(len(panel)),
        "n_observed": int(panel["vehicles"].notna().sum()),
        "n_missing": int(panel["vehicles"].isna().sum()),
        "n_city_specific": int(
            (city_diagnostics["fit_type"] == "city_specific").sum()
        ),
        "n_pooled": int((city_diagnostics["fit_type"] == "pooled").sum()),
        "n_all_missing_cities": int((city_diagnostics["n_obs"] == 0).sum()),
        "n_total_calibrated": int(completed["total_calibrated"].sum()),
        "warnings": warnings,
        "figures": figure_paths,
    }
    (out / "completion_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    args = parse_args()
    try:
        summary = run(args)
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
