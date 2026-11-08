#!/usr/bin/env python3
"""Two-parameter Weibull vehicle lifetime calculations.

CSV text is never silently decoded or rewritten as UTF-8. If the source
encoding cannot be identified reliably, pass --source-encoding explicitly.
"""

from __future__ import annotations

import argparse
import codecs
import csv
import datetime as dt
import io
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


BOM_ENCODINGS = (
    (codecs.BOM_UTF32_LE, "utf-32"),
    (codecs.BOM_UTF32_BE, "utf-32"),
    (codecs.BOM_UTF8, "utf-8-sig"),
    (codecs.BOM_UTF16_LE, "utf-16"),
    (codecs.BOM_UTF16_BE, "utf-16"),
)


@dataclass(frozen=True)
class WeibullLifetime:
    """Two-parameter Weibull lifetime model with eta in years and beta dimensionless."""

    eta: float
    beta: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.eta) or self.eta <= 0:
            raise ValueError("eta must be a finite positive number")
        if not math.isfinite(self.beta) or self.beta <= 0:
            raise ValueError("beta must be a finite positive number")

    def survival(self, age: float) -> float:
        if age < 0:
            raise ValueError("age must be non-negative")
        return math.exp(-((age / self.eta) ** self.beta))

    def cdf(self, age: float) -> float:
        return 1.0 - self.survival(age)

    def pdf(self, age: float) -> float:
        if age < 0:
            raise ValueError("age must be non-negative")
        if age == 0 and self.beta < 1:
            return math.inf
        return (
            (self.beta / self.eta)
            * ((age / self.eta) ** (self.beta - 1.0))
            * math.exp(-((age / self.eta) ** self.beta))
        )

    def hazard(self, age: float) -> float:
        if age < 0:
            raise ValueError("age must be non-negative")
        if age == 0 and self.beta < 1:
            return math.inf
        return (self.beta / self.eta) * ((age / self.eta) ** (self.beta - 1.0))

    def mean(self) -> float:
        return self.eta * math.gamma(1.0 + 1.0 / self.beta)

    def variance(self) -> float:
        g1 = math.gamma(1.0 + 1.0 / self.beta)
        g2 = math.gamma(1.0 + 2.0 / self.beta)
        return self.eta * self.eta * (g2 - g1 * g1)

    def median(self) -> float:
        return self.eta * (math.log(2.0) ** (1.0 / self.beta))

    def quantile(self, probability: float) -> float:
        if not 0.0 < probability < 1.0:
            raise ValueError("probability must be strictly between 0 and 1")
        return self.eta * ((-math.log(1.0 - probability)) ** (1.0 / self.beta))

    def scrap_probabilities(self, max_age: int) -> list[float]:
        if max_age < 1:
            raise ValueError("max_age must be at least 1")
        survival_prev = 1.0
        probabilities: list[float] = []
        for age in range(1, max_age + 1):
            survival_now = self.survival(float(age))
            probabilities.append(max(0.0, survival_prev - survival_now))
            survival_prev = survival_now
        return probabilities

    def survival_table(self, max_age: int) -> list[dict[str, float]]:
        if max_age < 0:
            raise ValueError("max_age must be non-negative")
        rows: list[dict[str, float]] = []
        for age in range(max_age + 1):
            rows.append(
                {
                    "age": float(age),
                    "survival": self.survival(float(age)),
                    "cdf": self.cdf(float(age)),
                    "pdf": self.pdf(float(age)),
                    "hazard": self.hazard(float(age)),
                }
            )
        return rows


def eta_from_mean(mean: float, beta: float) -> float:
    if not math.isfinite(mean) or mean <= 0:
        raise ValueError("mean must be a finite positive number")
    WeibullLifetime(1.0, beta)
    return mean / math.gamma(1.0 + 1.0 / beta)


def eta_from_median(median: float, beta: float) -> float:
    if not math.isfinite(median) or median <= 0:
        raise ValueError("median must be a finite positive number")
    WeibullLifetime(1.0, beta)
    return median / (math.log(2.0) ** (1.0 / beta))


def eta_from_cdf(age: float, cdf: float, beta: float) -> float:
    if not math.isfinite(age) or age <= 0:
        raise ValueError("age must be a finite positive number")
    if not 0.0 < cdf < 1.0:
        raise ValueError("cdf must be strictly between 0 and 1")
    WeibullLifetime(1.0, beta)
    return age / ((-math.log(1.0 - cdf)) ** (1.0 / beta))


def parameter_summary(eta: float, beta: float) -> dict[str, float]:
    model = WeibullLifetime(eta, beta)
    return {
        "eta_years": model.eta,
        "beta": model.beta,
        "mean_lifetime_years": model.mean(),
        "variance": model.variance(),
        "std_dev_years": math.sqrt(model.variance()),
        "median_lifetime_years": model.median(),
        "p10_lifetime_years": model.quantile(0.10),
        "p90_lifetime_years": model.quantile(0.90),
    }


def iter_survival_rows(eta: float, beta: float, max_age: int) -> Iterable[dict[str, float]]:
    return iter(WeibullLifetime(eta, beta).survival_table(max_age))


def normalize_encoding(label: str) -> str:
    """Return a strict codec name accepted by Python."""

    try:
        return codecs.lookup(label).name
    except LookupError as exc:
        raise ValueError(f"unknown encoding: {label!r}") from exc


def detect_source_encoding(data: bytes) -> str:
    """Detect only encodings that are unambiguous without guessing.

    BOMs are exact and ASCII is exact. Non-ASCII bytes without a BOM are
    ambiguous between UTF-8 and legacy single-byte encodings, so they are
    never guessed. Pass --source-encoding explicitly for those files.
    """

    for bom, encoding in BOM_ENCODINGS:
        if data.startswith(bom):
            return encoding

    if not data or all(byte < 0x80 for byte in data):
        return "ascii"

    raise ValueError(
        "Could not reliably detect the source encoding. Non-ASCII bytes "
        "without a BOM may be UTF-8 or a legacy encoding. Pass "
        "--source-encoding with the original encoding (for example CP1252, "
        "Windows-1252, ISO-8859-1, GBK, Shift_JIS, or UTF-8). No decoding "
        "or replacement was attempted."
    )


def read_source_text(path: Path, source_encoding: str | None) -> tuple[str, str]:
    """Read text with strict decoding and return text plus source encoding."""

    raw = path.read_bytes()
    encoding = (
        normalize_encoding(source_encoding)
        if source_encoding
        else detect_source_encoding(raw)
    )
    try:
        text = raw.decode(encoding, errors="strict")
    except UnicodeDecodeError as exc:
        raise ValueError(
            f"Could not decode {path} using {encoding!r} without replacement. "
            "Check --source-encoding; no fallback decoding was attempted."
        ) from exc

    # A BOM is not part of the first logical field.
    if text.startswith("\ufeff"):
        text = text[1:]
    return text, encoding


def _backup_path(path: Path) -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    candidate = path.with_name(f"{path.name}.bak-{stamp}")
    counter = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.name}.bak-{stamp}-{counter}")
        counter += 1
    return candidate


def backup_file(path: Path) -> Path:
    backup = _backup_path(path)
    shutil.copy2(path, backup)
    return backup


def read_cohort_csv(
    path: Path,
    year_col: str,
    sales_col: str,
    source_encoding: str | None = None,
) -> tuple[dict[int, float], str]:
    text, encoding = read_source_text(path, source_encoding)
    cohorts: dict[int, float] = {}
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames is None:
        raise ValueError("cohort CSV has no header")

    missing = [name for name in (year_col, sales_col) if name not in reader.fieldnames]
    if missing:
        raise ValueError(f"cohort CSV is missing required columns: {missing}")

    for line_no, row in enumerate(reader, start=2):
        try:
            year = int(float(row[year_col]))
            sales = float(row[sales_col])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid row {line_no}: {row}") from exc
        if sales < 0:
            raise ValueError(f"negative sales at row {line_no}")
        cohorts[year] = cohorts.get(year, 0.0) + sales
    return cohorts, encoding


def cohort_outflows(
    cohorts: dict[int, float],
    eta: float,
    beta: float,
    start_year: int,
    end_year: int,
) -> list[dict[str, float]]:
    if start_year > end_year:
        raise ValueError("start_year must not exceed end_year")
    if not cohorts:
        raise ValueError("cohort CSV contains no usable sales rows")

    model = WeibullLifetime(eta, beta)
    survival_by_age = {
        age: model.survival(float(age)) for age in range(0, end_year - min(cohorts) + 2)
    }
    rows: list[dict[str, float]] = []
    for year in range(start_year, end_year + 1):
        total = 0.0
        for cohort_year, sales in cohorts.items():
            age = year - cohort_year
            if age < 1:
                continue
            previous = survival_by_age.get(age - 1, model.survival(float(age - 1)))
            current = survival_by_age.get(age, model.survival(float(age)))
            total += sales * max(0.0, previous - current)
        rows.append({"year": float(year), "outflow_units": total})
    return rows


def write_rows(
    path: Path | None,
    rows: list[dict[str, float]],
    *,
    output_encoding: str | None = None,
    source_encoding: str | None = None,
    in_place: bool = False,
) -> None:
    """Write rows without silently changing encoding.

    If output_encoding is omitted and a source encoding is available, the
    output uses the source encoding. If there is no source encoding, a new
    output file requires an explicit --output-encoding.
    """

    if path is None:
        for row in rows:
            print(json.dumps(row, ensure_ascii=False))
        return

    if not rows:
        raise ValueError("no rows to write")

    if output_encoding:
        target_encoding = normalize_encoding(output_encoding)
        explicit_conversion = bool(
            source_encoding
            and normalize_encoding(source_encoding) != target_encoding
        )
    elif source_encoding:
        target_encoding = normalize_encoding(source_encoding)
        explicit_conversion = False
    else:
        raise ValueError(
            "A file output requires --output-encoding when no source encoding "
            "is available. The script will not assume UTF-8."
        )

    backup = None
    if path.exists() and (explicit_conversion or in_place):
        backup = backup_file(path)

    if explicit_conversion:
        message = f"Encoding conversion: {source_encoding} -> {output_encoding}"
        if backup is not None:
            message += f"; backup: {backup}"
        print(message, file=sys.stderr)

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open(
            "w",
            encoding=target_encoding,
            errors="strict",
            newline="",
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    except UnicodeEncodeError as exc:
        raise ValueError(
            f"Output cannot be encoded as {target_encoding!r} without replacement. "
            "Choose a different --output-encoding explicitly or correct the input data."
        ) from exc


def models_from_args(args: argparse.Namespace) -> WeibullLifetime:
    model = WeibullLifetime(args.eta, args.beta)
    if args.eta_multiplier is not None:
        model = WeibullLifetime(model.eta * args.eta_multiplier, model.beta)
    return model


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    def add_parameters(target: argparse.ArgumentParser) -> None:
        target.add_argument("--eta", type=float, required=True, help="Weibull scale parameter in years")
        target.add_argument("--beta", type=float, required=True, help="Weibull shape parameter")
        target.add_argument("--eta-multiplier", type=float, default=None, help="Optional scenario multiplier applied to eta")

    def add_output_encoding(target: argparse.ArgumentParser) -> None:
        target.add_argument(
            "--output-encoding",
            "--output_encoding",
            dest="output_encoding",
            default=None,
            help="Explicit output encoding. Omit to preserve the source encoding when available.",
        )

    mean_parser = sub.add_parser("mean", help="Compute mean lifetime")
    add_parameters(mean_parser)
    mean_parser.add_argument("--json", action="store_true")

    eta_parser = sub.add_parser("eta", help="Recover eta from mean lifetime")
    eta_parser.add_argument("--mean", type=float, required=True)
    eta_parser.add_argument("--beta", type=float, required=True)
    eta_parser.add_argument("--json", action="store_true")

    detect_parser = sub.add_parser(
        "detect-encoding",
        aliases=["detect_encoding", "detect"],
        help="Report a reliably detected source encoding without guessing",
    )
    detect_parser.add_argument("--input", type=Path, required=True)

    survival_parser = sub.add_parser("survival", help="Print survival table")
    add_parameters(survival_parser)
    survival_parser.add_argument("--max-age", type=int, default=30)
    survival_parser.add_argument("--out", type=Path, default=None)
    add_output_encoding(survival_parser)

    scrap_parser = sub.add_parser("scrap", help="Print annual scrappage probabilities")
    add_parameters(scrap_parser)
    scrap_parser.add_argument("--max-age", type=int, default=30)
    scrap_parser.add_argument("--out", type=Path, default=None)
    add_output_encoding(scrap_parser)

    table_parser = sub.add_parser("table", help="Compute annual outflows from cohort sales")
    add_parameters(table_parser)
    table_parser.add_argument("--input", type=Path, required=True)
    table_parser.add_argument("--year-col", default="year")
    table_parser.add_argument("--sales-col", default="sales")
    table_parser.add_argument("--start-year", type=int, required=True)
    table_parser.add_argument("--end-year", type=int, required=True)
    table_parser.add_argument("--out", type=Path, default=None)
    table_parser.add_argument(
        "--source-encoding",
        "--source_encoding",
        dest="source_encoding",
        default=None,
        help="Source encoding, for example CP1252, Windows-1252, ISO-8859-1, GBK, Shift_JIS, or UTF-8.",
    )
    add_output_encoding(table_parser)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "mean":
            summary = parameter_summary(args.eta, args.beta)
            if args.eta_multiplier is not None:
                summary = parameter_summary(args.eta * args.eta_multiplier, args.beta)
            if args.json:
                print(json.dumps(summary, ensure_ascii=False, indent=2))
            else:
                print(summary["mean_lifetime_years"])
            return 0
        if args.command == "eta":
            eta = eta_from_mean(args.mean, args.beta)
            payload = {"eta_years": eta, "beta": args.beta, "mean_lifetime_years": args.mean}
            if args.json:
                print(json.dumps(payload, ensure_ascii=False, indent=2))
            else:
                print(eta)
            return 0
        if args.command in {"detect-encoding", "detect_encoding", "detect"}:
            raw = args.input.read_bytes()
            encoding = detect_source_encoding(raw)
            print(json.dumps({"input": str(args.input), "encoding": encoding}, ensure_ascii=False, indent=2))
            return 0
        if args.command == "survival":
            model = models_from_args(args)
            rows = list(iter_survival_rows(model.eta, model.beta, args.max_age))
            write_rows(args.out, rows, output_encoding=args.output_encoding)
            return 0
        if args.command == "scrap":
            model = models_from_args(args)
            probabilities = model.scrap_probabilities(args.max_age)
            rows = [
                {"age": float(age), "scrap_probability": probability}
                for age, probability in enumerate(probabilities, start=1)
            ]
            write_rows(args.out, rows, output_encoding=args.output_encoding)
            return 0
        if args.command == "table":
            model = models_from_args(args)
            cohorts, source_encoding = read_cohort_csv(
                args.input,
                args.year_col,
                args.sales_col,
                args.source_encoding,
            )
            rows = cohort_outflows(cohorts, model.eta, model.beta, args.start_year, args.end_year)
            write_rows(
                args.out,
                rows,
                output_encoding=args.output_encoding,
                source_encoding=source_encoding,
                in_place=bool(args.out and args.out.resolve() == args.input.resolve()),
            )
            return 0
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())