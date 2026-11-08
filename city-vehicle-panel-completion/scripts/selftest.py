#!/usr/bin/env python3
"""Run an end-to-end smoke test for the panel completion script."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    script = Path(__file__).with_name("complete_vehicle_panel.py")
    years = np.arange(2010, 2016)
    rows = []
    missing = {
        "A": set(),
        "B": {2012},
        "C": {2011, 2012, 2013, 2014},
        "D": set(years.tolist()),
        "E": {2010, 2015},
        "F": set(),
    }
    truth = {}
    for index, city in enumerate("ABCDEF"):
        saturation = 1800 + index * 700
        growth = 0.55 + index * 0.04
        inflection = 2012.0 + index * 0.2
        values = saturation / (1 + np.exp(-growth * (years - inflection)))
        truth[city] = values
    totals = {
        int(year): float(
            sum(truth[city][position] for city in truth)
        )
        for position, year in enumerate(years)
    }
    for index, city in enumerate("ABCDEF"):
        for position, year in enumerate(years):
            value = np.nan if year in missing[city] else truth[city][position]
            rows.append(
                {
                    "city": city,
                    "year": int(year),
                    "vehicles": value,
                    "population": 2_000_000 + index * 300_000,
                    "gdp": 60_000 + index * 8_000,
                    "total": totals[int(year)],
                }
            )
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        input_path = root / "panel.csv"
        output_path = root / "out"
        pd.DataFrame(rows).to_csv(input_path, index=False, encoding="utf-8-sig")
        result = subprocess.run(
            [
                sys.executable,
                str(script),
                "--input",
                str(input_path),
                "--out",
                str(output_path),
                "--population-col",
                "population",
                "--gdp-col",
                "gdp",
                "--total-col",
                "total",
                "--model",
                "auto",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + "\n" + result.stderr)
        completed = pd.read_csv(
            output_path / "completed_panel_long.csv", encoding="utf-8-sig"
        )
        if completed["completed_value"].isna().any():
            raise AssertionError("Smoke test left unfilled cells.")
        observed = completed.dropna(subset=["observed_value"])
        if not np.allclose(
            observed["observed_value"], observed["completed_value"], rtol=0, atol=1e-9
        ):
            raise AssertionError("Observed values were modified.")
        annual = completed.groupby("year")["completed_value"].sum().to_dict()
        for year, target in totals.items():
            if not np.isclose(annual[year], target, rtol=1e-8, atol=1e-6):
                raise AssertionError(f"Total calibration failed in {year}.")
        city_d = completed[completed["city"] == "D"]
        if not city_d["method"].astype(str).str.contains("pooled_all_missing").all():
            raise AssertionError("All-missing city was not labelled as pooled fallback.")
        summary = json.loads(
            (output_path / "completion_summary.json").read_text(encoding="utf-8")
        )
        if not summary["figures"]:
            raise AssertionError("Expected diagnostic figures to be generated.")

        wide = pd.DataFrame(
            [
                {
                    "city": city,
                    "population": 2_000_000 + index * 300_000,
                    "gdp": 60_000 + index * 8_000,
                    **{
                        str(year): (
                            np.nan if year in missing[city] else truth[city][position]
                        )
                        for position, year in enumerate(years)
                    },
                }
                for index, city in enumerate("ABCDEF")
            ]
        )
        wide_input = root / "panel_wide.xlsx"
        wide_output = root / "out_wide"
        wide.to_excel(wide_input, index=False)
        wide_result = subprocess.run(
            [
                sys.executable,
                str(script),
                "--input",
                str(wide_input),
                "--out",
                str(wide_output),
                "--format",
                "wide",
                "--population-col",
                "population",
                "--gdp-col",
                "gdp",
                "--model",
                "gompertz",
                "--no-figures",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if wide_result.returncode != 0:
            raise RuntimeError(wide_result.stdout + "\n" + wide_result.stderr)
        wide_completed = pd.read_csv(
            wide_output / "completed_panel_long.csv", encoding="utf-8-sig"
        )
        if len(wide_completed) != len("ABCDEF") * len(years):
            raise AssertionError("Wide-format input did not produce the full grid.")
        if wide_completed["completed_value"].isna().any():
            raise AssertionError("Wide-format smoke test left unfilled cells.")

        norwegian_cities = ["Bærum", "Tromsø", "Ålesund", "Møre og Romsdal"]
        norwegian_rows = []
        for city_index, city in enumerate(norwegian_cities):
            for year_index, year in enumerate(years):
                value = (
                    np.nan
                    if (year_index == 2 and city_index % 2 == 0)
                    else 1200 + city_index * 300 + year_index * 90
                )
                norwegian_rows.append(
                    {
                        "city": city,
                        "year": int(year),
                        "vehicles": value,
                        "population": 1_500_000 + city_index * 200_000,
                    }
                )
        norwegian_input = root / "norwegian_cp1252.csv"
        norwegian_output = root / "out_norwegian"
        pd.DataFrame(norwegian_rows).to_csv(
            norwegian_input, index=False, encoding="cp1252"
        )
        norwegian_result = subprocess.run(
            [
                sys.executable,
                str(script),
                "--input",
                str(norwegian_input),
                "--out",
                str(norwegian_output),
                "--population-col",
                "population",
                "--no-figures",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if norwegian_result.returncode != 0:
            raise RuntimeError(norwegian_result.stdout + "\n" + norwegian_result.stderr)
        norwegian_summary = json.loads(
            (norwegian_output / "completion_summary.json").read_text(encoding="utf-8")
        )
        if norwegian_summary["input_encoding"]["encoding"] not in {
            "cp1252",
            "iso8859-1",
        }:
            raise AssertionError(
                "CP1252 Norwegian input was not detected: "
                f"{norwegian_summary['input_encoding']}"
            )
        detected_output_encoding = norwegian_summary["output_encoding"]
        norwegian_completed = pd.read_csv(
            norwegian_output / "completed_panel_long.csv",
            encoding=detected_output_encoding,
        )
        if set(norwegian_completed["city"]) != set(norwegian_cities):
            raise AssertionError("Norwegian characters were corrupted.")
        workbook = pd.read_excel(
            norwegian_output / "completion_workbook.xlsx",
            sheet_name="completed_long",
        )
        if set(workbook["city"]) != set(norwegian_cities):
            raise AssertionError("Norwegian characters were corrupted in XLSX.")
        print(
            json.dumps(
                {
                    "status": "ok",
                    "cities": summary["n_cities"],
                    "missing_filled": summary["n_missing"],
                    "total_calibrated": summary["n_total_calibrated"],
                    "wide_format": "ok",
                    "cp1252_norwegian": "ok",
                },
                ensure_ascii=False,
            )
        )


if __name__ == "__main__":
    main()
