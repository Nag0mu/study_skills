---
name: city-vehicle-panel-completion
description: Complete unbalanced city-year or region-year vehicle ownership panels with logistic or Gompertz saturation curves, pooled fallbacks for sparse cities, and optional yearly-total calibration. Use when vehicle stock/ownership data have missing years, partial city coverage, unbalanced panels, or known regional totals that must be reconciled.
metadata:
  short-description: Complete missing city vehicle ownership panels
---

# City Vehicle Panel Completion

Use this workflow when the user supplies vehicle ownership or vehicle stock data for cities, provinces, or countries over time and asks to fill missing years, construct a balanced panel, interpolate gaps, or produce forecast-ready inputs.

The central rule is: **preserve observed values exactly and label every filled value by its estimation method.** A completed cell is not an official observation.

## Workflow

1. Inspect the source workbook or table without editing it. Identify:
   - city/region identifier and year
   - vehicle variable and whether it is a count or already per 1,000 people
   - population, GDP per capita, and any yearly regional/national total
   - whether the data are long (`city, year, value`) or wide (city rows, year columns)
   - the source text encoding, especially for city names containing characters such as `æ`, `ø`, `å`, or Chinese characters
2. Resolve ambiguous units before fitting. Ask only when the source cannot establish whether the vehicle variable is total stock/count or a normalized ownership rate.
3. Copy or derive the input into a working file; never overwrite the source data.
4. Run `scripts/complete_vehicle_panel.py` with explicit column mappings. Prefer `--basis auto`; it models vehicles per 1,000 population when population and counts are available.
5. Inspect `completion_report.md`, `completion_summary.json`, `city_parameters.csv`, and `total_calibration.csv` before presenting results.
6. If a yearly total is inconsistent with the observed city values, report the conflict. Do not silently reduce observed values to force reconciliation.
7. Deliver the completed long and wide panels plus the report. State how many cells were filled, how many cities used pooled fallback, and how many values were total-calibrated.

## Input Layouts

Long format:

```text
city,year,vehicles,population,gdp,total
```

Wide format:

```text
city,population,gdp,2008,2009,...,2019
```

Blank cells and tokens such as `NA`, `N/A`, `-`, `--`, `.`, and `..` are treated as missing. Never convert a blank into zero.

## Encoding And Fonts

Text encoding and font rendering are different problems:

- CSV encoding controls whether bytes become the correct Unicode characters.
- XLSX stores Unicode internally; `æ`, `ø`, `å`, and Chinese characters do not depend on the workbook's font for storage.
- Fonts matter only when characters are displayed or rendered to PDF/images. A font lacking a glyph causes a visible box even when the underlying text is correct.

The script defaults to `--encoding auto`. Detection checks UTF BOMs, strict UTF-8, `charset-normalizer`, and strict fallback codecs. It explicitly corrects common Western-European misclassification of CP1252 files as CP1250 when characters such as `æ`, `ø`, and `å` indicate CP1252. If detection is ambiguous, require an explicit codec:

```text
--encoding cp1252
--encoding iso-8859-1
--encoding gb18030
--encoding shift_jis
```

CSV output defaults to a Unicode-safe encoding derived from the input. Override it with `--output-encoding`. The script also writes `completion_workbook.xlsx`, which preserves Unicode independently of CSV settings. Use `--no-excel` only when XLSX output is not needed.

Never use replacement decoding. If a source encoding cannot be established safely, stop and request the original encoding rather than silently turning characters into `�`.

Read [references/input-and-output.md](references/input-and-output.md) when mapping columns, resolving units, or explaining output fields and status codes.

## Run The Completion Script

Long-format example:

```bash
python scripts/complete_vehicle_panel.py \
  --input "city_vehicles.xlsx" \
  --out "city_vehicles_completed" \
  --city-col "城市" \
  --year-col "年份" \
  --value-col "汽车保有量" \
  --population-col "常住人口" \
  --total-col "区域总量" \
  --value-unit counts \
  --basis auto \
  --model auto
```

Wide-format example:

```bash
python scripts/complete_vehicle_panel.py \
  --input "city_vehicles_wide.xlsx" \
  --out "city_vehicles_completed" \
  --format wide \
  --city-col "城市" \
  --population-col "人口" \
  --value-unit counts \
  --model gompertz
```

If totals are in another file, use `--total-file`, `--total-year-col`, and `--total-value-col`.

Run the bundled smoke test after changing the script:

```bash
python scripts/selftest.py
```

## Modeling Rules

- Use a city-specific curve only when the city has enough informative observations. Defaults are at least 5 observed years and a span of at least 3 years.
- Use a pooled curve for sparse cities. The script normalizes city trajectories, fits a common shape, and estimates each city's scale.
- For cities with no observations, use population/GDP-nearest observed cities when available. Mark these values as low confidence.
- Use logistic as the transparent baseline and Gompertz as the main robustness alternative. `--model auto` compares AICc when both converge.
- Treat extrapolation beyond the observed year range as optional and explicit. Do not fill future years unless the user requested forecasts or set `--start-year`/`--end-year` with `--extrapolate`.
- If a known yearly total exists, calibrate only the missing cells. Allocate the target missing total according to preliminary model shares, with population as a fallback weight.
- If `observed_sum > reported_total`, skip calibration for that year and flag it. Never alter observed values.
- Treat 95% intervals as approximate. They do not include uncertainty from model selection, missingness, or total constraints.

Read [references/method.md](references/method.md) when choosing the modeling basis, assessing missingness, defending the method, or writing a thesis/paper methods section.

## Output Contract

The script writes:

- `completed_panel_long.csv`: auditable one-row-per-city-year panel
- `completed_panel_wide.csv`: completed city-by-year matrix
- `city_parameters.csv`: observations, selected model, parameters, RMSE, scale, and confidence
- `missingness_by_year.csv`: observed and missing counts by year
- `total_calibration.csv`: yearly reconciliation diagnostics when totals exist
- `completion_report.md`: methods, diagnostics, warnings, and limitations
- `completion_summary.json`: machine-readable run summary
- `completion_workbook.xlsx`: Unicode-safe workbook containing long, wide, parameter, missingness, and calibration sheets
- `figures/`: trajectory and total-reconciliation checks unless `--no-figures`

The key long-panel fields are:

- `observed_value`: original value; unchanged
- `completed_value`: observed value or model-filled value
- `method`: `observed`, `city_*`, `pooled_*`, `pooled_all_missing_*`, with `+total` when calibrated
- `confidence`: heuristic confidence score, not a formal probability
- `total_calibrated`: whether the filled value was reconciled to a known total

When presenting results, link filled values back to `method` and report the proportion of cells by method.

## Quality Checks

Before finishing:

- confirm no observed values changed
- confirm every city-year cell is present in the requested year range
- review all-missing cities and sparse cities
- review yearly total mismatches
- verify non-ASCII city names are unchanged between input and output
- compare logistic and Gompertz results if the model choice is close or materially changes policy conclusions
- run a holdout check when enough observed data exist: remove 15-25% of observed years, refit, and compare RMSE/MAE by city type
- do not describe completed values as measured or official statistics
