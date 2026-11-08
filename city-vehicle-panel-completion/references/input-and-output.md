# Input And Output Guide

## Accepted Inputs

The completion script accepts CSV, TSV, XLSX, XLS, and XLSM files.

Long format requires a city column, a year column, and a vehicle-value column. Optional columns are population, GDP per capita, and a repeated yearly total.

Wide format uses one row per city and year columns such as `2008`, `2009`, or `2010年`. Use `--format wide` if automatic detection is ambiguous.

## Text Encoding

`--encoding auto` is the default. Detection order:

1. UTF-8, UTF-16, or UTF-32 byte-order mark
2. strict UTF-8
3. `charset-normalizer`
4. a conservative Western-European correction for CP1250/CP1252 ambiguity
5. strict fallback candidates

For CSV files containing Norwegian, Swedish, Danish, German, French, or other Western-European names, `cp1252`, `iso-8859-1`, and `utf-8-sig` are supported. Chinese legacy files can be read with `gb18030`; Japanese files with `shift_jis`.

Use an explicit codec whenever auto-detection is uncertain:

```text
--encoding cp1252
--encoding iso-8859-1
--encoding gb18030
--encoding shift_jis
```

CSV output uses a Unicode-safe derivative of the source encoding by default. Change it with `--output-encoding`. The generated `completion_workbook.xlsx` stores Unicode directly and avoids CSV encoding problems.

Fonts do not determine how Unicode is stored in XLSX. Fonts matter when viewing or rendering the data: a font without Chinese or Norwegian glyphs can display boxes even when the underlying code points are correct. The figure generator selects an installed CJK-capable font when available and falls back to DejaVu Sans for Latin characters.

## Units

Use `--value-unit counts` when the vehicle column is a stock/count. If population is available, the default basis converts to vehicles per 1,000 people for curve fitting and converts predictions back to counts.

Use `--value-unit per_1000` when the source is already a normalized ownership rate. Do not divide by population again.

Vehicle counts, private passenger cars, and civil motor vehicles are different variables. Do not combine them without a documented conversion rule.

## Missing Values

The following tokens are treated as missing:

```text
blank, NA, N/A, NaN, None, Null, ., .., ..., -, --, —, –, n.a.
```

Zeros remain zeros. Do not infer that a missing observation means no vehicles.

## Duplicate Rows

Use `--aggregate-duplicates sum` when duplicate rows are components that should be added. Use `mean` when duplicates represent repeated measurements. Population, GDP, and yearly totals are averaged/medianed rather than summed.

## Yearly Totals

The total can be supplied in the same file with `--total-col`, or in a separate file with:

```text
--total-file
--total-year-col
--total-value-col
```

The total must use the same geographical coverage and vehicle definition as the city panel. A national total cannot be used to reconcile a partial set of cities without documenting the coverage difference.

## Method Field

- `observed`: original observed value
- `city_logistic`: city-specific logistic curve
- `city_gompertz`: city-specific Gompertz curve
- `pooled_logistic`: common logistic shape with a city-specific scale
- `pooled_gompertz`: common Gompertz shape with a city-specific scale
- `pooled_all_missing_logistic`: no observed city values; shape and scale borrowed from other cities
- `pooled_all_missing_gompertz`: no observed city values; shape and scale borrowed from other cities
- suffix `+total`: value was subsequently allocated under a yearly total constraint

## Confidence Field

`confidence` is a transparent heuristic used for triage:

- city-specific fit: `0.90`
- sparse city with observed values: approximately `0.75`
- all-missing city with covariate neighbours: approximately `0.45`
- all-missing city with only a pooled median: approximately `0.35`
- no usable cross-city information: approximately `0.10`

It is not a statistically calibrated probability. Use interval width, holdout error, and the report diagnostics for formal uncertainty.
