---
name: vehicle-weibull-lifetime
description: Derive, validate, and apply two-parameter Weibull vehicle lifetime parameters for end-of-life vehicle, fleet turnover, and material-flow modeling. Use when working with vehicle survival or scrappage rates, scale and shape parameters, mean-lifetime conversions, country, segment, or powertrain lifetime scenarios, or parameter tables derived from the supplied truck and passenger-car Weibull workbooks.
metadata:
  short-description: Vehicle Weibull lifetime parameters
---

# Vehicle Weibull Lifetime

Use two-parameter Weibull survival to translate vehicle lifetime assumptions into cohort survival and scrappage rates for dynamic fleet, end-of-life vehicle, retired battery, and material-flow models.

## Core Model

Let vehicle age be `k >= 0` years, scale `eta > 0` years, and shape `beta > 0`:

```text
Survival:        S(k) = exp[-(k/eta)^beta]
CDF:             F(k) = 1 - exp[-(k/eta)^beta]
Density:         f(k) = (beta/eta) * (k/eta)^(beta-1) * exp[-(k/eta)^beta]
Mean lifetime:   E[T] = eta * Gamma(1 + 1/beta)
Scale from mean: eta  = E[T] / Gamma(1 + 1/beta)
```

When a source reports a mean lifetime and shape parameter, recover `eta` with the scale-from-mean formula. Do not treat a reported mean lifetime as `eta`. When a source already reports `eta`, use it directly.

For annual cohort modeling, use the discrete-age scrappage probability implied by the distribution, for example the probability that a vehicle survives to age `k-1` but not to age `k`:

```text
phi(k) = S(k-1) - S(k)
```

Keep ages and `eta` in the same unit. The supplied vehicle tables use years.

## Workflow

1. Identify whether the source provides `eta`, mean lifetime, or both. Convert only when needed, and record whether the parameter is observed, transformed from a reported mean, or an analyst-assumed scenario.
2. Resolve the parameter hierarchy before calculating: country or region, vehicle class, segment, and powertrain. Prefer the most specific sourced parameter available; otherwise use the documented fallback rule.
3. Apply scenario multipliers multiplicatively to `eta` unless the source explicitly adjusts the mean lifetime, and keep `beta` unchanged unless the scenario says otherwise.
4. Generate annual survival and scrappage rates with `scripts/weibull_lifetime.py`, then reconcile them to the fleet balance:
   `F_out(t) = sum_k [F_in(t-k) * phi(k)]`.
5. Report the parameter value, unit, scenario, source, and whether each value is directly sourced or derived. Preserve the distinction between vehicle retirement and component failure, especially for batteries and fuel-cell stacks.
6. Validate that `S(0) = 1`, survival is non-increasing, `phi(k) >= 0`, `sum_k phi(k) <= 1`, and mean lifetime equals `eta * Gamma(1 + 1/beta)` within numerical tolerance.

Read `references/method.md` when deriving parameters, interpreting the source workbooks, or documenting assumptions. Read `references/parameters.md` when selecting or reporting supplied India, Europe, or commercial-vehicle parameters. Read `references/sources.md` when citing the underlying studies.

## Use The Calculation Script

Run the bundled script from this skill directory:

```bash
python scripts/weibull_lifetime.py mean --eta 15.7 --beta 2.6
python scripts/weibull_lifetime.py eta --mean 12 --beta 4
python scripts/weibull_lifetime.py survival --eta 15.7 --beta 2.6 --max-age 30
python scripts/weibull_lifetime.py scrap --eta 15.7 --beta 2.6 --max-age 30
python scripts/weibull_lifetime.py table --input cohort_sales.csv --year-col year --sales-col sales --eta 15.7 --beta 2.6 --start-year 2025 --end-year 2050
```

For repeated analytical work, import the functions from the script instead of reimplementing the formulas.

## Interpretation Rules

- `eta` controls the lifetime scale: a larger `eta` shifts retirement later.
- `beta` controls concentration and the hazard pattern: larger values generally concentrate scrappage around the characteristic age and imply increasing hazard; smaller values produce a wider retirement-age spread.
- A Weibull model fitted to vehicle retirement describes the whole vehicle exit process, not battery, fuel-cell-stack, engine, or motor reliability.
- Do not assume BEV, PHEV, FCV, HEV, and ICEV have systematically different whole-vehicle lifetimes without data. The supplied workbooks use the same baseline distribution within a vehicle class or segment and treat shorter new-energy-vehicle life as a scenario.
- Do not assign country-specific commercial-vehicle parameters when the source only supports a vehicle-type parameter. State the fallback grouping explicitly.

## Output Contract

When producing parameters or model inputs, include at minimum: grouping key, `eta`, `beta`, unit, baseline or scenario label, and source or derivation note. When producing cohort outputs, retain year, vehicle age, survival, and scrappage probability columns so the result can be audited and joined back to fleet data.

## Encoding Policy

Text and CSV inputs are not assumed to be UTF-8. Preserve the source encoding by default.

- Read with strict decoding. Never use replacement decoding or silently fall back to another codec.
- Detect encoding only when unambiguous: a BOM or pure ASCII. For other bytes, require `--source-encoding`; do not guess between UTF-8 and legacy encodings.
- Accept explicit encodings such as CP1252, Windows-1252, ISO-8859-1, GBK, Shift_JIS, and UTF-8.
- When an output is derived from an input, default the output to the detected or supplied source encoding.
- For a new output without a source file, require `--output-encoding`; do not assume UTF-8.
- Convert only when the caller explicitly supplies a different `--output-encoding`. Back up an existing destination first and report `source encoding -> output encoding`.
- Preserve non-ASCII characters exactly, including but not limited to ø, å, ö, ñ, ç, ä, and ü. Never replace undecodable characters with `�`.
- Do not treat Excel `.xlsx` internals as delimited text. Excel XML is normally UTF-8; the policy applies to external CSV and text data exchanged with the workbook, not to XML parsing inside the workbook.

Examples:

```bash
python scripts/weibull_lifetime.py detect-encoding --input cohort_sales.csv
python scripts/weibull_lifetime.py table --input cohort_sales.csv --source-encoding CP1252 --eta 15.7 --beta 2.6 --start-year 2025 --end-year 2050
python scripts/weibull_lifetime.py table --input cohort_sales.csv --source-encoding CP1252 --output-encoding UTF-8 --out cohort_outflows_utf8.csv
python scripts/weibull_lifetime.py survival --eta 15.7 --beta 2.6 --max-age 30 --out survival.csv --output-encoding Windows-1252
```

If encoding detection fails, stop and ask for the original encoding. Do not inspect a few characters and infer a codec from visual plausibility.
