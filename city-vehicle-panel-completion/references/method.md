# Method Notes

## Logistic Growth

For city `i` and year `t`:

```text
Y_it = K_i / (1 + exp(-r_i * (t - t0_i))) + e_it
```

- `K_i`: saturation level
- `r_i`: growth rate
- `t0_i`: inflection year

The model is estimated by nonlinear least squares with bounded parameters and multiple starting values.

## Gompertz Growth

```text
Y_it = K_i * exp(-exp(-b_i * (t - c_i))) + e_it
```

Gompertz allows a slower approach to saturation and is widely used in vehicle ownership forecasting. The script compares logistic and Gompertz by AICc when `--model auto` is used.

## Rate Basis

When population is available and the vehicle input is a count, the script fits:

```text
y_it = 1000 * vehicles_it / population_it
```

Predictions are converted back to counts after fitting. This avoids letting city size dominate a pooled curve.

## Sparse Cities

For cities with too few observations:

1. Each observed city trajectory is normalized by its median positive level.
2. A pooled logistic or Gompertz shape is fitted to all normalized observations.
3. A city-specific scale is estimated from observed values relative to the pooled shape.
4. Cities with no observations use the scale of population/GDP-near observed cities when covariates exist.

This is a practical fallback, not a substitute for a fully specified nonlinear mixed-effects model.

## Total Reconciliation

For year `t`, let `O_t` be observed cities, `M_t` missing cities, and `T_t` the reported total.

```text
R_t = T_t - sum_{i in O_t} observed_it
```

The missing cells are allocated in proportion to preliminary model values:

```text
completed_it = R_t * preliminary_it / sum_{j in M_t} preliminary_jt
```

Population is used as a fallback weight if a model value is missing or non-positive.

Observed values are never changed. If `R_t < 0`, calibration is skipped for that year and the condition is reported.

## Missingness Assumptions

Using partial city trajectories without imputation is valid under a missing-at-random working assumption. The model borrows strength across time and across cities.

If reporting is systematically related to the unobserved ownership level, results can be biased. For example, smaller or lower-capacity cities may be less likely to report. Use these checks:

- compare missingness by population, GDP, region, and year
- fit with complete cases only
- fit with the pooled fallback disabled or restricted
- conduct leave-one-city-out and leave-one-year-out validation
- report ranges rather than a single completed series when ambiguity is material

## Formal Alternatives

For an academic paper or high-stakes forecast, consider:

- nonlinear mixed-effects logistic/Gompertz models with random effects on `K`, `r`, and `t0`
- Bayesian hierarchical models with informative priors or literature-based saturation bounds
- state-space models with latent trajectories and missing observations
- multiple imputation, with all model uncertainty propagated before combining estimates
- constrained Bayesian or generalized least-squares reconciliation to yearly totals

## Validation

Preferred validation for a 2008-2019 city panel:

1. Hold out 2017-2019 for cities with late observations.
2. Randomly mask 15-25% of observed city-years.
3. Hold out entire cities, especially sparse and all-missing cities.
4. Compare city-specific, pooled, and total-calibrated predictions.
5. Report MAE/RMSE by population group and observation count.

## Selected Literature

- de Jong, G., Fox, J., Daly, A., Pieters, M., & Smit, R. (2004). Comparison of car ownership models. *Transport Reviews*, 24(4), 379-408. https://doi.org/10.1080/0144164032000138733
- Dargay, J., & Gately, D. (1999). Income's effect on car and vehicle ownership, worldwide: 1960-2015. *Transportation Research Part A*, 33(2), 101-138. https://doi.org/10.1016/S0965-8564(98)00026-3
- Clark, S. D. (2006). Estimating local car ownership models. *Journal of Transport Geography*, 15(3), 184-197. https://doi.org/10.1016/j.jtrangeo.2006.02.014
- Wu, T., Zhao, H., & Ou, X. (2014). Vehicle ownership analysis based on GDP per capita in China: 1963-2050. *Sustainability*, 6(8), 4877-4899. https://doi.org/10.3390/su6084877
- Li, Y., Lu, M., Chen, Y., & Hu, Y. (2019). Exploration of sustainable urban transportation development in China through the forecast of private vehicle ownership. *Sustainability*, 11(16), 4259. https://doi.org/10.3390/su11164259
- Lu, H., Ma, H., Sun, Z., & Wang, J. (2017). Analysis and prediction on vehicle ownership based on an improved stochastic Gompertz diffusion process. *Journal of Advanced Transportation*, 2017, 4013875. https://doi.org/10.1155/2017/4013875
