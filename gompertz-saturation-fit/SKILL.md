---
name: gompertz-saturation-fit
description: Fit and interpret Gompertz saturation (S-curve) models for stock/ownership data, with free or fixed saturation levels, 95% CIs, scenario comparison, chronological validation, extrapolation and spatial allocation weights. Use when the user mentions Gompertz 模型 / S 曲线 / 饱和点 / 餍足点 / 保有量预测 / 人均保有率与人均GDP 的关系 / 报废量或渗透率的长期预测, or asks to fit a saturation curve, choose a saturation ceiling, or reproduce Zhang et al. 2020's China end-of-life vehicle PGM recycling method. Works on CSV/XLSX year-by-region or year-by-country panels.
metadata:
  short-description: Fit Gompertz saturation curves with scenario comparison and validation
---

# Gompertz 饱和模型拟合与预测

用可复现的统计流程拟合 Gompertz S 型饱和曲线（人均保有量、渗透率、保有规模等），
并把"饱和点到底是数据给的还是假设给的"这件事说清楚。

模型形式（论文/主流写法）：

```
y = γ · exp( α · exp( β · x ) )      α < 0, β < 0，γ 为餍足点
```

## 何时用

- 要拟合"随收入/时间趋于饱和"的存量或比率（汽车保有量、空调普及率、EV 渗透率…）
- 要做长期外推并给出饱和上限，或要比较不同饱和点假设下的结论
- 要把全国总量按各省/各市分配（Gompertz 权重法）
- 要复现张安迎等(2020)《中国报废汽车中铂族金属回收潜力估算》的保有量建模部分

## 工作流

1. **确认数据结构**：x 是"人均 GDP"还是"年份"？y 是人均率还是总量？
   有存量+人口时用 `--stock-col/--pop-col` 让脚本自动算人均率。
   确认 x 的单位（万元/人 vs 元/人）——单位改变 β 的数量级。
2. **先外推自变量**（若 x 是人均 GDP/人口）：用 `scripts/forecast_trend.py` 做
   Brown/Holt 线性趋势外推，脚本会按一步预测 RMSE 自动选模型。
3. **拟合**：`scripts/fit_gompertz.py` 同时跑自由 γ 与若干固定 γ 情景。
4. **判读可识别性**：看拐点 x\* 是否落在观测 x 范围内、γ 的 95% CI 是否过宽、
   留出验证误差是否随时间放大。
5. **出结论**：参数表 + 情景对比表 + 逐期预测误差 + 预测区间，标明外推前提。

## 常用命令

```bash
# 拟合 + 多情景 + 验证 + 预测带
python scripts/fit_gompertz.py \
  --input data.csv --x-col gdp_pc --y-col cars_pc \
  --gamma 0.62 --gamma 0.45 \
  --bootstrap 500 --backtest --compare-logistic \
  --predict-range 7:25:0.5 \
  --out out_dir

# 由存量+人口自动构造人均率
python scripts/fit_gompertz.py --input data.csv \
  --x-col gdp_pc --stock-col vehicles --pop-col population --out out_dir

# 空间分配权重（论文第 1.2 节方法）
python scripts/fit_gompertz.py --input province_panel.csv \
  --x-col gdp_pc --y-col cars_pc --weight-col province \
  --weight-year 2018 --allocate 2535.05 --out out_dir

# 外推人均GDP与人口（论文第 1.1 节 Brown 线性趋势）
python scripts/forecast_trend.py --input drivers.csv \
  --year-col year --value-col gdp_pc --value-col population \
  --horizon 32 --out trend_dir
```

改过脚本后跑一遍自检：`python scripts/selftest.py`

## 选饱和点 γ 的决策规则

1. 自由 γ 的 95% CI 不宽、且拐点落在数据范围内 → 用自由 γ
2. CI 极宽或拐点在数据范围外 → 饱和值是外推假设，改用外部基准固定 γ
   （OECD 通用 0.62 辆/人；低收入国家 0.3~0.45；文献常用 0.4~0.7），并做敏感性分析
3. 不同 γ 情景结论相反 → 结论必须写成"在 γ=X 情景下……"，不得给单一数字
4. 任何情况下都要同时报自由解与替代情景（论文表1/表2 的做法）

## 输出契约

`--out` 目录内：

| 文件 | 内容 |
|---|---|
| `fit_report.md` | 方法、情景对比表、参数与 95% CI、诊断、警告、权重表 |
| `parameters.csv` | 每个情景的 γ/α/β、R²、adj R²、RMSE、MAPE、AIC、拐点 |
| `predictions.csv` | 预测值 + 95% 预测带（`y_lower95`/`y_upper95`） |
| `spatial_weights.csv` | 归一化分配权重（及 `--allocate` 的分配量） |
| `backtest.csv` | 按时间顺序的分折留出验证 |
| `fit_summary.json` | 机器可读摘要 |
| `gompertz_fit.png` | 观测点、拟合曲线、预测带、饱和线、拐点 |

## 硬性规则

- y 必须为正；点数 < 5 直接报错，不要强行拟合
- x 是年份量级时脚本内部会自动平移（x' = x − min(x)）；报告里 α 若显示
  "(shifted)"，就在文中写明平移量，**不要**引用溢出前的巨大数值
- 外推值与观测值必须区分标注；留出验证按时间顺序做，不用随机抽折
- γ 的 CI、拐点位置、预测区间要一并给出，不能只给 R²
- 结构性变化（电动化、政策退坡、共享出行）会打破单饱和假设：留出误差后期放大时，
  明确写出外推前提已失效
- 权重法只分配份额，不替代各地区的独立模型

## 参考资料

- [references/model.md](references/model.md) — 数学推导、剖面估计法、数值陷阱、
  可识别性判读、不确定度、模型选择、x 的选择与单位、空间权重法、常见坑清单
- [references/case-zhang2020.md](references/case-zhang2020.md) — 论文参数表、复现命令、
  下游链条（报废量→Weibull→物质流）、论文自述局限与写作模板
- [references/example-norway.md](references/example-norway.md) — 真实数据示例（挪威私人汽车
  人均保有量），含结果解读与"结构变化导致留出误差放大"的读法
- `assets/example_norway_cars.csv` — 可直接运行的示例数据
