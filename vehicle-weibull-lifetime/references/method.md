# 方法说明：车辆两参数 Weibull 寿命模型

## 1. 适用范围

本方法用于整车退出、报废、车队更新、在用车保有量和材料流模型中的车辆寿命假设。模型对象是“整车退出”，不是动力电池、燃料电池堆、发动机、驱动电机等零部件的失效寿命。

原始材料：

- `D:\FILES\USST_study\Dong_qiuxia\xXiong\印度汽车_Weibull寿命参数.xlsx`
- `D:\FILES\USST_study\Dong_qiuxia\xXiong\欧洲威布尔模型参数确定.xlsx`
- `D:\FILES\USST_study\Dong_qiuxia\xXiong\Resource, economic, and carbon benefits of end-of-life trucks’ urban mining in China.pdf`
- `D:\FILES\USST_study\Dong_qiuxia\xXiong\材料含量_Sheet2加入PickupTruck.xlsx` 仅用于后续材料流与质量组成联动，不提供 Weibull 寿命参数。

## 2. 基本公式

设车辆使用年龄为 `k >= 0`（年），尺度参数为 `eta > 0`（年），形状参数为 `beta > 0`。

```text
生存函数：       S(k) = exp[-(k/eta)^beta]
累积报废函数：   F(k) = 1 - S(k) = 1 - exp[-(k/eta)^beta]
概率密度：       f(k) = (beta/eta) * (k/eta)^(beta-1) * exp[-(k/eta)^beta]
风险函数：       h(k) = f(k)/S(k) = (beta/eta) * (k/eta)^(beta-1)
平均寿命：       E[T] = eta * Gamma(1 + 1/beta)
寿命方差：       Var(T) = eta^2 * [Gamma(1 + 2/beta) - Gamma(1 + 1/beta)^2]
中位寿命：       Median(T) = eta * (ln 2)^(1/beta)
q 分位寿命：      T_q = eta * [-ln(1-q)]^(1/beta)
```

注意：这里的 `Gamma(.)` 是 Gamma 函数；文献表格中的平均寿命参数有时也用 `gamma` 表示，两者含义不同。

## 3. 参数换算

### 3.1 已知平均寿命与形状参数

如果文献报告的是平均寿命 `mu = E[T]`，而不是 Weibull 尺度参数，则：

```text
eta = mu / Gamma(1 + 1/beta)
```

例如欧洲客车取平均寿命 12 年、`beta = 4`：

```text
eta = 12 / Gamma(1.25) = 13.239... ≈ 13.24 年
```

### 3.2 已知中位寿命与形状参数

```text
eta = median(T) / (ln 2)^(1/beta)
```

### 3.3 已知某年龄的累计报废概率

如果车辆在年龄 `t` 的累计报废比例为 `p`，则：

```text
eta = t / [-ln(1-p)]^(1/beta)
```

### 3.4 已知尺度参数与形状参数

直接使用：

```text
E[T] = eta * Gamma(1 + 1/beta)
```

不要再次乘除 Gamma 函数，否则会造成重复换算。

## 4. 年尺度队列模型

对每个生产年或投放年队列 `c`，设当年销量或投放量为 `N_c`。若按整数年龄建模，年龄 `k` 的队列退出比例建议使用两个时点生存率之差：

```text
phi(k) = S(k-1) - S(k),  k = 1, 2, ...
```

其中 `S(0) = 1`。则：

```text
队列在年龄 k 的存活量 = N_c * S(k)
队列在年龄 k 的报废量 = N_c * [S(k-1) - S(k)]
```

在日历年 `t` 的报废量可以汇总为：

```text
F_out(t) = sum_c [N_c * phi(t-c)]
```

原始论文的方法部分给出了存量平衡和队列退出关系：

```text
F_in(t) = S(t) - S(t-1) + F_out(t)
F_out(t) = sum_k [F_in(t-k) * phi(k)]
```

论文和 Excel 中的 `phi(k)` 对应“特定年龄的退出/报废份额”。若直接使用连续概率密度 `f(k)` 作为整数年龄的退出权重，权重之和通常不等于 1。实际做年度队列模型时，优先使用 `S(k-1) - S(k)`，或先把离散密度在年龄区间上归一化，再作为 `phi(k)`。

## 5. 参数口径与选择顺序

1. 优先使用最具体的可用参数：国家或地区、车辆大类、segment、动力类型。
2. 如果只有国家组参数，使用国家组 × segment 组合。
3. 如果只有车型参数而没有国家参数，则按车型统一使用，不人为细分国家差异。
4. 如果基准情景没有证据支持动力类型差异，则同一车型或 segment 内各动力类型使用同一套整车寿命分布；动力类型差异放在短寿命或长寿命情景中测试。
5. 每一种参数都应记录为“文献直接参数”“由平均寿命换算的尺度参数”“国家组中位数”或“情景假设”，不得把二次计算结果伪装成文献直接拟合值。

## 6. 情景与敏感性

### 6.1 尺度参数情景

若只调整 `eta`：

```text
eta_scenario = eta_baseline * multiplier
beta_scenario = beta_baseline
```

情景倍数通常直接作用于 `eta`，而不是先作用于平均寿命。除非情景明确指定“平均寿命降低 X%”，此时应先调整 `mu`，再用 `eta = mu / Gamma(1 + 1/beta)` 换算。

### 6.2 形状参数情景

如果只调整 `beta`，平均寿命会同时变化，因为：

```text
E[T] = eta * Gamma(1 + 1/beta)
```

若希望“固定平均寿命、只改变分布离散程度”，必须先给定目标平均寿命 `mu_target`，再计算：

```text
eta_scenario = mu_target / Gamma(1 + 1/beta_scenario)
```

### 6.3 材料中的敏感性设置

- 印度乘用车：`eta` 做 `±10%`，`beta` 做 `±0.3`。
- 欧洲乘用车：基准、轻度短寿命、中度短寿命；PHEV 使用 `0.95η`、`0.90η`，BEV/FCV 使用 `0.90η`、`0.80η`，HEV/ICEV 基准保持 `1.00η`。
- 欧洲商用车：BEV、PHEV、FCV 的 `eta` 可分别降低 10% 和 20%；`beta` 不变。

### 6.4 平均寿命报告

情景结果同时报告：

```text
E[T_scenario] = eta_scenario * Gamma(1 + 1/beta_scenario)
```

否则只报告 `eta` 会掩盖 `beta` 变化对平均寿命和分布形态的影响。

## 7. 校验清单

- `eta > 0`，`beta > 0`。
- `S(0) = 1`。
- `S(k)` 随年龄单调不增。
- `phi(k) >= 0`。
- 在所采用的最大年龄范围内，`sum_k phi(k) <= 1`；延长最大年龄后应逐步趋近 1。
- 平均寿命数值校验：`E[T]` 与 `eta * Gamma(1 + 1/beta)` 的误差小于数值容差。
- 情景换算后重新计算 `E[T]`，不要只检查 `eta`。
- 队列汇总结果与存量平衡关系一致。
- 国家、车型、segment、动力类型和单位的映射与原始工作簿一致。
- 新能源整车寿命情景不能解释为动力电池寿命情景。

## 8. Excel 公式对应关系

印度工作簿中的主要公式：

```text
平均寿命：     =B4*GAMMA(1+1/C4)
eta 低情景：   =B4*0.9
eta 高情景：   =B4*1.1
beta 敏感范围：=C4-0.3 与 C4+0.3
```

欧洲乘用车工作簿中的主要公式：

```text
SUV 尺度参数： =Car尺度参数*1.05
MPV 尺度参数： =Car尺度参数*1.08
形状参数：     与 Car 基准保持一致
```