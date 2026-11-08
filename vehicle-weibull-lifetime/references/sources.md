# 来源与引用

## 1. 本研究论文

Xu, G., Xiong, X., Liu, W., Yang, S., Miao, Y., Liu, F., Liu, H., Feng, Y., Liang, X., He, Z., & Nakajima, K. (2025). Resource, economic, and carbon benefits of end-of-life trucks' urban mining in China. Communications Earth & Environment, 6, 805. https://doi.org/10.1038/s43247-025-02832-x

论文方法部分给出：

- 使用动态存量驱动的队列平衡模型估算在用车和报废车。
- 报废量由队列投放量、车辆年龄和寿命分布共同决定。
- 两参数 Weibull 分布用于拟合整车寿命分布。
- 形状参数 `beta` 反映函数形态，尺度参数 `eta` 反映数据分散程度。
- 概率密度和累积分布用于计算特定年龄的退出份额。
- 形状参数与尺度参数的数据来源见论文 Supplementary Table 4。

注意：论文 PDF 正文提到 Supplementary Table 4 提供形状和尺度参数，但该 PDF 本身不包含全部补充表数值。当前工作簿中的欧洲乘用车、欧洲商用车和印度参数属于后续整理的参数表或二次计算结果。

## 2. 欧洲乘用车基础参数

Held, M., Rosat, N., Georges, G., Pengg, H., & Boulouchos, K. (2021). Lifespans of passenger cars in Europe: Empirical modelling of fleet turnover dynamics. European Transport Research Review, 13, 9. https://doi.org/10.1186/s12544-020-00464-0

用途：

- 提供欧洲 25 个国家的国家级平均寿命参数和形状参数。
- 工作簿将 Held 等报告的平均寿命参数 `gamma` 转换为标准 Weibull 尺度参数 `eta`，公式为 `eta = gamma / Gamma(1 + 1/beta)`。
- 未覆盖国家的参数按车队更新较慢且车龄结构相近国家的中位数插补。

## 3. 欧洲商用车和客车参数

Xu, G., Xiong, X., Liu, W., et al. (2025). Resource, economic, and carbon benefits of end-of-life trucks' urban mining in China. Communications Earth & Environment, 6, 805. https://doi.org/10.1038/s43247-025-02832-x

Zeng, X., Ali, S. H., Tian, J., & Li, J. (2020). Mapping anthropogenic mineral generation in China and its implications for a circular economy. Nature Communications, 11, 1544. https://doi.org/10.1038/s41467-020-15246-4

Harris, A., Soban, D., Smyth, B. M., & Best, R. (2018). Assessing life cycle impacts and the risk and uncertainty of alternative bus technologies. Renewable and Sustainable Energy Reviews, 97, 569-579. https://doi.org/10.1016/j.rser.2018.08.045

Sokic, M. D., Ilic, I., Manojlovic, V., Markovic, B., Gulisija, Z., Pavlovic, M., & Strbac, N. (2016). Modeling and prediction of the end of life vehicles generation based on a two-parameter Weibull distribution. Acta Polytechnica Hungarica, 13, 41-60.

用途：

- Van、Pickup Truck、Light Trucks 使用 `eta = 8.02`、`beta = 5.58`。
- Medium Trucks 使用 `eta = 10.09`、`beta = 5.58`。
- Heavy Trucks 使用 `eta = 12.80`、`beta = 5.58`。
- Buses 的平均寿命取 12 年，形状参数取 `beta = 4.00`，换算得到 `eta = 13.24` 年。

## 4. 车型修正依据

Zheng, J., Zhou, Y., Yu, R., Zhao, D., Lu, Z., & Zhang, P. (2019). Survival rate of China passenger vehicles: A data-driven approach. Energy Policy, 129, 587-597. https://doi.org/10.1016/j.enpol.2019.02.037

Rith, M., Fillone, A., & Biona, J. B. M. M. (2020). Development and application of the Weibull distribution-based vehicle survivorship models for a metropolis of a developing country. Songklanakarin Journal of Science and Technology, 42(6), 1233-1238. https://doi.org/10.14456/sjst-psu.2020.160

用途：

- 支持 Car、SUV、MPV 生存特征不同的假设。
- 工作簿据此对 SUV 采用 `eta * 1.05`、对 MPV 采用 `eta * 1.08`，并保持 `beta` 不变。
- 上述 5% 和 8% 是国家组基准参数上的研究设置修正，不是文献直接拟合值。

## 5. 新能源整车寿命依据

Nguyen-Tien, V., Zhang, C., Strobl, E., & Elliott, R. J. R. (2025). The closing longevity gap between battery electric vehicles and internal combustion vehicles in Great Britain. Nature Energy, 10, 354-364. https://doi.org/10.1038/s41560-024-01698-1

用途：

- 表明 BEV 预计寿命已接近汽油车。
- 支持基准情景下同一 segment 内不同动力类型采用相同整车寿命参数。
- 新能源短寿命情景只作为敏感性测试，不应被表述为已证实差异。

## 6. 国家分组和车龄数据

- European Automobile Manufacturers' Association (ACEA), 2025. Vehicles on European roads 2025.
- Autoalan Tiedotuskeskus, 2025. Average age of passenger cars in European countries.
- Eurostat, 2026. Transport equipment statistics.
- United Nations Economic Commission for Europe (UNECE), 2026. Road vehicle fleet at 31 December by vehicle category and age group.
- Turkey: TurkStat, 2025.
- Russia: AUTOSTAT, 2024.
- Kazakhstan: Bureau of National Statistics, 2024.
- Ukraine: Institute for Automotive Market Research (Eauto), 2025a, 2025b.

这些来源用于按平均车龄、车龄结构和车队更新速度划分 I、II、III 类国家组，以及为缺少统一统计数据的国家作代理归类。

## 7. 材料流联动

材料含量文件：

`D:\FILES\USST_study\Dong_qiuxia\xXiong\材料含量_Sheet2加入PickupTruck.xlsx`

该文件提供车型和动力类型的整车质量与材料质量占比，可用于把 Weibull 队列报废量转换为材料流。它不提供新的 Weibull 寿命参数。后续计算时：

```text
材料流(t, m) = 报废车辆数(t) * 整车质量 * 材料质量占比(m)
```

动力电池、燃料电池等部件材料流应使用部件寿命模型或退役比例，不得直接等同整车 Weibull 退出。