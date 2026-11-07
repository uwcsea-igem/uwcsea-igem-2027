# uwcsea-igem-2027

UWCSEA iGEM 2027 赛季仓库。软件包名为 **TipCraft**；项目名称仍由 Kevin 与全队讨论决定。Plug-n-Phage、TailorMade、SOCKET 仅为建议。

## 当前交付物

- [长期平台路线](docs/10-platform-roadmap.md)：以可复用识别接口标准、公开数据与 TipCraft 为产物，分阶段验证宿主预测、结合和组装能力。
- [官方规则与评审证据](docs/11-igem-rule-evidence.md)、[前沿先例与创新假设](docs/12-prior-art-and-innovation.md)：按现行规则准备 2027 项目，以可验证贡献支撑奖项目标。

- [英文角色申报稿](docs/role-answer.md)：软件管线、建模、传感器算法与硬件、wiki 与数据基础设施；等待本人核对后提交。
- [一页计划](docs/00-project-plan.md)：目标、预算假设、分级变化和拟议决策门；英文首页同时展示全文。
- [四项文献核验](docs/02-literature.md)：有出处、有限定条件，不把群聊陈述当成事实。
- [可运行回溯评估](software/tipcraft/benchmarks/retro-t7-t3/README.md)：文献观测、基线和独立预测评估接口。

371–553 是已发表结构片段，不能据此认定整个片段可以随意替换。T7/T3 交换有发表先例；不等于通用插座已成立。天然 T7 的 LPS 识别与未来蛋白受体定向必须区分。BindCraft 的结合成功率不能当作尾丝组装成功率，后者目前未知。来源见核验表。

[最新受体讨论核查](docs/21-receptor-dependence-and-source-check.md)已对照 Classroom 和原始论文：团队澄清了两类问题的意图；waaC/F/G/L 不能统一标成 T7 受体阴性对照，敲除后的噬菌斑变化也不能单独证明直接结合或组装成功。

## 运行

在仓库根目录使用 Python 3.11+ 与 uv：

```bash
uv run --project software/tipcraft tipcraft benchmark --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --output data/processed/retro-t7-t3.json
uv run --project software/tipcraft python -m unittest discover -s software/tipcraft/tests -v
```

当前版本是证据与评估工具。已实际运行独立 WIsH 基因组评分：覆盖六条天然噬菌体观测、四条交换案例明确弃权；T3 的已知阳性排在两个阴性之后。详见 [12 月决策门报告](docs/07-december-gate-report.md)。匿名重跑结果一致，但方法选择时已见过标签，不能宣称盲测验证。结构预测、设计和 DNA 订单接口仍未实现。没有伪造模型得分、实验记录或订单。

下一方法 GenoPHI 的 94 个噬菌体 ID 和全部 402 个训练宿主 ID，现已通过原始实验记录对应到公开存档文件元数据。[训练来源准入与分区](docs/20-source-admission-and-partitions.md)记录了命名映射及当前缺少完整宿主文件／ANI 排除的 HOLD 状态；实验菌株身份和预测能力仍需各自的证据。

## 决策与协作

12 月和 2 月决策门的具体日期、负责人、预算均待全队确认。没有确认的实验室、导师、材料资格和经费时，继续推进计算交付物。传感器备选线独立放在 `sensor/`。

代码使用 MIT，原创文档使用 CC BY 4.0；第三方资料遵守原许可。私人群聊与附件不进入仓库。角色申报以实际题目中的 9 号职责为主，协作承担 2 号分析工作；1 号已由 Kyungmin 认领。Classroom 提交状态见角色申报记录。
