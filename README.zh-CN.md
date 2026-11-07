# uwcsea-igem-2027

UWCSEA iGEM 2027 赛季仓库。软件包名为 **TipCraft**；项目名称仍由 Kevin 与全队讨论决定。Plug-n-Phage、TailorMade、SOCKET 仅为建议。

## 当前交付物

- [英文角色申报稿](docs/role-answer.md)：软件管线、建模、传感器算法与硬件、wiki 与数据基础设施；等待本人核对后提交。
- [一页计划](docs/00-project-plan.md)：目标、预算假设、分级变化和拟议决策门；英文首页同时展示全文。
- [四项文献核验](docs/02-literature.md)：有出处、有限定条件，不把群聊陈述当成事实。
- [可运行回溯评估](software/tipcraft/benchmarks/retro-t7-t3/README.md)：文献观测、基线和独立预测评估接口。

371–553 是已发表结构片段，不能据此认定整个片段可以随意替换。T7/T3 交换有发表先例；不等于通用插座已成立。天然 T7 的 LPS 识别与未来蛋白受体定向必须区分。BindCraft 的结合成功率不能当作尾丝组装成功率，后者目前未知。来源见核验表。

## 运行

在仓库根目录使用 Python 3.11+ 与 uv：

```bash
uv run --project software/tipcraft tipcraft benchmark --cases software/tipcraft/benchmarks/retro-t7-t3/cases.csv --output data/processed/retro-t7-t3.json
uv run --project software/tipcraft python -m unittest discover -s software/tipcraft/tests -v
```

当前版本是证据与评估工具。已有基线不代表设计管线验证通过；结构预测、设计和 DNA 订单接口均未实现。没有伪造模型得分、实验记录或订单。

## 决策与协作

12 月和 2 月决策门的具体日期、负责人、预算均待全队确认。没有确认的实验室、导师、材料资格和经费时，继续推进计算交付物。传感器备选线独立放在 `sensor/`。

代码使用 MIT，原创文档使用 CC BY 4.0；第三方资料遵守原许可。私人群聊与附件不进入仓库。角色稿和仓库链接尚未发送到群聊或 Classroom。
