# 命令行基准覆盖

`是`表示已有永久回归与 upstream 2.4.0 对照；`部分`表示只覆盖了代表性
场景；`否`表示尚无永久基准。

| 模型 | 基础输出 | 链方向 | 多记录 | 阈值/Best-N | 次优路径 | 低内存 |
| --- | --- | --- | --- | --- | --- | --- |
| `ungapped` | 是 | 部分 | 否 | 部分 | 是 | 不适用 |
| `ungapped:trans` | 是 | 部分 | 是 | 是 | 部分 | 不适用 |
| `affine:global` | 是 | 部分 | 否 | 否 | 否 | 部分 |
| `affine:bestfit` | 是 | 部分 | 否 | 否 | 否 | 部分 |
| `affine:local` | 是 | 部分 | 是 | 是 | 是 | 部分 |
| `affine:overlap` | 是 | 部分 | 否 | 否 | 否 | 部分 |
| `protein2dna` | 是 | 是 | 否 | 部分 | 部分 | 不适用 |
| `protein2dna:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 不适用 |
| `est2genome` | 是 | 部分 | 是 | 是 | 部分 | 是 |
| `protein2genome` | 是 | 是 | 是 | 部分 | 部分 | 是 |
| `protein2genome:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 部分 |
| `coding2coding` | 是 | 不适用 | 否 | 否 | 是 | 不适用 |
| `coding2genome` | 是 | 是 | 是 | 是 | 是 | 是 |
| `cdna2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `genome2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `ner` | 是 | 部分 | 是 | 是 | 部分 | 不适用 |

## 跨模型检查

| 功能 | 状态 |
| --- | --- |
| help、version、别名 | 是 |
| FASTA byte chunk | 是 |
| 反序 ID 的多记录排序 | 是 |
| 默认 score/subopt | 部分 |
| 非法参数与类型冲突 | 部分 |
| pretty 与通用 RYO | 部分 |
| 页眉/页脚 | 否 |
| GFF2 | 不在范围内 |

下一批优先补 affine 非 local scope、各模型 tie order、pretty/RYO 和
header/footer。所有 fixture 必须小、确定且可在 CI 中稳定运行。
