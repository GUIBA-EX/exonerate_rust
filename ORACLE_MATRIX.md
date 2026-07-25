# 命令行基准覆盖

`是`表示已有永久 upstream 2.4.0 基准；`部分`表示代表性覆盖；`否`表示暂无基准。

| 模型 | 基础输出 | 链方向 | 多记录 | 阈值/Best-N | 次优路径 | 低内存 |
| --- | --- | --- | --- | --- | --- | --- |
| `ungapped` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `ungapped:trans` | 是 | 部分 | 是 | 是 | 是 | 是 |
| `affine:global` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `affine:bestfit` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `affine:local` | 是 | 部分 | 是 | 是 | 是 | 是 |
| `affine:overlap` | 是 | 部分 | 否 | 否 | 是 | 是 |
| `protein2dna` | 是 | 是 | 否 | 部分 | 部分 | 不适用 |
| `protein2dna:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 不适用 |
| `est2genome` | 是 | 部分 | 是 | 是 | 部分 | 是 |
| `protein2genome` | 是 | 是 | 是 | 部分 | 部分 | 是 |
| `protein2genome:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 部分 |
| `coding2coding` | 是 | 不适用 | 否 | 否 | 是 | 不适用 |
| `coding2genome` | 是 | 是 | 是 | 是 | 是 | 是 |
| `cdna2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `genome2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `ner` | 是 | 部分 | 是 | 是 | 是 | 是 |

## 跨模型检查

| 功能 | 状态 |
| --- | --- |
| help、version、别名 | 是 |
| FASTA byte chunk | 是 |
| 反序 ID 的多记录排序 | 是 |
| 默认 score/subopt | 是 |
| 非法参数与类型冲突 | 部分 |
| pretty 与通用 RYO | 是 |
| 页眉/页脚 | 是 |
| GFF2 | 不在范围内 |

低内存基准包含零分 NER 和 splice-boundary `genome2genome` 的原子 `%P`
transition 等价性。
