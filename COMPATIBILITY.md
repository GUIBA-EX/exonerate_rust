# 兼容范围

行为基准为仓库中的 Exonerate 2.4.0：
`upstream/src/program/exonerate`。

| 范围 | 状态 |
| --- | --- |
| FASTA、多记录和分块读取 | 已实现 |
| DNA IUPAC 与蛋白质 BLOSUM62 打分 | 已实现 |
| 16 个 CLI 模型及短别名 | 已实现 |
| 正反链、正向参考坐标与反向坐标模式 | 已实现 |
| intron、phase、split codon、frameshift | 已实现 |
| sugar、cigar、vulgar、pretty alignment | 已实现 |
| 常用 RYO 与原子 `%P` 转移 | 已实现，仍在扩充基准 |
| `--subopt` 与 `--bestn` | 主要模型已实现；部分 affine scope 待补 |
| 启发式搜索与 `--exhaustive` | 已实现 |
| `--dpmemory` checkpoint traceback | 主要模型已实现 |
| GFF3 | 项目格式，已实现 |
| GFF2 | 不在范围内 |

“已实现”表示存在可执行路径和回归测试，不等于所有参数组合均已证明
逐字节兼容。精确证据见 [命令行基准覆盖](ORACLE_MATRIX.md)。

兼容性优先级依次为：

1. 得分与 traceback 路径；
2. 坐标、链方向和排序；
3. sugar/cigar/vulgar/RYO；
4. 错误、header/footer 等 CLI 外围行为。
