# 兼容范围

行为基准为仓库内的 Exonerate 2.4.0（`upstream/src/program/exonerate`）。

| 范围 | 支持 |
| --- | --- |
| 输入 | FASTA、多记录、分块读取 |
| 模型 | 16 个 CLI 模型及短别名 |
| 比对 | DNA IUPAC、BLOSUM62、intron、phase、split codon、frameshift |
| 坐标 | 正反链、正向参考坐标与反向坐标模式 |
| 输出 | sugar、cigar、vulgar、pretty、RYO（含原子 `%P`）、GFF3 |
| 搜索 | 启发式、`--exhaustive`、`--subopt`、`--bestn` |
| 内存 | 专用模型与通用 C4 的 checkpoint traceback |
| CLI | upstream 风格的 verbose、header/footer 与 hostname |

GFF2 不在范围内。

兼容性优先保证得分与 traceback，其次是坐标、排序和报告格式。表中“支持”
表示存在执行路径和回归测试；逐项证据见[命令行基准覆盖](ORACLE_MATRIX.md)。
