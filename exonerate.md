# 架构

## 目标

exonerate-rs 以 upstream 2.4.0 为行为基准，不做逐行 C 翻译。核心设计是：

> 模型描述语义，DP 执行模型，traceback 保存事实，报告层负责格式。

兼容模式优先保证得分、路径、坐标和排序；启发式模式只负责缩小候选
区域，最终路径仍由精确 DP 产生。

## 工作区

```text
crates/
├── exonerate-core/  # 模型、DP、traceback、启发式和报告
└── exonerate-cli/   # 参数、FASTA、排序和输出
```

## 核心数据流

```text
FASTA
  → Sequence
  → Model / ModelIr
  → exhaustive 或 heuristic region
  → full-memory 或 checkpoint DP
  → Alignment + RawStep + TraceRun
  → sugar / cigar / vulgar / pretty / RYO / GFF3
```

### 模型

`ModelIr` 由状态、转移、坐标推进量、打分 kernel、scope 和 tie policy
组成。模型关闭时必须验证：

- start/end 可达；
- epsilon 无环；
- 循环会推进坐标；
- span 与 intron 边界合法；
- 状态和转移顺序稳定。

转移顺序属于兼容行为，因为同分路径依赖 tie-breaking。

### 动态规划

执行器只接受模型和序列，不负责格式化。主要路径：

- full-memory：保存完整 parent；
- checkpoint：保存分段 score，按需重算 parent；
- heuristic：先找候选区域，再调用相同精确 DP。

`--dpmemory` 选择 checkpoint 计划。它限制 DP 规划内存，不承诺限制整个
进程 RSS。

### 回溯

`RawStep` 保存原子转移及得分，供 `%P` 使用；`TraceRun` 保存合并后的
生物事件，供普通报告使用。两者不能互相替代。

支持的事件包括 match、gap、NER、splice 5′/3′、intron、split codon
和 frameshift。

### 坐标

内部统一使用零基、半开区间。反向链保留遍历方向；报告层根据
`--forwardcoordinates` 决定是否投影到正向参考坐标。

### 次优路径

Waterman–Eggert 枚举通过禁止已使用的 equivalenced pairs 生成下一条
路径。`bestn` 按 query 排序，并保留 cutoff 处分数相同的全部路径。

## 不变量

- score、trace 和坐标必须守恒；
- full-memory 与 checkpoint 输出相同；
- heuristic 命中区域时与 exhaustive 输出相同；
- 反向链不能改变生物事件顺序；
- 输出层不能重新推断 DP 状态；
- GFF3 是项目格式，GFF2 不在范围内。

## 验证

测试分三层：

1. core：短序列穷举、状态图、打分和坐标不变量；
2. 命令行：固定 fixture 与 upstream 的字节级基准；
3. 资源：full/checkpoint 等价、时间和 peak RSS。

CI 执行 rustfmt、严格 clippy、测试、doctest 和 rustdoc；大型性能 fixture
不进入普通 CI。
