# 通用 C4 checkpoint traceback 设计

本设计以 Exonerate 2.4.0 的 `upstream/src/c4/viterbi.c` 和
`upstream/src/c4/optimal.c` 为依据。它适用于含 NER、intron 与 phase
长状态的通用 C4 图；现有模型专用 checkpoint 不替代这里的通用机制。

## 上游不变量

1. 是否使用 reduced space 由完整 traceback 表加滚动行的估计内存决定；
   尺寸溢出也必须进入 reduced space。
2. checkpoint 保留的是完整 C4 cell：得分、所有 shadow，以及额外的
   checkpoint shadow。只保存得分不足以恢复长状态路径。
3. checkpoint shadow 编码 `(state, rolling-row, query-position)`，用于在
   相邻 checkpoint 之间恢复 continuation 的开始状态和 cell。
4. continuation 使用 corner scope，在两个 checkpoint 间重算子区域；
   递归持续到子区域可容纳完整 traceback 表为止。
5. 最终 path 由子区域 traceback 顺序拼接，transition 必须映射回原始
   C4 图的 transition ID。

## Rust 侧实现

1. `GenericScoreMatrix`、`GenericShadowMatrix` 和 `GenericParentMatrix`
   使用带逻辑行标签的全表/滚动行统一接口，物理槽复用不会暴露为错误逻辑行。
2. checkpoint 同时保存过去依赖行、已由普通边写入的未来 frontier、
   frontier parent continuation、phase donor shadow 和
   `GenericLongStatePayload`。
3. 前向扫描只保留依赖窗口并周期性创建 checkpoint；traceback 从终点开始，
   逐块恢复边界、重算 parent，再按原 transition ID 拼接 trace/raw trace。
4. replay 始终使用完整输入序列，只限制计算行范围，避免 splice PSSM 在子序列
   边界上改变分数。
5. full/checkpoint 等价回归覆盖 NER、query/target/joint intron、phase 0/1/2、
   反向目标和 suboptimal forbidden-pair 枚举。

`align_model_ir_with_dp_memory` 会计算完整 generic 表的计划量；若全表超过请求
预算则自动选择 checkpoint backend。内建 local `ungapped` 和
`ungapped:trans` 继续使用更小的专用线性执行器。CLI 的 NER 与
`ungapped:trans` 普通及 suboptimal 路径均可在 `--dpmemory 0` 下保持精确
输出。

## 验收

- `--dpmemory 0` 不得分配全矩阵 parent 表。
- 结果必须与 full generic executor 完全一致，包括同分 tie、原子转移和
  split-codon/phase/intron 顺序。
- 预算计划必须计入 checkpoint payload、递归子区域 parent 表和滚动行；
  它是计划预算，不宣称为 RSS 硬上限。
