# 通用 C4 checkpoint traceback

通用 C4 图可包含 NER、intron、phase 等跨行状态。完整 parent 表超过
`--dpmemory` 预算时，执行器使用 checkpoint 回放，而不改变模型语义。

## 不变量

- checkpoint 保存 score、phase shadow、continuation parent frontier 和长状态队列；
- 回放使用完整输入，只限制重算的行区间，保留 splice PSSM 上下文；
- 每个区块按原 transition ID 重建 trace 与 raw trace；
- 预算同时计入 checkpoint payload、滚动行和 replay parent block；
- full 与 checkpoint 必须在得分、tie、坐标、trace 和原子 transition 上一致。

内建 local `ungapped` 与 `ungapped:trans` 保留专用线性执行器；其余通用图在
完整表超出预算时自动选择 checkpoint backend。`--dpmemory 0` 强制该路径。
