# 开发状态

更新时间：2026-07-26。当前发布范围内没有已知阻塞项。

## 已完成

- upstream CLI 列出的 16 个模型均有 Rust 执行路径。
- 通用 C4 状态图支持有限状态、epsilon、NER、query/target/joint intron
  和 phase 0/1/2 split codon。
- 复杂模型支持 codon gap、双侧 frameshift、splice PSSM 与原子 traceback。
- 主要模型支持 Waterman–Eggert pair exclusion、`--subopt` 和 `--bestn`。
- affine、EST、protein/coding/cDNA/genome 模型已接入 checkpoint traceback。
- 通用 C4 API 在完整表超过 `--dpmemory` 计划时自动使用 checkpoint backend；
  NER、query/target/joint intron 和 phase 0/1/2 均保留精确 traceback。
- 通用 C4 checkpoint 联合保存 score/shadow 历史、未来 continuation parent
  frontier 和 NER/intron 持久候选队列；预算规划会计入全部常驻快照、
  rolling rows 和 replay parent block。
- NER 与 `ungapped:trans` 的普通及 suboptimal CLI 路径可在
  `--dpmemory 0` 下保持与全表相同的输出。
- 通用 `ungapped` 与 `ungapped:trans` IR 的 local 常规路径已有精确线性
  空间分派，作为 C4 reduced-space 路线的简单图基线。
- generic full/linear/checkpoint 执行共用带逻辑行标签的 score、phase shadow
  和 parent matrix 接口。
- genome2genome checkpoint 回放会恢复 score history 与 intron queue，并在完整
  query 上限制计算行，保留 checkpoint 边界后的 splice PSSM 上下文。
- DNA 与翻译模型具有启发式候选区域，并可用 `--exhaustive` 回退。
- `scripts/validate_peak_rss.sh` 会验证 full/checkpoint 输出逐字节一致并比较
  release peak RSS；500×500 NER 实测约从 331608 KiB 降至 68700 KiB。
- CLI 的默认、显式和负数 verbose 行为及真实 hostname 已与 upstream 对齐。
- 当前验证门槛：117 个 core 测试、72 个 CLI 测试、严格 clippy、
  rustfmt、rustdoc 和 `git diff --check`。
- GitHub Actions 固定 Rust 1.87.0，执行上述检查。

## 范围边界

- `--dpmemory` 是 DP 规划预算，不是包含分配器和进程固定开销的 RSS 硬上限。
- Oracle 矩阵中的“部分”表示尚未穷举所有参数组合，不是已知功能故障。
- GFF2 明确不做；项目输出 GFF3。
