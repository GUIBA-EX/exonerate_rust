# 开发状态

更新时间：2026-07-23。完整兼容目标尚未完成。

## 已完成

- upstream CLI 列出的 16 个模型均有 Rust 执行路径。
- 通用 C4 状态图支持有限状态、epsilon、NER、query/target/joint intron
  和 phase 0/1/2 split codon。
- 复杂模型支持 codon gap、双侧 frameshift、splice PSSM 与原子 traceback。
- 主要模型支持 Waterman–Eggert pair exclusion、`--subopt` 和 `--bestn`。
- affine、EST、protein/coding/cDNA/genome 模型已接入 checkpoint traceback。
- DNA 与翻译模型具有启发式候选区域，并可用 `--exhaustive` 回退。
- 当前验证门槛：100 个 core 测试、64 个 CLI 测试、严格 clippy、
  rustfmt、rustdoc 和 `git diff --check`。
- GitHub Actions 固定 Rust 1.87.0，执行上述检查。

## 未完成

- DNA/protein `affine:bestfit` 与 `affine:overlap` 的完整
  suboptimal/`bestn` 枚举。
- `affine:global` 的 suboptimal 行为确认。
- 通用 C4 executor 的 checkpoint traceback。
- `genome2genome` checkpoint 恢复仍会重复计算前缀。
- `--dpmemory` 的分配器开销和 peak-RSS 验证。
- 部分模型的 pretty、RYO、错误、页眉/页脚基准。

GFF2 明确不做。

## 后续顺序

1. 完成模型级 suboptimal/`bestn` 功能。
2. 泛化 C4 checkpoint，并消除 genome2genome 前缀重算。
3. 验证内存预算和 peak RSS。
4. 补齐命令行基准，完成发布审计。
