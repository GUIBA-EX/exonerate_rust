# exonerate-rs

[![CI](https://github.com/GUIBA-EX/exonerate_rust/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/exonerate_rust/actions/workflows/ci.yml)

Exonerate 2.4.0 的 Rust 重实现。目标是在保持得分、坐标、路径和主要
CLI 输出兼容的前提下，提供安全、可测试的精确比对与低内存执行路径。

## 支持范围

- 基础模型：`ungapped`、`ungapped:trans`、`affine:{global,bestfit,local,overlap}`。
- 生物模型：`ner`、`est2genome`、`protein2dna`、
  `protein2genome`、`coding2coding`、`coding2genome`、
  `cdna2genome`、`genome2genome`，以及两个 protein best-fit 变体。
- DNA、蛋白质、翻译、frameshift、phase 0/1/2 intron、split codon。
- sugar、cigar、vulgar、pretty alignment、RYO 和项目自有 GFF3。
- exact DP、启发式候选区域和 checkpoint traceback。

GFF2 不在项目范围内。

## 构建与检查

项目固定使用 Rust 1.87.0：

```bash
git clone --recurse-submodules https://github.com/GUIBA-EX/exonerate_rust.git
cd exonerate_rust
cargo build --release --locked
cargo test --workspace --all-targets --locked
```

本地执行与 CI 相同的检查：

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps --locked
```

## 使用

```bash
# DNA 局部比对
cargo run --release -p exonerate -- \
  --model affine:local query.fa target.fa

# 蛋白质到基因组
cargo run --release -p exonerate -- \
  --model protein2genome protein.fa genome.fa

# cDNA 到基因组
cargo run --release -p exonerate -- \
  --model cdna2genome --minintron 30 cdna.fa genome.fa

# 强制精确 DP，并使用低内存计划
cargo run --release -p exonerate -- \
  --model genome2genome --exhaustive yes --dpmemory 16 query.fa target.fa
```

CLI 默认 `--score 100`、`--subopt yes`。`-E`/`--exhaustive` 强制完整
DP；`--dpmemory` 是 DP 规划预算，不是进程 RSS 硬上限。内部坐标为
零基、半开区间；默认报告正向参考坐标。

## 文档

- [架构](exonerate.md)
- [兼容范围](COMPATIBILITY.md)
- [开发状态](DEVELOPMENT_STATUS.md)
- [命令行基准覆盖](ORACLE_MATRIX.md)

许可证：GPL-3.0-only。
