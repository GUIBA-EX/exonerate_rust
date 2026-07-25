# exonerate-rs

[![CI](https://github.com/GUIBA-EX/exonerate_rust/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/exonerate_rust/actions/workflows/ci.yml)

<p align="center">
  <img src="docs/images/ferris-flat.svg" alt="Ferris，Rust 社区吉祥物的修订版平面原图" width="480">
  <br>
  <sub>Ferris by Karen Rustad Tölva · <a href="https://rustacean.net/">rustacean.net</a> · CC0</sub>
</p>

Exonerate 2.4.0 的 Rust 重实现。项目以得分、traceback、坐标和主要 CLI
输出兼容为目标，同时提供精确的低内存执行路径。

## 快速开始

需要 Rust 1.87.0；仓库中的 `rust-toolchain.toml` 会自动选择该版本。

```bash
git clone --recurse-submodules https://github.com/GUIBA-EX/exonerate_rust.git
cd exonerate_rust
cargo build --release --locked

target/release/exonerate-rs --help
target/release/exonerate-rs \
  --model affine:local \
  --verbose 0 \
  query.fa target.fa
```

输入为 FASTA 文件。使用 `-q QUERY.fa -t TARGET.fa` 与使用最后两个位置参数
等价。

## 支持的模型

```text
ungapped             ungapped:trans
affine:global        affine:bestfit
affine:local         affine:overlap
ner                  est2genome
protein2dna          protein2dna:bestfit
protein2genome       protein2genome:bestfit
coding2coding        coding2genome
cdna2genome          genome2genome
```

示例：

```bash
# 蛋白质到基因组
target/release/exonerate-rs \
  --model protein2genome protein.fa genome.fa

# cDNA 到基因组
target/release/exonerate-rs \
  --model cdna2genome --minintron 30 cdna.fa genome.fa

# 完整搜索，并在需要时使用 checkpoint traceback
target/release/exonerate-rs \
  --model genome2genome \
  --exhaustive yes \
  --dpmemory 16 \
  query.fa target.fa
```

## 行为约定

- 默认值与 upstream 一致：`--verbose 1`、`--score 100`、`--subopt yes`。
- `-E` 或 `--exhaustive yes` 跳过启发式候选区域，执行完整搜索。
- `--dpmemory` 是 DP 规划预算，不是整个进程的 RSS 硬上限；在支持
  checkpoint 的模型中，设为 `0` 可强制走低内存路径。
- 支持 sugar、cigar、vulgar、pretty alignment、RYO 和项目自有 GFF3。
- 内部坐标为零基半开区间；默认输出正向参考坐标。
- GFF2 不在项目范围内。

## 批量候选审计

`--tasks` 将多个显式任务并发执行；每行固定选择一个 FASTA 记录，适合把
上游候选序列的选择和比对证据分开保存。结果仍按清单顺序写出，不随线程数改变。

```text
task_id	model	query_fasta	query_id	target_fasta	target_id
gene_001	protein2genome	proteins.fa	gene_001	candidates.fa	contig_42
```

```bash
target/release/exonerate-rs \
  --tasks tasks.tsv --threads 4 \
  --audit protein-candidate \
  --result-tsv evidence.tsv \
  --evidence-gff3 evidence.gff3
```

`--audit protein-candidate` 是 `protein2dna` / `protein2genome` 的薄预设：关闭
人读报告、执行完整搜索，并要求 `--result-tsv`。TSV 包含覆盖度、缺口、移码和
内含子计数；GFF3 是比对证据（`match` / `match_part`），不是基因结构预测。

### 参考面板正交性审计

在蛋白候选审计中，可使用扩展任务清单让同一候选与多个 family 竞争：

```text
task_id	sample_id	locus_id	candidate_id	expected_family	reference_family	model	query_fasta	query_id	target_fasta	target_id
S1.L1.expected	S1	Locus_001	contig_42	Fam_001	Fam_001	protein2genome	proteins.fa	Fam_001	candidates.fa	contig_42
S1.L1.other	S1	Locus_001	contig_42	Fam_001	Fam_014	protein2genome	proteins.fa	Fam_014	candidates.fa	contig_42
```

```bash
target/release/exonerate-rs \
  --tasks orthology_tasks.tsv --audit protein-candidate --threads 4 \
  --result-tsv alignments.tsv --orthology-report orthology.tsv \
  --orthology-min-coverage 0.80 --orthology-min-delta 20
```

`orthology.tsv` 每个候选一行，保留预期与最强竞争 family 的分数、差值、覆盖度、
移码、内含子数和可回溯的 task ID。状态只描述给定面板中的证据：`accepted`、
`ambiguous`、`fragment`、`unexpected_family`、`no_hit` 或 `failed`；不用于声明
真实拷贝数或取代 gene tree。

完整的模型和参数兼容证据见[兼容范围](COMPATIBILITY.md)与
[命令行基准覆盖](ORACLE_MATRIX.md)。

## 开发与验证

CI 执行以下检查：

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps --locked
```

Linux 上可额外比较 full/checkpoint 的输出和 peak RSS：

```bash
scripts/validate_peak_rss.sh
```

## 文档

- [兼容范围](COMPATIBILITY.md)：实现边界与兼容优先级
- [命令行基准覆盖](ORACLE_MATRIX.md)：各模型的 upstream 证据
- [架构](exonerate.md)：工作区、数据流和设计不变量
- [通用 C4 checkpoint 设计](C4_CHECKPOINT_DESIGN.md)：低内存回放设计
- [DP 内存验证](MEMORY_VALIDATION.md)：预算语义与 RSS 验证
- [引用信息](CITATION.cff)：软件引用元数据

许可证：[GPL-3.0-only](LICENSE)。
