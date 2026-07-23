# 总体判断

这个项目不适合做逐文件、逐函数的 C→Rust 机械翻译。合理路线是：

> **以现有 Exonerate 可执行文件作为行为判定基准，围绕“模型中间表示 + 通用 DP 引擎 + 启发式流水线”重新设计 Rust 实现。**

Rust-Bio 应当作为序列算法、索引和基础比对组件库使用，但不能替代 Exonerate 的 C4 核心。Exonerate 的独特价值并不在 FASTA 解析或普通 Smith–Waterman，而在以下能力：

- 用状态、转移、坐标推进量和打分函数描述通用比对模型；
- 从同一模型产生精确 Viterbi、低内存 traceback 和启发式搜索；
- 支持 `est2genome`、`protein2genome`、frameshift、splice site、intron phase 等复杂语义；
- 输出包含 `M/C/G/N/5/3/I/S/F` 等生物学事件的 vulgar/GFF/RYO 格式。

C4 本身是状态、转移、打分器、scope、shadow、portal 和 span 组成的模型图；现有实现还会动态生成并编译 C 代码。
Exonerate 的精确 DP 使用 reduced-space/checkpoint traceback，启发式路径则来自 BSDP/SDP，并不保证全局最优。

因此建议采用两套执行策略：

- **compat 引擎**：严格复刻模型、默认参数、tie-breaking、坐标和输出，用于替换现有 Exonerate。
- **fast 引擎**：允许使用现代 seed、稀疏 chaining、SIMD 和并行调度；保证路径合法，但启发式结果不承诺与旧版逐字节相同。

------

# 1. 当前代码库的实际拆分

仓库把源码分为 `struct/general/sequence/comparison/database/c4/bsdp/sdp/model/hub/program/util`。

真正需要重写的核心大致是：

| 现有模块               | 职责                                                  | Rust 化判断        |
| ---------------------- | ----------------------------------------------------- | ------------------ |
| `sequence`, `database` | 序列、翻译、矩阵、splice predictor、FASTA 和索引      | 大量复用现成 crate |
| `comparison`, `seeder` | word neighbourhood、FSM/VFSM、HSP、query multiplexing | 部分复用，核心自研 |
| `c4`, `model`          | 模型图、Viterbi、traceback、codegen、全部比对模型     | 必须重写           |
| `bsdp`, `sdp`          | bounded/seeded sparse DP、启发式连接和扩展            | 必须重写           |
| `hub`, `program`       | 任务编排、best-N、输出、server、CLI                   | Rust 原生重构      |
| `util`, `ipcress`      | FASTA 工具和 in-silico PCR                            | 低风险独立迁移     |

现有模型包括 affine、EST-to-genome、protein-to-DNA、protein-to-genome、coding-to-genome、cDNA-to-genome、genome-to-genome、frameshift、phase 和 intron 等。

该 fork 仍以 2.4.0 为稳定基线，并明确说明旧多线程代码存在 data race，后续 master 已禁用多线程。 2026 年 7 月仓库补充了 GitHub Actions 构建和测试，但核心仍是 GLib、Autotools 和生成式 C4 架构。

------

# 2. 建议采用的 2026 Rust 技术栈

截至 2026 年 7 月 22 日，当前 Rust stable 是 1.97.1，Rust 2024 edition 已稳定。([Rust博客](https://blog.rust-lang.org/releases/latest/))

建议配置：

```toml
[workspace]
resolver = "3"

[workspace.package]
edition = "2024"
rust-version = "1.87"
license = "GPL-3.0-only"
```

`rust-version = "1.87"` 是一个合理起点，因为当前 Rust-Bio 4.0.1 的 MSRV 是 1.87；CI 同时测试 MSRV 和最新 stable。([Docs.rs](https://docs.rs/crate/bio/latest))

核心依赖建议如下：

- **`bio 4.x`**：alphabet、reverse complement、translation、q-gram、FM/FMD index、PSSM、简单 pairwise、sparse alignment。
- **`noodles-fasta`**：流式 FASTA 和 FAI 索引；当前版本为 0.63.0。
- **`noodles-gff`**：标准 GFF3 数据结构；旧 Exonerate 的 GFF2 输出不在当前兼容目标内。
- **`clap`**：复刻原有大量长短选项。
- **`rayon`**：query×target 或 query×database 粗粒度并行。
- **`tracing`**：结构化诊断和性能阶段统计。
- **`memmap2`**：大型数据库和持久化 seed index。
- **`serde`、`thiserror`、`bitvec`、`smallvec`**：模型配置、错误、mask 和短操作序列。

Rust-Bio 4.0.1 提供 pairwise、sparse alignment、FM/FMD、q-gram、PSSM、FASTA 等基础能力。([Docs.rs](https://docs.rs/bio/latest/bio/alignment/index.html)) `noodles` 则覆盖 FASTA、GFF3、SAM/BAM/CRAM 等格式；直接依赖所需的 format crate，通常比引入整个 meta crate更合适。([Docs.rs](https://docs.rs/crate/noodles/latest))

不建议默认引入 `rust-htslib`。Exonerate 的现有核心输入输出并不需要 BAM/CRAM，而 `rust-htslib` 会重新引入 HTSlib、C 工具链和 native build。未来确实增加 BAM/CRAM 输出时，可优先评估纯 Rust 的 noodles。([Docs.rs](https://docs.rs/crate/rust-htslib/latest))

------

# 3. Rust-Bio 能用在哪里，不能用在哪里

## 可以直接或间接使用

Rust-Bio 适合承担：

- DNA/protein alphabet 和 IUPAC 检查；
- reverse complement；
- 翻译和 ORF 辅助逻辑；
- substitution matrix 的基础表示；
- q-gram、FM/FMD index；
- exact k-mer anchor；
- affine 模型的参考实现和 differential oracle；
- fast 模式中的 sparse backbone 与 banded refinement。

Rust-Bio 的 sparse 模块从 exact k-mer matches 构建 LCSk++/SDP backbone，适合现代快速模式。([Docs.rs](https://docs.rs/bio/latest/bio/alignment/sparse/index.html))

## 不能直接替代

Rust-Bio 的标准 alignment operation 只有：

```text
Match, Subst, Del, Ins, Xclip, Yclip
```

([Docs.rs](https://docs.rs/bio/latest/bio/alignment/enum.AlignmentOperation.html))

而 Exonerate 的 vulgar/模型路径还需要 codon、NER、5′/3′ splice site、intron、split codon 和 frameshift。 因此不能把 `bio::alignment::Alignment` 作为内部权威表示。

另外，Rust-Bio sparse 默认基于 exact k-mer；Exonerate 的 `WordHood` 会基于 substitution/codon matrix 枚举达到阈值的近邻词，并通过 FSM/VFSM 对多个 query 进行 multiplexed scanning。 这部分仍需专门实现。

Rust-Bio 4.0 调整了 gap penalty 语义：长度为 1 的 gap 只收 gap-open，之后才收 gap-extension。([Docs.rs](https://docs.rs/bio/latest/bio/alignment/pairwise/index.html)) 从 Exonerate affine 模型的转移结构看，它同样是首个 gap symbol 收 open、后续循环收 extend，因此两者当前语义看起来是一致的；这是基于源码的推断，仍必须验证 tie-breaking、boundary clipping 和 traceback。

------

# 4. 推荐的 Cargo workspace

```text
exonerate-rs/
├── crates/
│   ├── exonerate-core
│   ├── exonerate-model
│   ├── exonerate-dp
│   ├── exonerate-seed
│   ├── exonerate-heuristic
│   ├── exonerate-splice
│   ├── exonerate-io
│   ├── exonerate-runtime
│   ├── exonerate-cli
│   ├── ipcress
│   └── exonerate-compat-tests
├── models/
│   ├── affine.toml
│   ├── est2genome.toml
│   └── protein2genome.toml
├── benches/
└── testdata/
```

职责边界应当严格：

- `core`：坐标、strand、sequence view、score、alignment trace。
- `model`：模型 IR、验证、模型组合、内置模型。
- `dp`：解释执行、score-only、traceback、checkpoint、band/region。
- `seed`：legacy wordhood、exact k-mer、FM/FMD、现代 minimizer/syncmer backend。
- `heuristic`：HSP、chaining、BSDP、SDP、refinement。
- `splice`：splice PSSM、intron 长度和 phase 打分。
- `io`：FASTA、annotation、sugar/cigar/vulgar/GFF/RYO。
- `runtime`：并行调度、缓存、best-N、确定性排序。
- `cli`：只负责参数和生命周期，不含比对算法。

------

# 5. C4 应如何重新设计

## 5.1 不要在 DP 内层使用动态 closure

第一版可以有解释执行器，但内部模型应使用密集 ID 和枚举式 scoring kernel：

```rust
pub struct Transition {
    pub from: StateId,
    pub to: StateId,
    pub query_advance: u32,
    pub target_advance: u32,
    pub kernel: KernelId,
    pub label: Label,
}

pub enum ScoreKernel {
    Constant(Score),
    Substitution(MatrixId),
    CodonSubstitution(MatrixId),
    SpliceSite(SpliceKind),
    IntronOpen,
    IntronClose {
        chain: Chain,
        min_len: u32,
        max_len: u32,
    },
    Frameshift,
    Phase,
}

pub enum Label {
    None,
    Match,
    Gap,
    Ner,
    Splice5,
    Splice3,
    Intron,
    SplitCodon,
    Frameshift,
}
```

不要在每个 cell 上调用 `Box<dyn Fn>`。解释执行阶段可以对 `ScoreKernel` 做 `match`；性能阶段再对内置模型做 AOT specialization。

## 5.2 模型关闭时做静态验证

对应原来的 `C4_Model_close`，Rust 版本应当验证：

- start/end 是否可达；
- epsilon transition 是否形成环；
- 每个循环是否至少推进一个坐标；
- transition 的最大 query/target advance；
- scope 和边界条件是否合法；
- portal 是否对应可 seed 的 match transition；
- span/intron 约束是否一致；
- state/transition 声明顺序是否固定。

最后一点很重要：**tie-breaking 经常隐式依赖转移迭代顺序**。只匹配最终 score 不足以保证 vulgar 输出兼容。

## 5.3 先解释执行，后 AOT，不先做 JIT

推荐顺序：

1. 通用解释执行器，作为可读、可验证的 reference engine；
2. `build.rs` 或 procedural macro 生成内置模型的 Rust kernel；
3. 对热点模型手写或生成 specialized loops；
4. 只有在确实需要运行时自定义模型时，才考虑 Cranelift JIT。

直接复刻旧 C4 的“运行时写 C 文件、调用编译器、加载 object”会继续带来部署、安全和可复现性问题，不值得保留。

------

# 6. DP 引擎设计

Exonerate 当前 Viterbi 区分 score、path、region 和 checkpoint 模式，并支持 continuation traceback。 Rust 实现应保留这四种计算模式，而不是只写一个全矩阵 traceback。

建议实现层次：

1. **scalar exhaustive reference kernel**
   完全安全 Rust，连续 `Vec<Score>`，用于正确性判定。
2. **rolling-row score kernel**
   只保留 `max_target_advance` 所需行。
3. **checkpoint traceback**
   移植原 reduced-space 行为，避免大序列完整 traceback matrix。
4. **region/banded kernel**
   仅计算 heuristic 给出的任意形状或分段 region。
5. **specialized affine kernel**
   可选接入 Rust-Bio 或 `block_aligner`。

`block_aligner 0.5.1` 提供 SSE2、AVX2、NEON 和 WASM SIMD 的 global/X-drop affine alignment，适合作为简单 affine 和 seed extension 的可选后端，但不能执行 protein-to-genome 等复杂 C4 模型。([Docs.rs](https://docs.rs/block-aligner/latest/block_aligner/))

当前稳定 Rust 的 `std::simd` 仍是 nightly-only，因此核心 crate 不应要求 nightly。稳定 scalar 实现应始终存在，SIMD 后端使用成熟 crate 或隔离的 `std::arch` 实现。([Rust文档](https://doc.rust-lang.org/std/simd/prelude/index.html))

------

# 7. Alignment 内部表示

不要直接使用普通 CIGAR。内部应继续保存“转移 + 重复次数”：

```rust
pub struct TraceRun {
    pub transition: TransitionId,
    pub repeats: u32,
}

pub struct Alignment {
    pub query_region: Region,
    pub target_region: Region,
    pub score: Score,
    pub model: ModelId,
    pub trace: Box<[TraceRun]>,
}
```

这样才能从同一条 trace 生成：

- human-readable alignment；
- sugar；
- legacy cigar；
- vulgar；
- query/target GFF；
- RYO 的 `%P...` per-transition 字段。

Exonerate 内部使用零起点、半开区间式的 inter-base coordinates，只有 human-readable 和 GFF 使用其他约定；反向链坐标也有自身规则。 因此应采用显式的坐标适配器，禁止在 formatter 中临时做零散的 `+1/-1`。

外部位置建议使用 `u64`，DP region 转换为 `usize` 前做 checked conversion；score 则保留 `i32` newtype 以匹配旧版 raw score，并实现受控的 `NEG_INF` 加法。

序列本体使用 `&[u8]`/`Arc<[u8]>`，不能假设 UTF-8。大小写必须保留，因为 soft masking 依赖大小写；mask 最好拆成独立 bitset，而不是在算法各处重复判断 ASCII case。

------

# 8. 启发式搜索必须双轨化

## Compat seed backend

为了复刻旧行为，需要移植：

1. substitution/codon word-neighbourhood 枚举；
2. query seed multiplexing；
3. FSM/VFSM 扫描；
4. saturation threshold；
5. HSP extension 和 HSP set；
6. portal-aware HSP 连接；
7. BSDP lazy bound/edge confirmation；
8. SDP drop-off、boundary 和 suboptimal traceback；
9. exact model-aware refinement。

BSDP 不是普通的 anchor chaining。它维护带 bound 的 node/edge、延迟确认 edge cost，并逐条生成高分路径。 SDP 则围绕 HSP seed 调度正反向 traceback 和 boundary。

## Fast seed backend

现代模式可以采用：

- exact q-gram/FMD；
- minimizer 或 syncmer；
- Rust-Bio sparse LCSk++/SDP backbone；
- X-drop extension；
- banded model-aware refinement。

Rust-Bio sparse 和 banded aligner非常适合生成快速 backbone，但它们基于 exact k-mer，且不保证恢复完整 Smith–Waterman 路径，所以只能作为 fast backend 或候选 region 生成器。([Docs.rs](https://docs.rs/bio/latest/bio/alignment/sparse/index.html))

最终输出路径仍应经过 Exonerate 模型的精确 DP 验证。这样 fast 模式可以改变召回候选的方式，却不会产生不符合模型的 alignment。

------

# 9. 并行模型

旧版多线程存在 data race，Rust 重写最直接的收益就是把并行边界重新定义。

推荐只在粗粒度上并行：

```text
共享只读：
  Arc<Model>
  Arc<TargetIndex>
  Arc<SubstitutionMatrices>
  Arc<SpliceModels>

每个 worker 私有：
  DpScratch
  SeedScratch
  HspArena
  TracebackBuffers
  ResultBuffer
```

调度策略：

- query 或 query block 作为 Rayon task；
- target index 只读共享；
- DP 矩阵内部初期不并行；
- worker 不直接写 stdout；
- 每个任务分配顺序号；
- 单独 aggregator 按兼容顺序合并输出；
- `bestn` 和相同分数 tie 使用明确、稳定的 comparator。

Rayon 当前提供 work-stealing data parallelism，适合这种独立 comparison 的 CPU 调度。([Docs.rs](https://docs.rs/crate/rayon/latest?utm_source=chatgpt.com))

`exonerate-server` 若重写，可用 async runtime 处理连接和背压，但 DP 仍应放在专用 CPU pool；不要把 CPU-bound alignment 直接放进 async executor。

------

# 10. 迁移顺序

## 阶段 0：冻结行为规范

先固定两个 C reference：

- `v2.4.0`：稳定用户行为；
- 当前 master：已合入的 bug fix。

建立 `COMPATIBILITY.md`，记录：

- 每个 CLI option、默认值和别名；
- model 自动选择规则；
- gap、splice、frameshift 和 intron 打分；
- 坐标和 reverse-complement 规则；
- tie-breaking；
- sugar/cigar/vulgar/GFF/RYO 的精确格式。

现有仓库虽然各模块有很多 C unit test，但 `test/exonerate` 下只有一个 shell integration test，因此必须扩充端到端 corpus。

## 阶段 1：第一个纵向切片

只实现：

- CLI 关键参数；
- FASTA streaming；
- DNA/protein matrix；
- `ungapped`；
- `affine:global/local/bestfit/overlap`；
- exhaustive；
- sugar/cigar/vulgar；
- 单线程。

这一步必须真正端到端运行，而不是先移植几十个底层容器。

## 阶段 2：通用 Model IR 和精确 C4

实现：

- scope；
- state/transition；
- portals/spans/shadows；
- layout；
- score/path/region/checkpoint；
- suboptimal alignment；
- RYO 和 GFF。

## 阶段 3：复杂生物学模型

建议顺序：

1. `est2genome`；
2. `protein2dna`；
3. `protein2genome`；
4. `coding2coding`；
5. `coding2genome`；
6. `cdna2genome`；
7. `genome2genome`。

原 intron 模型默认约束包括最小 30、最大 200000 和 opening penalty -30，并结合 splice predictor score；这些细节必须作为兼容规范，而不是仅实现 GT–AG 检查。

## 阶段 4：Legacy heuristic

依次迁移：

- WordHood；
- FSM/VFSM；
- HSPSet；
- Heuristic portal/span；
- BSDP；
- SDP；
- refinement。

这是工作量和风险最大的阶段。

## 阶段 5：并行、索引和外围工具

加入：

- deterministic Rayon 调度；
- mmap target index；
- database chunking；
- `ipcress`；
- FASTA utilities；
- legacy server protocol adapter。

## 阶段 6：性能专门化

在 benchmark 证明必要后再加入：

- AOT-generated model kernels；
- affine SIMD；
- protein substitution vectorization；
- PGO/LTO；
- fast seed backend；
- NUMA-aware database sharding。

------

# 11. 测试策略

最关键的是 differential testing，而不是普通单元测试数量。

每个生成样例同时运行 C 和 Rust：

```text
C exonerate ──> canonical result
Rust engine ──> candidate result
               ↓
        score / coordinates /
        transition trace /
        textual output comparison
```

应覆盖：

- 全部模型；
- 正反链；
- ambiguous IUPAC；
- soft masking；
- query/target clipping；
- gap 长度 1、2、N；
- split codon；
- frameshift；
- splice phase；
- intron 长度边界；
- equal-score ties；
- suboptimal overlap；
- best-N；
- RYO escape 和 per-transition 字段。

属性测试至少包括：

- trace 的 query/target advance 等于报告区间长度；
- 重新按 transition 求和得到原 score；
- exhaustive score 不低于 heuristic score；
- global 消耗全部序列；
- bestfit 消耗完整 query；
- reverse-complement 对称性；
- vulgar serialize/parse round-trip；
- 单线程和多线程输出一致。

C 可执行文件应作为测试 oracle，而不是通过细粒度 FFI 嵌进 Rust。现有接口大量使用 GLib 容器、全局 argument set 和生成代码，细粒度 FFI 会保留几乎所有旧复杂度。现有构建仍直接依赖 GLib，并通过 Autoconf 选择 compiled C4 models。

------

# 最终建议

最适合这个项目的实施原则是：

1. **Rust-Bio 做基础设施，不做 C4 替代品。**
2. **先建立兼容 oracle，再写 Rust。**
3. **先完成 affine 纵向切片，再抽象通用模型。**
4. **精确 DP 先于启发式搜索。**
5. **compat 和 fast 分离，避免性能优化破坏可复现性。**
6. **解释执行器作为永久 reference backend，AOT/SIMD 只是优化 backend。**
7. **并行只放在 comparison 粗粒度边界，并集中排序输出。**
8. **保留 GPLv3 和原作者版权声明。** 当前源码明确采用 GPL version 3。

第一个工程里程碑不应是“移植了多少 C 文件”，而应是：

> Rust 版本能够对 `ungapped` 和四种 affine model 读取真实 FASTA，以 exhaustive 模式生成与 C 版本一致的 score、坐标、vulgar 和 sugar，并在随机 differential corpus 上稳定通过。

这个纵向切片一旦成立，后续 C4、splice、BSDP/SDP 和并行化才有可靠的演进基础。
