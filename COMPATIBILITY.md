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

## 2026-10-11 扩展检查

原有固定回归用例通过不代表全输入兼容。新增确定性差分工具：

```sh
cargo build --locked
python3 scripts/audit_compatibility.py --output /tmp/exonerate-audit.json
```

工具使用种子 711，覆盖 16 个模型，以及 ungapped 和 4 个 affine scope 的
protein/protein 输入；包含随机短序列、插入、完全匹配和反向互补输入。
默认执行 exhaustive、subopt=no、score=1，DNA 输入分别检查 revcomp=no/yes。
每个输入分别运行普通 Rust、dpmemory=0 和仓库内 upstream 2.4.0。
发现差异时退出码为 1，并保存输入、完整参数、双方输出及低内存输出。
低内存参数对部分模型不适用，此时一致性结果不证明选择了 checkpoint backend。

首次 592 组检查发现 117 组报告差异；该历史统计与 26 个代表性反例保存在
[修复前证据](docs/audit/compatibility-counterexamples.json)。根据原版源码修复后，
同一批 592 组及另一种子 20261011 的 2,368 组均无 upstream 差异、
无普通/低内存输出差异、无执行错误。见[修复后验证](docs/audit/compatibility-validation.json)。
蛋白输入在差分工具中显式指定类型，避免短序列被原版自动识别为 DNA。

### 依据原始代码完成的修正

| 原版依据 | 修正 |
| --- | --- |
| `upstream/src/database/fastapipe.c` | translated/coding 枚举四种链组合；EST 和 NER 反转 query；EST 分别保留各链结果 |
| `upstream/src/c4/viterbi.c`、`c4.c` | 同分终点按 target 优先顺序选择；保持原版零分路径与局部重新起始的优先级 |
| `upstream/src/model/ner.c` | NER 在精确搜索中可使用短跨度；minner/maxner 是 span 提示，不是硬性 DP 边界；恢复开启、循环和关闭的原子路径 |
| `upstream/src/model/genome2genome.c` 与闭合模型转移顺序 | 保留零分 gap 接入 CDS 前的路径，防止截短同分命中 |
| `upstream/src/c4/alignment.c` | CIGAR 前导插入/删除空格与原版一致 |

永久测试保存了 120 个原版报告反例，分别以普通和低内存参数运行；
另有 NER 短跨度的原子 RYO 基准，验证三种 minner/maxner 设置。
以上结论限于已执行用例，不表示所有输入与参数组合完全兼容。

### 本轮已修复

- protein2dna 局部比对在 target 第 1、2 个偏移起始时漏掉内部 query 残基：
  query=`HYR`、target=`ACGATAATG` 的正向最高分应为 5，原实现仅为 1。
- 启发式同分 seed 区域选择依赖 HashMap 随机顺序；极大 seedpadding 会溢出。
- 正交性判定误用 target 覆盖度，且 expected_task_id 不一定指向最高分参考任务。
- 内含子和剪接位点被计入普通 gap 数。
- GFF3 反向 query 的 Target 范围倒置/越界、非 ASCII 标识符损坏、重复 ID；
  现在 CLI GFF3 使用独立 ID，并始终输出正向参考坐标，不受 forwardcoordinates 影响。
- 零长度命中生成非法 GFF3 区间，机器报告将其误记为 aligned。
- 批量失败行比 TSV 表头多一列；重复 task_id 未拒绝。
- 缓冲输出刷新失败被忽略，磁盘写满时错误返回成功。
- 显式序列类型会被后出现的 model 选项覆盖，导致参数顺序改变校验结果。
- 输出路径与输入或另一份报告相同时会覆盖文件；现在在执行前拒绝这类路径冲突。

测试同时覆盖单/双工作线程输出顺序、中文与保留字符 ID、反向坐标、
空比对、磁盘写入失败和输入文件保护。大规模真实数据性能、全部非法参数组合、
随机次优路径差分不在这 2,960 组统计内；次优路径依靠现有固定回归用例验证。
可用 `--suboptimal` 扩展差分，但本地 upstream 在部分短序列次优搜索上超过
15 秒超时，需将参考程序超时与 Rust 错误分开判断。

## 代码审查后修正

- Unix（包括 Linux/macOS）输出路径保护同时核对设备号和 inode，拒绝输入或报告的硬链接别名。
- 批量 evidence GFF3 的 ID 使用带长度的任务/目标组件，避免含点名称造成组合冲突。
- 可读报告的范围遵循 forwardcoordinates，序列提取仍使用原始参考坐标。
- 按原版 `translate.c` 处理 IUPAC 模糊密码子；结果唯一时保留氨基酸，否则为 X。
  原版生成的 5,832 种密码子基准覆盖全部 IUPAC 字符、U、X 和缺口，并核对大小写。

新增的模糊碱基差分由 150 组不一致降为零，普通/低内存输出一致，执行错误为零。
本轮重跑普通基准 592 组，以及模糊碱基基准 1,184 组；208 项回归测试通过。
详见[审查修正验证](docs/audit/review-fixes-validation.json)。新增检查可复现为：

```sh
python3 scripts/audit_compatibility.py --ambiguous --seed 42 --cases 32 --output /tmp/exonerate-iupac-audit.json
```

## 真实数据验证：优化前基线（2026-10-11）

使用 upstream 自带的 CALM、p53、hTRT 和果蝇 tube 蛋白/cDNA 原始记录，
以及 NCBI RefSeqGene [CALM1 NG_013338.1](https://www.ncbi.nlm.nih.gov/nuccore/NG_013338.1)
和 [TP53 NG_017013.2](https://www.ncbi.nlm.nih.gov/nuccore/NG_017013.2)。
下载原文、来源 URL 和 SHA-256 保存在 `docs/audit/real-data/inputs/`。
蛋白长度 149–1,132 aa，cDNA 长度 1,760–4,015 nt；基因组记录分别为
18,293 和 32,772 nt。未使用人工插入内含子的替代序列。

使用 release 构建、显式序列类型、score=100、subopt=no，检查
protein2dna、protein2dna:bestfit、protein2genome、ungapped、est2genome。
基因组任务同时检查 exhaustive=yes/no；每组运行 upstream、Rust dpmemory=64
和 Rust dpmemory=0。每个进程限定 60 秒及 8 GiB 虚拟地址空间，并关闭 core dump。

- 24 组 / 72 次执行；21 组全部完成，sugar/cigar/vulgar 与原版逐字一致，
  两种 Rust 内存配置之间也一致。
- CALM 和 p53 的精确蛋白到基因组比对均通过，包含真实剪接与拆分密码子；
  p53 命中覆盖 393 aa、9 个内含子。
- 2 组启发式 EST（subopt=no）分别申请单块 8,916,934,656 / 12,927,768,672 字节，
  在 8 GiB 上限下中止；两种 dpmemory 参数均发生。
- 1 组 p53 精确 EST 在两种 Rust 配置下均超过 60 秒，原版 7.8 秒完成。
  超时/资源失败不计作一致，也不能据此判定其最终比对结果。
- 另运行 16 个真实参考蛋白/候选 cDNA 面板任务，1/4 线程报告逐字一致，
  TSV 列数、GFF3 ID 唯一性与 Parent 引用均通过；四个候选在该面板中均为 accepted。

### 性能与阻塞项

以下为单次测量，不是统计性能基准；dpmemory 本来就不是 RSS 硬上限。

| 任务 | 原版耗时 | Rust 64 MiB 配置耗时 | 观察 |
| --- | --- | --- | --- |
| hTRT 4,015 nt ungapped 自比对 | 0.21 s | 56.99 s | 输出一致；约 271 倍耗时 |
| CALM 蛋白→基因组，精确 | 0.35 s | 4.19 s | 输出一致 |
| p53 蛋白→基因组，精确 | 1.59 s | 20.14 s | 输出一致 |
| CALM cDNA→基因组，精确 EST | 5.86 s | 42.66 s | 输出一致 |
| p53 蛋白→基因组，启发式 | 0.02 s | 15.24 s | 输出一致；Rust RSS 6.31 GiB，原版约 10.2 MiB |

优先处理：

1. EST/蛋白到基因组的启发式入口未传递 DP 内存预算，区域精修直接分配完整矩阵。
   CALM 启发式蛋白比对也达到约 1.34 GiB；EST 在上述限制下无法完成。
2. `ungapped()` 从每个格点重新扫描后续对角线，存在
   O(n × m × min(n,m)) 重复计算；这是长序列自比对的明显瓶颈。
3. 精确 EST 的运行时间仍需优化；不能将小样本兼容性通过视作生产性能达标。

完整证据见 [逐次结果](docs/audit/real-data/results.json)、
[耗时与内存表](docs/audit/real-data/results.tsv)、
[批量验证](docs/audit/real-data/batch-results.json)。可复现为：

```sh
cargo build --release --locked
python3 scripts/test_real_data.py
python3 scripts/test_real_data.py --batch-only --output docs/audit/real-data/batch-results.json
```

这是基因级试验，未覆盖全基因组吞吐、大面板、随机次优路径或所有 16 个模型。
该基线保留优化前发现的问题；后续修复与复测见下节。历史证据以 JSON 内嵌的 stdout/stderr 为准；早期 `runs/` 目录曾被中间复测复用，新版工具已按报告名隔离日志。


## 对照原版的修复与性能优化：第一阶段（2026-10-11）

本轮实现改动：

- `ungapped` 每条对角线只扫描一次，移除从每个格点重复扫描的三次方工作量；保留零分前缀、禁止匹配格点与同分终点规则。
- 参考 `upstream/src/sequence/splice.c` 预计算剪接分值；密码子采用原版 IUPAC 展开语义的查表翻译，连续密码子与蛋白替换分值在 DP 外缓存。
- EST、protein2genome 回溯仅保存紧凑前驱，在选中路径上展开报告步骤；启发式 `subopt=no` 入口传递 `dpmemory` 并使用 checkpoint 回溯。
- 当整个 target 不超过 maxintron 时，候选不可能过期，因此仅保留最佳候选；仍保留原有同分选择规则。较小 maxintron 继续使用有界队列。
- protein2dna 评分改用两行滚动缓冲，同时选择终点；protein2dna/coding2coding 缩小单步前驱记录。
- 参考 `upstream/src/model/est2genome.c` 的 `C4_Model_make_stereo` 补齐独立的反向剪接分支，避免将两种剪接方向混入同一路径。TP53 精确 EST 反向结果从错误的 1008 分恢复为原版的 1017 分；对应 sugar/cigar/vulgar 完全一致。
- RYO 展开 EST 的逐碱基 intron 循环，反向分支使用原版 transition 名称和剪接标签。新增固定反向 EST 报告及原子 RYO 基准。

209 项回归测试、格式检查、严格 Clippy 及无警告文档构建通过。release 构建的 1,184 组模糊碱基差分无原版差异、无普通/低内存差异、无执行错误，见[差分结果](docs/audit/performance-compatibility.json)。

性能复测使用同一组真实序列，每个配置运行三次，比较包含启动开销的高精度墙钟时间中位数；执行期间不同时编译或运行其他测试。结果保留每次 stdout/stderr、命令、RSS、二进制及输入 SHA-256。`--require-faster` 要求每组输出与原版一致、两种 Rust 内存模式一致，且两种模式的中位时间都小于原版，否则返回失败。三次中位数仅用于本机这批输入的验收，不构成所有输入或硬件的速度保证。

复现：

```sh
cargo build --release --locked
python3 scripts/test_real_data.py --repeat 3 --require-faster --output /tmp/exonerate-performance/results.json
```

完整复测共 24 组、216 次执行：全部完成，无超时或内存中止；22 组与原版输出一致，两种 Rust 内存配置在全部 24 组中一致。20 组精确搜索均通过严格性能条件，4 组启发式搜索未通过，**整体验收失败（退出码 1）**。

| 任务 | 原版中位耗时 | Rust 64 MiB 中位耗时 | Rust 峰值 RSS |
| --- | --- | --- | --- |
| hTRT ungapped 自比对 | 0.222 s | 0.046 s | 3.2 MiB |
| hTRT protein2dna | 0.243 s | 0.192 s | 81.3 MiB |
| CALM protein2genome 精确 | 0.364 s | 0.273 s | 11.6 MiB |
| TP53 protein2genome 精确 | 1.607 s | 1.265 s | 28.8 MiB |
| CALM EST 精确 | 5.885 s | 3.038 s | 32.6 MiB |
| TP53 EST 精确 | 7.805 s | 4.363 s | 55.3 MiB |
| TP53 protein2genome 启发式 | 0.023 s | 1.284 s | 28.0 MiB |

第一阶段证据见[逐次结果](docs/audit/real-data/performance-validated/results.json)、[结果表](docs/audit/real-data/performance-validated/results.tsv)。其中 `samples` 保留全部三次测量，表中为中位数；`dpmemory` 不是进程 RSS 上限。

16 个真实批量任务的 1/4 线程复测也通过：报告逐字一致，TSV 列数、GFF3 ID 与 Parent 引用有效，详见[批量复测](docs/audit/real-data/performance-validated/batch-results.json)。

第一阶段发现的问题（后续进展见下一节）：

- 启发式 EST 的两个真实基因组输入仍与原版不同。重新检查原版入口后确认，默认采用 `upstream/src/sdp/sdp.c` 的稀疏扩展（`bsdp/heuristic.c` 是另一条可选路径）；当前 seed 区域按 maxintron 扩展，在这些输入上覆盖整个 target，实际退化为精确搜索。内存问题已修复，但搜索算法仍需移植，不能以更高分替代原版启发式输出。
- 两个启发式 protein2genome 输入虽输出一致，耗时仍高于原版。需要原版默认 HSP 边界与稀疏搜索，继续压缩矩阵或评分缓存不能弥补搜索量差距。
- 所有上述真实数据测试使用 `subopt=no`；`subopt=yes` 专用路径的完整矩阵内存问题、其他模型的大数据性能和全基因组吞吐尚未验收。
- 当前不满足“全部功能、全部输入只能比原版快”的要求；严格验收失败必须保留，不能将已通过的精确搜索用例泛化为全面达标。


## 默认 SDP 启发式修复：第二阶段（2026-10-11）

重新沿 `upstream/src/hub/gam.c` 的默认入口核对，实际启用的是 `sdp/sdp.c`，而非可选 BSDP。现为默认配置、`subopt=no` 的 EST 和局部 protein2genome 实现对应的稀疏扩展：

- HSP 采用原版 DNA 12 字长、蛋白 6 字长及 BLOSUM 邻域，按累计分值选择中心；无合格 HSP 时不再返回精确搜索产生的额外命中。
- 从 HSP 中心进行反向扩展，跳过 intron span 循环并记录起点／span 边界；正向仅在这些边界恢复冻结的 intron 候选，保留原版 50 分下降阈值和反向 shadow 计分语义。
- EST 保持独立的正反剪接分支；protein2genome 覆盖三个剪接相位、拆分密码子、移码与仿射间隙。采用滚动评分行、活动位置位图、延迟剪接评分和稀疏回溯，避免 seed 区域扩展成整个矩阵。
- 按 `analysis.c` / `HSPset_revcomp` 保留反向 EST 转换前的 HSP 分值中心，修复扩大随机审计后发现的额外边界反例。
- 原子 RYO 差分还发现并修复了 protein2genome 移码名称：原版继承 protein2dna 的 `p2d` 后缀。

固定回归增加 25 个原版生成案例，涵盖剪接相位、反向 target、移码、缺失和较小 maxintron，同时比较 sugar/cigar/vulgar 及逐步骤 RYO，两种内存配置均验证。

214 项回归测试、格式检查、严格 Clippy 和无警告文档构建通过。新增 800 组启发式差分同时比较三个汇总报告及原子 RYO，覆盖 forcegtag 和 maxintron 选项，零原版差异、零内存配置差异、零执行错误；见 [逐步骤差分结果](docs/audit/seeded-raw-validated.json)。发现问题时的 [原始报告](docs/audit/seeded-raw-compatibility.json) 保留用于追溯。

复现差分：

```sh
python3 scripts/audit_compatibility.py --rust target/release/exonerate-rs --heuristic --raw-steps --cases 200 --output /tmp/seeded-audit.json
```

真实数据完整复测：24 组、216 次执行，全部完成且输出与原版一致；两种 Rust 内存配置的三次中位耗时均低于原版，严格验收通过。证据见 [SDP 完整结果](docs/audit/real-data/sdp-final/results.json) 和 [结果表](docs/audit/real-data/sdp-final/results.tsv)。第一阶段失败记录保留，不能与本阶段混用。

| 启发式任务 | 原版中位耗时 | Rust 64 MiB | 加速比 |
| --- | --- | --- | --- |
| calm.human.genomic.protein2genome.exhaustive-no | 0.007973 s | 0.005641 s | 1.41× |
| calm.human.genomic.est2genome.exhaustive-no | 0.013218 s | 0.008859 s | 1.49× |
| p53.human.genomic.protein2genome.exhaustive-no | 0.023092 s | 0.015249 s | 1.51× |
| p53.human.genomic.est2genome.exhaustive-no | 0.079125 s | 0.032060 s | 2.47× |

16 个真实批量任务的 1/4 线程报告也逐字一致，TSV 与 GFF3 校验通过；见 [最终批量结果](docs/audit/real-data/sdp-final/batch-results.json)。

仍未验收的范围包括 `subopt=yes` 的大数据内存、其他模型的大规模性能及全基因组吞吐。自定义 seed 配置、protein2genome:bestfit 和非标准蛋白 query 字母仍保留原有后端；上述 SDP 兼容性结论不能推广至这些分支。稀疏回溯保存访问过的路径，`dpmemory` 不是进程内存硬上限，也不代表全基因组内存已经达标。速度结论限于本机、本次输入和测量方式，不保证所有输入都比原版快。


## protein2dna 回溯压缩：第三阶段（2026-10-11）

局部及 bestfit 蛋白／DNA 比对的三个前驱矩阵改用单字节编码；编码保存转移类型和来源状态，回溯时还原操作及推进长度。评分递推、同分选择和报告路径保持原有逻辑。本次没有改动 coding2coding 使用的前驱表示。

后台验证完成：214 项回归测试、严格 Clippy、格式检查通过；1,184 组包含模糊碱基的差分发现 1 个 genome2genome 反向同分路径报告差异，已用修改前二进制复现，属于既有问题；本次涉及的 protein2dna / bestfit 无差异，所有配置均无两种内存配置差异、无执行错误，见 [差分报告](docs/audit/packed-trace-compatibility.json)。

8 个真实 protein2dna / bestfit 案例，每个配置重复五次，共 120 次执行；全部输出一致，普通及低内存配置均快于原版，见 [原版对比](docs/audit/real-data/packed-trace/results.json)。另外对本次修改前后可执行文件交替运行五轮，共 80 次执行，8 组全部结果一致且中位耗时降低；保存二进制哈希、每次命令与结果，见 [修改前后交替对比](docs/audit/real-data/packed-trace/paired.json)。

hTRT bestfit 修改前后中位耗时从 0.218775 s 降至 0.180416 s，降低约 17.5%；峰值 RSS 的中位数从 83,104 KiB 降至 16,776 KiB，降低约 79.8%。单独的原版对比中，原版为 0.225049 s，Rust 64 MiB 为 0.186163 s。两轮测量分别报告，不混用其时间计算收益。8 组修改前后加速比为 1.13–1.21×；结论仍限于这批数据和本机测量。


## genome2genome 同分路径与逐步骤报告修复（v0.2.1 之后）

修复第三阶段保留的反向同分路径差异。原版 `C4_Model_close` 的 transition 顺序先处理 post-CDS 间隙返回，再处理 CDS-to-UTR 空转移；本实现现在保留该优先级。原版与 Rust 均得到 204 分，反向路径统一为 `M 27 27 C 6 6 G 0 2 M 27 27`。

逐步骤 RYO 同步修正：post-CDS 入口使用 `start to match` 名称，补齐间隙闭合步骤，删除重复终止步骤及空 pre-CDS 前缀；checkpoint 按实际终止步骤数截取回放片段，避免丢失末尾匹配。四个原版反例固定为回归案例，在 0/64 MiB 配置下比较完整报告。

验证：216 项回归测试、严格 Clippy、格式及文档检查通过。原来的 1,184 组差分现在零原版差异、零内存配置差异、零执行错误，见 [完整复测](docs/audit/genome-tie-validation.json)。另增加 400 组 genome2genome 原子 RYO 差分，全部一致，见 [逐步骤复测](docs/audit/genome-raw-validated.json)。保留 [修复前逐步骤反例](docs/audit/genome-raw-validation.json) 供追溯。

```sh
python3 scripts/audit_compatibility.py --rust target/release/exonerate-rs --cases 32 --ambiguous --output /tmp/genome-tie.json
python3 scripts/audit_compatibility.py --rust target/release/exonerate-rs --models genome2genome --cases 200 --ambiguous --raw-steps --output /tmp/genome-raw.json
```

以上修复位于 v0.2.1 之后，已发布的 v0.2.1 内容不变。此次是正确性修复，没有新增全基因组性能保证。
