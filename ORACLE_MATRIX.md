# 命令行基准覆盖

`是`表示已有永久 upstream 2.4.0 基准；`部分`表示代表性覆盖；`否`表示暂无基准。

这些标记仅表示固定用例覆盖，不表示所有输入兼容。2026-10-11 扩展差分
检查发现的链枚举和同分路径问题已经修复，详见[兼容性验证](COMPATIBILITY.md#2026-10-11-扩展检查)。

| 模型 | 基础输出 | 链方向 | 多记录 | 阈值/Best-N | 次优路径 | 低内存 |
| --- | --- | --- | --- | --- | --- | --- |
| `ungapped` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `ungapped:trans` | 是 | 部分 | 是 | 是 | 是 | 是 |
| `affine:global` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `affine:bestfit` | 是 | 部分 | 否 | 部分 | 是 | 是 |
| `affine:local` | 是 | 部分 | 是 | 是 | 是 | 是 |
| `affine:overlap` | 是 | 部分 | 否 | 否 | 是 | 是 |
| `protein2dna` | 是 | 是 | 否 | 部分 | 部分 | 不适用 |
| `protein2dna:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 不适用 |
| `est2genome` | 是 | 部分 | 是 | 是 | 部分 | 是 |
| `protein2genome` | 是 | 是 | 是 | 部分 | 部分 | 是 |
| `protein2genome:bestfit` | 是 | 部分 | 否 | 否 | 部分 | 部分 |
| `coding2coding` | 是 | 是 | 否 | 否 | 是 | 不适用 |
| `coding2genome` | 是 | 是 | 是 | 是 | 是 | 是 |
| `cdna2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `genome2genome` | 是 | 是 | 是 | 是 | 部分 | 是 |
| `ner` | 是 | 部分 | 是 | 是 | 是 | 是 |

## 跨模型检查

| 功能 | 状态 |
| --- | --- |
| help、version、别名 | 是 |
| FASTA byte chunk | 是 |
| 反序 ID 的多记录排序 | 是 |
| 默认 score/subopt | 是 |
| 非法参数与类型冲突 | 部分 |
| pretty 与通用 RYO | 是 |
| 页眉/页脚 | 是 |
| GFF2 | 不在范围内 |

低内存基准包含零分 NER 和 splice-boundary `genome2genome` 的原子 `%P`
transition 等价性。

固定基准新增 120 个原版差分反例，以及 NER 短跨度原子路径/跨度提示基准。
种子 711 的 592 组与种子 20261011 的 2,368 组精确比对均无差异；
范围、复现命令和限制见兼容性文档。

真实生物序列试验覆盖 5 个模型、24 组配置，以及 16 个批量候选任务。
优化前 21 组完成并与原版一致，3 组受内存/时间限制未完成；后续优化修复资源失败与 EST 反向剪接，启发式差异和严格性能验收状态见兼容性文档。
另增加反向 EST 的 sugar/cigar/vulgar 和逐原子 RYO 原版基准，普通/低内存后端均验证。

默认 SDP 启发式另固定 25 个原版案例，覆盖剪接相位、移码、反向 target 和 maxintron 边界；包含原子 RYO 名称／分值，0/64 MiB 配置均验证。
后续 24 组真实数据严格验收已通过，完整证据见 `docs/audit/real-data/sdp-final/results.json`；限制见兼容性文档第二阶段。

后续补齐 genome2genome 编码区前后边界和反向同分路径的 4 个固定原版反例，比较完整 sugar/cigar/vulgar 及原子 RYO；0/64 MiB 配置均验证。1,184 组全模型差分和 400 组 genome2genome 原子差分均通过，详见兼容性文档后续修复章节。
