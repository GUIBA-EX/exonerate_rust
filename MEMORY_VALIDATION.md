# DP 内存验证

`--dpmemory` 是 DP 算法选择和规划预算，不是进程 RSS 硬上限。进程峰值还包含
可执行文件、分配器元数据、输入、输出和 checkpoint 队列。

仓库提供 `scripts/validate_peak_rss.sh`，它会：

1. 用固定种子的伪随机序列构造 500×500 确定性 DNA fixture，避免只验证
   候选队列容易退化的简单重复序列；
2. 用已知足以容纳该 fixture 完整 DP 表的 4096 MiB 计划和
   `--dpmemory 0` checkpoint 分别运行 NER；
3. 用 `/usr/bin/time` 读取 maximum resident set size；
4. 验证两次 sugar/cigar/vulgar 和原子 `%P` transition 输出逐字节相同，
   且 checkpoint 峰值低于全表。

脚本不接受自定义序列长度，避免所谓 full reference 因输入变大而静默切换到
checkpoint backend。

运行：

```sh
scripts/validate_peak_rss.sh
```

2026-07-26 在当前 Linux x86-64 开发环境、每条序列 500 bases 的结果为：

| 后端 | peak RSS |
| --- | ---: |
| full generic table (`--dpmemory 4096`) | 331608 KiB |
| generic checkpoint (`--dpmemory 0`) | 68700 KiB |

绝对值会随工具链、链接器和分配器变化，因此自动验证只断言完整报告相同和
checkpoint RSS 较低，不固定某个机器相关阈值。
