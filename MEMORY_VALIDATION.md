# DP 内存验证

`--dpmemory` 是 DP 规划预算，不是进程 RSS 硬上限。RSS 还包括可执行文件、
分配器、输入、输出和运行时开销。

`scripts/validate_peak_rss.sh` 使用固定的 500×500 伪随机 DNA fixture：

1. 以 4096 MiB 运行完整 NER DP；
2. 以 `--dpmemory 0` 运行 checkpoint DP；
3. 比较 sugar、cigar、vulgar 和原子 `%P` 输出；
4. 断言 checkpoint 的 peak RSS 更低。

```sh
scripts/validate_peak_rss.sh
```

fixture 长度固定，确保 reference run 不会因输入扩大而改用 checkpoint backend。
