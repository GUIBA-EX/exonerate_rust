# Development status (paused 2026-07-23)

Development is intentionally paused at this checkpoint. The full compatibility goal is **not complete** yet.

## Completed and verified

- Exact scored atomic `%P` traceback was added for the implemented complex biological models, including target/query/joint introns, phase 0/1/2 split codons, bilateral frameshifts, epsilon transitions, and composed `cdna2genome` / `genome2genome` state changes.
- Exhaustive `%P` transition oracles match upstream line-for-line for the composed fixtures:
  - `cdna2genome`: 430 transitions.
  - `genome2genome`: 236 transitions.
- Waterman–Eggert-style equivalenced-pair exclusion now provides suboptimal enumeration for the implemented non-affine models: `coding2coding`, `protein2dna`, `est2genome`, `protein2genome` (local and bestfit), `coding2genome`, `cdna2genome`, `genome2genome`, `ner`, and `ungapped:trans`.
- CLI oracle fixes completed during the current audit:
  - RYO `%C` leading-space compatibility.
  - upstream `%m` value `protein2genome:local`.
  - `coding2genome` query strand is reported as `+`.
  - codon matches use vulgar `C` where upstream does.
  - joint `genome2genome` intron loops remain separate target/query `I` operations instead of being incorrectly collapsed.
- The following complete sugar/cigar/vulgar/common-RYO fixture records are byte-identical to the locally built upstream 2.4.0 oracle: `protein2genome`, `coding2coding`, `coding2genome`, and `genome2genome`.
- Exact checkpointed traceback is implemented for all affine scopes and DNA/protein scoring. It uses rolling score rows, periodic score checkpoints, and per-section parent recomputation.
- `-D/--dpmemory Mb` is accepted and selects the affine checkpoint path when the estimated full affine DP exceeds the limit. A forced `--dpmemory 0` run is byte-identical to the full-memory CLI run on the regression fixture.
- Current validation gate: 72 core tests + 3 CLI tests pass; strict workspace clippy passes with `-D warnings`; `git diff --check` passes.

## Still incomplete

- CLI compatibility is not yet globally proven. The audit still needs a systematic oracle matrix for every supported model, both strands, multi-record ordering, thresholds, `bestn`, suboptimal output order, pretty alignment, errors, headers/footer, and all supported RYO tokens.
- `--showquerygff` / `--showtargetgff` are **not upstream-byte-compatible**. Rust currently emits project-specific GFF3 records, while upstream emits wrapped GFF2 dumps with model-specific similarity/gene/exon/intron attributes. This needs a separate upstream-compatible renderer; the existing GFF3 API should not be mislabeled as equivalent.
- Checkpointed low-memory traceback currently covers affine DNA/protein models only. Complex intron/split-codon/frameshift models still allocate full parent matrices and must be migrated to checkpoint recomputation before the low-memory goal is complete.
- The current affine checkpoint section size is chosen near `sqrt(query_len)`. `--dpmemory` decides when to enter checkpoint mode, but the section/checkpoint layout is not yet derived from the exact requested MiB cap. Memory-limit accounting and peak-RSS tests remain to be added.
- Non-affine suboptimal enumeration has regression coverage, but full upstream path-order/tie-order oracle coverage remains to be completed.
- README and `COMPATIBILITY.md` need a final documentation pass after the remaining compatibility and memory work is finished.

## Recommended continuation order

1. Add upstream-compatible GFF2 and process wrapper output, while keeping explicit GFF3 output as a separate option/API.
2. Finish the supported-model CLI oracle matrix and turn every fixed discrepancy into a permanent regression.
3. Generalize checkpoint traceback to the generic C4 state executor, then route specialized complex models through equivalent checkpoint state snapshots.
4. Make checkpoint spacing obey `--dpmemory` quantitatively and verify peak RSS on large fixtures.
5. Run a requirement-by-requirement completion audit before declaring the compatibility goal complete.

The authoritative upstream executable for local oracle runs is `upstream/src/program/exonerate`. The conda executable was observed to hang in this environment and should not be used for the audit.
