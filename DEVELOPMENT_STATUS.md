# Development status (in progress 2026-07-23)

The full compatibility goal is **not complete** yet.

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
  - reverse-target reports now retain reverse traversal coordinates (high to low) in the generic, protein-to-DNA, protein-to-genome, and coding-to-genome full and checkpointed paths.
  - `--forwardcoordinates no` restores reverse-oriented coordinates, while the upstream-default forward-reference projection remains the default.
  - RYO `%qd` and `%td` now preserve each record's full FASTA definition line rather than aliasing the identifier.
  - RYO equivalenced `%e[tism]` and percent `%p[ciIsS]` fields now follow upstream's MATCH-only accounting, including translated codon comparisons and the BLAST-identity gap denominator.
  - RYO `%qc...` and `%tc...` now extract coding-only ranges and sequences, including codon gaps and split codons, instead of aliasing the whole aligned range.
  - All non-transition RYO sequence fields now use the upstream 70-column FASTA block layout, including its terminating newline.
  - Human-readable composed/translated reports now include the upstream `Model:` header line; affine reports retain the upstream header layout without it.
  - Protein-affine pretty reports now retain their upstream `Model: affine:…:protein2protein` header (DNA affine reports omit it).
  - Invalid `--querytype` / `--targettype` values, model/type conflicts, and repeated `--model` options now fail before FASTA processing, rather than silently selecting an incompatible scoring mode.
  - Generic affine runs now infer DNA versus protein alphabets from FASTA records when neither type option is explicit, matching upstream's default input handling.
  - Upstream byte-range input partitioning is available through `--querychunkid` / `--querychunktotal` and `--targetchunkid` / `--targetchunktotal`; boundaries advance to complete FASTA records.
  - For model families with an implemented suboptimal executor, `--bestn` now invokes the same candidate enumeration as upstream even without an explicit `--subopt`, then retains every path tied at the per-query cutoff in deterministic coordinate order.
  - Protein-to-protein `affine:local` now routes CLI `--subopt` and `--bestn` requests through its existing Waterman–Eggert executor instead of silently returning one optimal path.
  - Default `ungapped` DNA, protein, and protein-to-DNA models now enumerate pair-disjoint HSPs for `--subopt` / `--bestn`, including all paths tied at the cutoff.
  - CLI defaults now match upstream's `--score 100` and enabled `--subopt` behavior; scopes without suboptimal enumeration retain their single optimal path without failing.
  - Phase-1/2 intron split-codon exclusion now checks every participating nucleotide pair, rather than only the two boundary pairs.
  - codon matches use vulgar `C` where upstream does.
  - joint `genome2genome` intron loops remain separate target/query `I` operations instead of being incorrectly collapsed.
- The following complete sugar/cigar/vulgar/common-RYO fixture records are byte-identical to the locally built upstream 2.4.0 oracle: `protein2genome`, `coding2coding`, `coding2genome`, and `genome2genome`.
- Exact checkpointed traceback is implemented for all affine scopes and DNA/protein scoring, `est2genome`, `protein2genome`, `coding2genome`, and composed `cdna2genome` (including UTR/CDS epsilon edges, phase 0/1/2 introns, and frameshifts). It uses rolling score rows, periodic score checkpoints, and per-section parent recomputation. Frameshift models preserve the preceding five score rows at each checkpoint so section boundaries remain exact.
- `-D/--dpmemory Mb` is accepted and selects the checkpoint path when the estimated full DP exceeds the limit for affine scopes, `est2genome`, `protein2genome`, `coding2genome`, and `cdna2genome`. Exhaustive `genome2genome` now uses the same full-DP estimate to select exact block replay and feeds the requested byte budget into its checkpoint-block planner; its forced `-D 0` and budget-triggered `-D 1` regression outputs are byte-identical to the full-memory CLI run. The heuristic path retains its explicit `-D 0` low-memory selector.
- Current validation gate: 100 core tests + 64 CLI tests pass; strict workspace clippy passes with `-D warnings`; `git diff --check` passes. The CLI oracle regressions cover help/version aliases, both input chunk directions, `ungapped` DNA/protein/protein-to-DNA tied HSPs and score boundaries, all affine scopes (including protein-affine tied suboptimal traceback, default suboptimal behavior, and `bestn` multi-record ordering), NER and translated multi-record ordering and score boundaries, EST multi-record ordering and score boundaries, both `protein2dna` and `protein2genome` scopes (including translated multi-record ordering, reverse-target reports, `--forwardcoordinates no`, the protein-to-genome tie path, and reverse RYO ranges/sequences), composed-model multi-record ordering and score boundaries, `coding2coding` and `coding2genome` suboptimal ordering, `bestn` implicit-suboptimal and tie ordering, affine plus translated/composed percent-filter ordering, boolean exhaustive syntax, full FASTA RYO definitions, coding regions, equivalenced statistics, human-readable affine and composed headers, and explicit invalid-option/model-type/RYO failures; affine, EST, and all four query/target strand combinations for `coding2genome`, `cdna2genome`, and `genome2genome` (including the default heuristic and suboptimal genome2genome paths, plus forced and budget-triggered low-memory equivalence where implemented), plus complete spliced/split-codon reports for `protein2genome`, `coding2genome`, `cdna2genome`, and `genome2genome`; the latter keeps exhaustive forward-only/all-strand and heuristic all-strand `-D 0` CLI output pinned.

## Still incomplete

- CLI compatibility is not yet globally proven. The audit still needs a systematic oracle matrix for every supported model, both strands, multi-record ordering, thresholds, `bestn`, suboptimal output order, pretty alignment, errors, headers/footer, and all supported RYO tokens.
- `--showquerygff` / `--showtargetgff` are **not upstream-byte-compatible**. Rust emits project-specific GFF3 records, while upstream emits wrapped GFF2 dumps. GFF2 compatibility is explicitly out of scope and will not be implemented.
- Checkpointed low-memory traceback covers affine DNA/protein models, `est2genome`, `protein2genome`, `coding2genome`, `cdna2genome`, and `genome2genome`. The latter currently rebuilds each prefix to recreate its parent block, because query/joint intron candidates are monotonic deques spanning an unbounded number of query rows. Its saved score/queue checkpoint format is ready, but restoring it to avoid prefix recomputation remains a performance optimization rather than a correctness gap.
- Checkpoint spacing minimizes an explicit estimate covering saved score checkpoints, temporary parent blocks, rolling score rows, and (where required) five predecessor rows for frameshifts. The `genome2genome` planner includes its per-checkpoint live-intron-queue snapshot payload and receives the requested `--dpmemory` MiB budget. An infeasible budget, including zero, deterministically selects the minimum-footprint checkpoint plan. Allocation-overhead accounting and peak-RSS tests are still needed before this can be called a strict cap.
- Non-affine suboptimal enumeration now expands 3:3 codon MATCH transitions into all three blocked nucleotide pairs, fixing shifted-path re-enumeration. `coding2coding` and `coding2genome` have complete upstream order oracles; the remaining model/tie-order matrix still needs completion.
- README and `COMPATIBILITY.md` need a final documentation pass after the remaining compatibility and memory work is finished.

## Recommended continuation order

1. Finish the supported-model CLI oracle matrix and turn every fixed discrepancy into a permanent regression.
2. Generalize checkpoint traceback to the generic C4 state executor, then route the remaining complex models through equivalent checkpoint state snapshots.
3. Make checkpoint spacing obey `--dpmemory` quantitatively and verify peak RSS on large fixtures.
4. Run a requirement-by-requirement completion audit before declaring the compatibility goal complete.

The authoritative upstream executable for local oracle runs is `upstream/src/program/exonerate`. The conda executable was observed to hang in this environment and should not be used for the audit.
