# CLI oracle matrix

This matrix tracks byte-for-byte comparisons against the locally built
upstream 2.4.0 executable at `upstream/src/program/exonerate`. It is a coverage
inventory, not a claim that a model is globally compatible.

Legend:

- **yes**: a permanent CLI regression pins the upstream result.
- **partial**: at least one case is pinned, but important variants remain.
- **no**: no permanent CLI oracle covers this dimension yet.
- **n/a**: the dimension does not apply to the current implementation.

| Model | Base reports | Strands | Multi-record / order | Thresholds / `bestn` | Suboptimal order | Pretty / RYO | `-D` equivalence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ungapped` | yes | no | no | partial | yes | no | n/a |
| `ungapped:trans` | yes | partial | yes | yes | no | partial | n/a |
| `affine:global` | yes | partial | no | no | n/a | partial | partial |
| `affine:bestfit` | yes | partial | no | no | n/a | partial | partial |
| `affine:local` | yes | partial | yes | yes | yes | partial | partial |
| `affine:overlap` | yes | partial | no | no | n/a | partial | partial |
| `protein2dna` | yes | yes | no | no | partial | partial | n/a |
| `protein2dna:bestfit` | yes | partial | no | no | partial | partial | n/a |
| `est2genome` | yes | partial | yes | yes | no | partial | yes |
| `protein2genome` | yes | yes | yes | partial | partial | yes | yes |
| `protein2genome:bestfit` | yes | partial | no | no | partial | partial | partial |
| `coding2coding` | yes | n/a | no | no | yes | partial | n/a |
| `coding2genome` | yes | yes | yes | yes | yes | yes | yes |
| `cdna2genome` | yes | yes | yes | yes | no | yes | yes |
| `genome2genome` | yes | yes | yes | yes | partial | yes | yes |
| `ner` | yes | partial | yes | yes | no | partial | n/a |

## Cross-model CLI surfaces

| Surface | Coverage | Permanent regression |
| --- | --- | --- |
| Help, version, and aliases | yes | `help_and_version_aliases_are_available_without_fasta_inputs` |
| FASTA byte chunking | yes | `fasta_byte_chunks_match_the_upstream_record_boundaries` |
| Reversed-ID multi-record `bestn` order | yes | `bestn_multi_record_order_matches_the_upstream_oracle` |
| Translated multi-record `bestn` order | yes | `protein2genome_bestn_multi_record_order_matches_the_upstream_oracle` |
| Composed multi-record order and score boundaries | yes | `composed_multi_record_order_and_score_boundaries_match_the_upstream_oracle` |
| EST multi-record order and score boundary | yes | `est2genome_multi_record_order_and_score_boundary_match_the_upstream_oracle` |
| NER and translated multi-record order | yes | `ner_and_ungapped_trans_multi_record_order_match_the_upstream_oracle` |
| Invalid model/type combinations | partial | `invalid_cli_options_and_ryo_tokens_fail_explicitly` |
| Default score and suboptimal mode | partial | `subopt_is_enabled_by_default_like_the_upstream_cli` |
| Boolean syntax | partial | `exhaustive_accepts_upstream_boolean_syntax` |
| GFF2 wrapper output | out of scope | current output remains project-specific GFF3 |
| Headers and footer | no | none |

## Next oracle batches

1. Exact score and percent boundary cases for the remaining scoring families.
2. `bestn` cutoff ties and suboptimal ordering for `est2genome`,
   `protein2genome`, `cdna2genome`, `genome2genome`, `ner`, and
   `ungapped:trans`.
3. Pretty alignment and common RYO fields for every model family.
4. Upstream process headers/footer. GFF2 compatibility is explicitly out of
   scope.
