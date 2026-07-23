# Compatibility status

The upstream C executable in `upstream/` is the behavioural oracle.  This Rust
workspace currently provides the first end-to-end compatibility slice.

| Surface | Status |
| --- | --- |
| FASTA records (streaming parser, first identifier token) | implemented |
| DNA scoring (upstream 24×24 `nucleic` matrix, including IUPAC) | implemented |
| `ungapped` exhaustive local search | implemented |
| `affine:global`, `affine:bestfit`, `affine:local`, `affine:overlap` | implemented |
| `--gapopen`, `--gapextend`, deterministic tie precedence | implemented |
| forward and reverse-complement query strand | implemented |
| `sugar:`, `cigar:`, `vulgar:`, human-readable `--showalignment`, and query/target GFF3 output | implemented |
| protein scoring (upstream BLOSUM62), protein:protein affine alignment | implemented |
| protein:dna direct codon-affine alignment, 1/2/4/5-nt frameshifts | implemented (no intron) |
| splice predictor (upstream primate PSSM) | implemented |
| `est2genome` local affine/intron Viterbi, target-strand stereo, canonical `5/I/3` traceback and CLI | implemented; upstream score and 44 nt/238 nt trace fixture verified |
| `protein2genome` phase 0/1/2 intron Viterbi, split-codon traceback, target-strand stereo and CLI | implemented; upstream 125-point phase-intron oracle, coordinates and vulgar trace verified |
| `coding2coding` codon-affine, dual-side frameshift Viterbi and traceback | implemented as library and CLI API; upstream 169-point oracle and coordinate-conserving frameshift trace verified |
| `coding2genome` coding-DNA/genome Viterbi, target phase 0/1/2 introns, dual-side frameshifts and CLI | implemented; upstream 194-point phase-2 intron oracle and full-span coordinates verified |
| `cdna2genome` unified 5′ UTR/CDS/3′ UTR state lattice, target introns, CDS phases/frameshifts and CLI | implemented; upstream 1281-point composed-model oracle and 270 nt/432 nt traceback coordinates verified |
| `genome2genome` UTR and CDS query/target/joint introns, CDS phase 0/1/2 split codons, codon gaps, and bilateral frameshifts | implemented; upstream 557-point joint-intron oracle plus phase-specific query/joint trace and coordinate-conservation regressions verified |
| generic C4 Viterbi/traceback with finite states and bounded query/target/joint long states | implemented; scored atomic raw transitions now include epsilon, NER, and query/target/joint introns; affine short-sequence, protein2dna, intron, and phase oracles verified |
| `ner` bounded 2-D spans | implemented; upstream 208-point oracle and `N` traceback verified |
| `ungapped:trans` | implemented; upstream 22-point codon oracle and CLI aliases verified |
| `protein2genome:bestfit` | implemented; upstream C/Rust -8 score, full-query boundaries, and vulgar path verified |
| common `--ryo` identifiers, ranges, score, rank, model, sugar/cigar/vulgar blocks, escapes, and scored raw `%P...` transitions | implemented for current models; composed cDNA/genome and genome/genome `%P` transition oracles are exact |
| `--subopt yes` equivalenced-pair exclusion | implemented for affine local and current non-affine models; full upstream tie-order audit remains open |
| DNA `affine:local` k-mer seed/diagonal cluster/exact refinement | implemented with exhaustive fallback and `-E`; exact 800-nt planted path verified, measured 10x wall-time and 7.9x RSS improvement on the checked fixture |
| generic target, query, and joint split-codon phase shadows | implemented; exact phase-1 query/joint scoring, traceback, and coordinate regressions verified |
| intron-aware seed-region refinement for `est2genome`, `coding2genome`, `cdna2genome`, and `genome2genome` | implemented with seedless exhaustive fallback and `-E`; planted 200-nt intron fixtures match exhaustive vulgar output |
| translated three-frame-per-strand seeding for `protein2genome` and `protein2genome:bestfit` | implemented; forward, reverse, intron, local, and bestfit outputs verified against exhaustive DP |
| affine checkpoint traceback and `-D`/`--dpmemory` dispatch | implemented for DNA/protein affine scopes; complex-model checkpointing and exact MiB layout remain open |

Coordinates are zero-based, half-open inter-base coordinates internally and in
the report fields, matching Exonerate sugar/vulgar conventions.  The report
writer is deliberately separate from the trace representation so that later
biological labels can be added without changing the DP API.
