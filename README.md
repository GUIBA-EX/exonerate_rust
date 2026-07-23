# exonerate-rs

GPL-3.0-only Rust reimplementation of Exonerate.  The current vertical slice
implements exhaustive DNA and protein `ungapped` and `affine:{global,bestfit,local,overlap}`
alignment for FASTA input, with deterministic traceback and sugar/cigar/vulgar/GFF3 reports.  Protein queries can also be searched directly against DNA (`--querytype protein
--targettype dna`) with codon-affine gaps and target-side frameshift transitions;
`protein2dna` is the intron-free model.  `protein2genome` is available for protein-to-genome local alignment with phase 0/1/2 introns, the upstream splice PSSM, and canonical split-codon/intron traceback.  `est2genome` is also available for local DNA-to-genome alignment with affine gaps, bounded target introns, upstream splice PSSMs, both target strands, and canonical intron traceback. A generic C4 graph executor is available as a library API, including bounded query-, target-, and joint-intron long states and target phase 0/1/2 codon shadows. NER spans are supported, as are `ungapped:trans` and `protein2genome:bestfit`. Target, query, and joint split-codon shadows are implemented in the generic executor. DNA `affine:local` uses k-mer seed/diagonal-cluster/exact-refine search by default. DNA genomic models (`est2genome`, `coding2genome`, `cdna2genome`, and `genome2genome`) use intron-aware padded target regions before unchanged exact DP. `protein2genome` and `protein2genome:bestfit` use translated amino-acid k-mers across all three frames of each requested target strand. `-E` selects exhaustive DP.

Clone with the upstream behavioural oracle used by the regression tests:

```bash
git clone --recurse-submodules https://github.com/GUIBA-EX/exonerate_rust.git
cd exonerate_rust
cargo test --workspace
```

```bash
cargo run --release -p exonerate -- --model affine:local query.fa target.fa

# direct protein-to-DNA alignment (codon gaps and frameshifts)
cargo run --release -p exonerate -- --model affine:local \
  --querytype protein --targettype dna protein.fa target-dna.fa

# cDNA-to-genome alignment (5′ UTR, coding region, 3′ UTR in one DP lattice)
cargo run --release -p exonerate -- --model cdna2genome \
  --minintron 30 --maxintron 200000 cdna.fa genome.fa

# coding-DNA-to-genome alignment (codon gaps, frameshifts, phase 0/1/2 introns)
cargo run --release -p exonerate -- --model coding2genome \
  --minintron 30 --maxintron 200000 coding.fa genome.fa

# intron-aware protein-to-genome alignment (phase 0/1/2 introns)
cargo run --release -p exonerate -- --model protein2genome \
  --minintron 30 --maxintron 200000 protein.fa genome.fa

# local EST-to-genome alignment; intron controls are optional
cargo run --release -p exonerate -- --model est2genome \
  --minintron 30 --maxintron 200000 est.fa genome.fa

# emit GFF3 match / match_part records from the traceback
cargo run --release -p exonerate -- --model est2genome \
  --showvulgar no --showgff yes est.fa genome.fa
```

## Heuristic DNA search

`affine:local` uses exact k-mer seeds to select a padded candidate rectangle, then runs the unchanged exact Viterbi/traceback kernel inside that rectangle. Use `-E` or `--exhaustive` to force the full matrix. `--wordlen`, `--seedpadding`, and `--seedrepeat` control seeding; ambiguous or seedless inputs fall back to exhaustive alignment.

On the checked 3,000 x 5,000 nt planted-alignment fixture, the release build produced the identical `M 800 800` vulgar path in 0.02 s and 28,116 KiB maximum RSS, versus 0.20 s and 221,816 KiB for `-E`. These figures are regression evidence, not a cross-machine performance guarantee.

See [`exonerate.md`](exonerate.md) for the staged architecture and
[`COMPATIBILITY.md`](COMPATIBILITY.md) for the implemented compatibility surface.
