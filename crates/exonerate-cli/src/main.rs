use exonerate_core::{
    Alignment, HeuristicConfig, IntronScoring, Model, RawStep, Scoring, Sequence, Strand,
    align_cdna_to_genome_database_heuristic, align_cdna_to_genome_database_suboptimal,
    align_cdna_to_genome_database_with_dp_memory_stranded,
    align_coding_to_genome_database_heuristic, align_coding_to_genome_database_suboptimal,
    align_coding_to_genome_database_with_dp_memory, align_coding2coding_database,
    align_coding2coding_database_suboptimal, align_database_heuristic, align_database_suboptimal,
    align_database_with_dp_memory, align_est2genome_database_heuristic,
    align_est2genome_database_suboptimal, align_est2genome_database_with_dp_memory,
    align_genome_to_genome_database_heuristic_stranded,
    align_genome_to_genome_database_heuristic_stranded_with_rolling_scores,
    align_genome_to_genome_database_stranded, align_genome_to_genome_database_suboptimal_stranded,
    align_ner_database, align_ner_database_suboptimal, align_protein_database_suboptimal,
    align_protein_database_with_dp_memory, align_protein_to_dna_database,
    align_protein_to_dna_database_suboptimal, align_protein_to_genome_bestfit_database_heuristic,
    align_protein_to_genome_database_heuristic, align_protein_to_genome_database_suboptimal,
    align_protein_to_genome_database_with_dp_memory, align_protein_ungapped_database_suboptimal,
    align_ungapped_database_suboptimal, align_ungapped_translated_database,
    align_ungapped_translated_database_suboptimal, dna_self_score, dna_substitution_score,
    protein_self_score, protein_substitution_score, read_fasta, reverse_complement, translate_dna,
    translated_self_score,
};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::process::ExitCode;

fn usage() -> &'static str {
    "Usage: exonerate-rs [--shorthelp|--help] [--version] [--model MODEL] [--querytype dna|protein] [--targettype dna|protein] [--querychunkid N --querychunktotal N] [--targetchunkid N --targetchunktotal N] [--gapopen N] [--gapextend N] [--codongapopen N] [--codongapextend N] [--frameshift N] [--minintron N] [--maxintron N] [--intronpenalty N] [--forcegtag yes|no] [--minner N] [--maxner N] [--neropen N] [--wordlen N] [--seedpadding N] [--seedrepeat N] [-D N|--dpmemory N] [--score N] [--percent N] [--bestn N] [--ryo FORMAT] [-q QUERY.fa] [-t TARGET.fa] [--subopt yes|no] [--exhaustive [yes|no]] [--revcomp yes|no] [--forwardcoordinates yes|no] [--forwardonly] [--showsugar yes|no] [--showcigar yes|no] [--showvulgar yes|no] [--showgff yes|no] [--showquerygff yes|no] [--showtargetgff yes|no] QUERY.fa TARGET.fa\n\nImplemented models: ungapped, ungapped:trans, affine:global, affine:bestfit, affine:local, affine:overlap, coding2coding, coding2genome, cdna2genome, protein2dna, protein2dna:bestfit, protein2genome, protein2genome:bestfit, est2genome, genome2genome, ner"
}

fn version() -> String {
    format!("exonerate-rs {}", env!("CARGO_PKG_VERSION"))
}
fn yes_no(value: &str) -> Result<bool, String> {
    match value {
        "yes" | "y" | "true" => Ok(true),
        "no" | "n" | "false" => Ok(false),
        _ => Err(format!("expected yes or no, got {value:?}")),
    }
}
fn strand_symbol(strand: Strand) -> char {
    match strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        Strand::Unknown => '.',
    }
}

fn sequence_type(sequence: &Sequence) -> &'static str {
    if sequence.bases.iter().all(|base| {
        matches!(
            base.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'U'
                | b'N'
                | b'R'
                | b'Y'
                | b'M'
                | b'K'
                | b'S'
                | b'W'
                | b'H'
                | b'B'
                | b'V'
                | b'D'
        )
    }) {
        "dna"
    } else {
        "protein"
    }
}

/// Upstream's `Sequence_print_fasta_block`: sequence RYO fields have no
/// header, but they are still FASTA-wrapped at 70 columns and end in `\n`.
fn fasta_block_text(bases: &[u8]) -> String {
    if bases.is_empty() {
        return String::new();
    }
    let mut output = String::with_capacity(bases.len() + bases.len() / 70 + 1);
    for chunk in bases.chunks(70) {
        output.push_str(&String::from_utf8_lossy(chunk));
        output.push('\n');
    }
    output
}

fn aligned_sequence_text(sequence: &Sequence, start: u64, end: u64, strand: Strand) -> String {
    let begin = start.min(end) as usize;
    let finish = start.max(end) as usize;
    let mut bases =
        sequence.bases[begin.min(sequence.bases.len())..finish.min(sequence.bases.len())].to_vec();
    if strand == Strand::Reverse {
        bases = reverse_complement(&bases);
    }
    String::from_utf8_lossy(&bases).into_owned()
}

#[derive(Default)]
struct RyoStats {
    total: u64,
    identical: u64,
    similar: u64,
    gaps: u64,
    match_score: i64,
    self_score: i64,
}

struct CodingRegion {
    begin: u64,
    end: u64,
    bases: Vec<u8>,
}

fn ryo_symbol(sequence: &Sequence, position: u64, advance: u32, strand: Strand) -> Option<u8> {
    let advance = advance as usize;
    if !matches!(advance, 1 | 3) {
        return None;
    }
    let start = if strand == Strand::Reverse {
        position.checked_sub(advance as u64)? as usize
    } else {
        position as usize
    };
    let end = start.checked_add(advance)?;
    let bases = sequence.bases.get(start..end)?;
    let oriented = if strand == Strand::Reverse {
        reverse_complement(bases)
    } else {
        bases.to_vec()
    };
    Some(if advance == 3 {
        translate_dna(&oriented, 0)[0]
    } else {
        oriented[0]
    })
}

fn ryo_stats(alignment: &Alignment, query: &Sequence, target: &Sequence) -> RyoStats {
    let mut stats = RyoStats::default();
    let mut query_position = if alignment.query_strand == Strand::Reverse {
        alignment.query_start.max(alignment.query_end)
    } else {
        alignment.query_start.min(alignment.query_end)
    };
    let mut target_position = if alignment.target_strand == Strand::Reverse {
        alignment.target_start.max(alignment.target_end)
    } else {
        alignment.target_start.min(alignment.target_end)
    };
    for run in &alignment.trace {
        for _ in 0..run.repeats {
            if run.op == exonerate_core::Op::Match {
                let query_symbol = ryo_symbol(
                    query,
                    query_position,
                    run.query_advance,
                    alignment.query_strand,
                );
                let target_symbol = ryo_symbol(
                    target,
                    target_position,
                    run.target_advance,
                    alignment.target_strand,
                );
                if let (Some(query_symbol), Some(target_symbol)) = (query_symbol, target_symbol) {
                    let protein = run.query_advance == 3
                        || run.target_advance == 3
                        || sequence_type(query) == "protein"
                        || sequence_type(target) == "protein";
                    let score = if protein {
                        protein_substitution_score(query_symbol, target_symbol)
                    } else {
                        dna_substitution_score(query_symbol, target_symbol)
                    };
                    let self_score = if protein {
                        protein_substitution_score(query_symbol, query_symbol)
                    } else {
                        dna_substitution_score(query_symbol, query_symbol)
                    };
                    stats.total += 1;
                    stats.identical += u64::from(query_symbol.eq_ignore_ascii_case(&target_symbol));
                    stats.similar += u64::from(score > 0);
                    stats.match_score += i64::from(score);
                    stats.self_score += i64::from(self_score);
                }
            } else if matches!(
                run.op,
                exonerate_core::Op::Insert | exonerate_core::Op::Delete
            ) {
                stats.gaps += 1;
            }
            let query_advance = u64::from(run.query_advance);
            let target_advance = u64::from(run.target_advance);
            query_position = if alignment.query_strand == Strand::Reverse {
                query_position.saturating_sub(query_advance)
            } else {
                query_position + query_advance
            };
            target_position = if alignment.target_strand == Strand::Reverse {
                target_position.saturating_sub(target_advance)
            } else {
                target_position + target_advance
            };
        }
    }
    stats
}

fn ryo_oriented_bases(
    sequence: &Sequence,
    position: u64,
    advance: u32,
    strand: Strand,
) -> Option<Vec<u8>> {
    let advance = advance as usize;
    let start = if strand == Strand::Reverse {
        position.checked_sub(advance as u64)? as usize
    } else {
        position as usize
    };
    let end = start.checked_add(advance)?;
    let bases = sequence.bases.get(start..end)?;
    Some(if strand == Strand::Reverse {
        reverse_complement(bases)
    } else {
        bases.to_vec()
    })
}

fn ryo_coding_region(
    alignment: &Alignment,
    sequence: &Sequence,
    on_query: bool,
) -> Option<CodingRegion> {
    let strand = if on_query {
        alignment.query_strand
    } else {
        alignment.target_strand
    };
    let (start, end) = if on_query {
        (alignment.query_start, alignment.query_end)
    } else {
        (alignment.target_start, alignment.target_end)
    };
    let mut position = if strand == Strand::Reverse {
        start.max(end)
    } else {
        start.min(end)
    };
    let mut region: Option<CodingRegion> = None;
    for run in &alignment.trace {
        let advance = if on_query {
            run.query_advance
        } else {
            run.target_advance
        };
        for _ in 0..run.repeats {
            let include = match run.op {
                exonerate_core::Op::Match => advance == 3,
                exonerate_core::Op::SplitCodon => advance > 0,
                exonerate_core::Op::Insert | exonerate_core::Op::Delete => advance == 3,
                _ => false,
            };
            if include {
                let bases = ryo_oriented_bases(sequence, position, advance, strand)?;
                if let Some(region) = &mut region {
                    region.bases.extend(bases);
                    if run.op == exonerate_core::Op::Match && advance == 3 {
                        region.end = position;
                    }
                } else if run.op == exonerate_core::Op::Match && advance == 3 {
                    region = Some(CodingRegion {
                        begin: position,
                        end: position,
                        bases,
                    });
                }
            }
            position = if strand == Strand::Reverse {
                position.saturating_sub(u64::from(advance))
            } else {
                position + u64::from(advance)
            };
        }
    }
    region
}

fn percent(numerator: i64, denominator: i64) -> String {
    if denominator == 0 {
        "0.00".to_owned()
    } else {
        format!("{:.2}", numerator as f64 * 100.0 / denominator as f64)
    }
}

fn report_tail(report: &str) -> String {
    report
        .split_whitespace()
        .skip(10)
        .collect::<Vec<_>>()
        .join(" ")
}

/// Return the complete definition text for each FASTA record, excluding `>`.
///
/// The core reader intentionally retains only the identifier because that is
/// all alignment needs.  RYO's `%qd` and `%td` are reporting fields, however,
/// and upstream prints the full definition line there.
fn read_fasta_definitions(path: &str) -> Result<HashMap<String, String>, String> {
    let file = File::open(path).map_err(|error| format!("failed to open {path}: {error}"))?;
    let mut definitions = HashMap::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|error| format!("failed to read {path}: {error}"))?;
        let Some(definition) = line.strip_prefix('>') else {
            continue;
        };
        let definition = definition.trim_end_matches('\r');
        let id = definition
            .split_whitespace()
            .next()
            .filter(|id| !id.is_empty())
            .ok_or_else(|| format!("empty FASTA definition in {path}"))?;
        definitions.insert(id.to_owned(), definition.to_owned());
    }
    Ok(definitions)
}

/// Select the FASTA records in the same byte-range chunk used by upstream's
/// `FastaDB_open_list_with_limit`.  Boundaries are advanced to the next FASTA
/// header so a record is never split between workers.
fn fasta_chunk(
    path: &str,
    records: Vec<Sequence>,
    chunk_id: usize,
    chunk_total: usize,
) -> Result<Vec<Sequence>, String> {
    if chunk_total == 0 {
        return Ok(records);
    }
    if chunk_id == 0 || chunk_id > chunk_total {
        return Err(format!("chunk id should be between 1 and {chunk_total}"));
    }
    let bytes = fs::read(path).map_err(|error| format!("failed to read {path}: {error}"))?;
    if bytes.is_empty() {
        return Ok(records);
    }
    let find_next_start = |position: usize| {
        let mut previous = b'\n';
        for (offset, byte) in bytes.iter().enumerate().skip(position) {
            if *byte == b'>' && previous == b'\n' {
                return offset;
            }
            previous = *byte;
        }
        bytes.len() - 1
    };
    let chunk_size = bytes.len() / chunk_total;
    let start = find_next_start((chunk_id - 1) * chunk_size);
    let stop = if chunk_id == chunk_total {
        bytes.len() - 1
    } else {
        find_next_start(chunk_id * chunk_size)
    };
    let headers = bytes
        .iter()
        .enumerate()
        .filter_map(|(offset, byte)| {
            (*byte == b'>' && (offset == 0 || bytes[offset - 1] == b'\n')).then_some(offset)
        })
        .collect::<Vec<_>>();
    if headers.len() != records.len() {
        return Err(format!(
            "failed to locate FASTA record boundaries in {path}"
        ));
    }
    Ok(records
        .into_iter()
        .zip(headers)
        .filter_map(|(record, header)| ((start..stop).contains(&header)).then_some(record))
        .collect())
}

#[allow(clippy::too_many_arguments)]
fn render_ryo_plain(
    format: &str,
    alignment: &Alignment,
    query: &Sequence,
    target: &Sequence,
    query_definition: &str,
    target_definition: &str,
    model: &str,
    rank: usize,
) -> Result<String, String> {
    let chars = format.chars().collect::<Vec<_>>();
    let stats = ryo_stats(alignment, query, target);
    let mut output = String::new();
    let mut index = 0;
    while index < chars.len() {
        if chars[index] == char::from(92) {
            index += 1;
            if index >= chars.len() {
                output.push(char::from(92));
                break;
            }
            output.push(match chars[index] {
                'n' => char::from(10),
                't' => char::from(9),
                other => other,
            });
            index += 1;
            continue;
        }
        if chars[index] != '%' {
            output.push(chars[index]);
            index += 1;
            continue;
        }
        index += 1;
        let Some(&token) = chars.get(index) else {
            return Err("trailing percent in ryo format".into());
        };
        index += 1;
        match token {
            '%' => output.push('%'),
            's' => output.push_str(&alignment.score.to_string()),
            'r' => output.push_str(&rank.to_string()),
            'm' => output.push_str(model),
            'e' => {
                let Some(&field) = chars.get(index) else {
                    return Err("incomplete equivalenced ryo token".into());
                };
                index += 1;
                match field {
                    't' => output.push_str(&stats.total.to_string()),
                    'i' => output.push_str(&stats.identical.to_string()),
                    's' => output.push_str(&stats.similar.to_string()),
                    'm' => output.push_str(&(stats.total - stats.identical).to_string()),
                    other => return Err(format!("unknown equivalenced ryo token %e{other}")),
                }
            }
            'p' => {
                let Some(&field) = chars.get(index) else {
                    return Err("incomplete percent ryo token".into());
                };
                index += 1;
                let total = stats.total as i64;
                match field {
                    'c' => output.push_str(&percent(total, query.bases.len() as i64)),
                    'i' => output.push_str(&percent(stats.identical as i64, total)),
                    'I' => {
                        output.push_str(&percent(stats.identical as i64, total + stats.gaps as i64))
                    }
                    's' => output.push_str(&percent(stats.similar as i64, total)),
                    'S' => output.push_str(&percent(stats.match_score, stats.self_score)),
                    other => return Err(format!("unknown percent ryo token %p{other}")),
                }
            }
            'g' => output.push(strand_symbol(alignment.target_strand)),
            'S' => output.push_str(alignment.sugar().strip_prefix("sugar: ").unwrap()),
            'C' => {
                output.push(' ');
                output.push_str(&report_tail(&alignment.cigar()));
            }
            'V' => output.push_str(&report_tail(&alignment.vulgar())),
            'q' | 't' => {
                let is_query = token == 'q';
                let sequence = if is_query { query } else { target };
                let definition = if is_query {
                    query_definition
                } else {
                    target_definition
                };
                let (start, end, strand) = if is_query {
                    (
                        alignment.query_start,
                        alignment.query_end,
                        alignment.query_strand,
                    )
                } else {
                    (
                        alignment.target_start,
                        alignment.target_end,
                        alignment.target_strand,
                    )
                };
                let Some(&field) = chars.get(index) else {
                    return Err("incomplete query or target ryo token".into());
                };
                index += 1;
                match field {
                    'i' => output.push_str(&sequence.id),
                    'd' => output.push_str(definition),
                    'l' => output.push_str(&sequence.bases.len().to_string()),
                    's' => {
                        let bases = if strand == Strand::Reverse {
                            reverse_complement(&sequence.bases)
                        } else {
                            sequence.bases.clone()
                        };
                        output.push_str(&fasta_block_text(&bases));
                    }
                    'S' => output.push(strand_symbol(strand)),
                    't' => output.push_str(sequence_type(sequence)),
                    'a' | 'c' => {
                        let Some(&region_field) = chars.get(index) else {
                            return Err("incomplete aligned-region ryo token".into());
                        };
                        index += 1;
                        let coding_region = (field == 'c')
                            .then(|| ryo_coding_region(alignment, sequence, is_query))
                            .flatten();
                        if field == 'c' && coding_region.is_none() {
                            return Err(format!(
                                "%{token}c requires a coding region in the alignment"
                            ));
                        }
                        let (region_start, region_end, region_sequence) =
                            if let Some(region) = coding_region {
                                (
                                    region.begin,
                                    region.end,
                                    String::from_utf8_lossy(&region.bases).into_owned(),
                                )
                            } else {
                                (
                                    start,
                                    end,
                                    aligned_sequence_text(sequence, start, end, strand),
                                )
                            };
                        match region_field {
                            'b' => output.push_str(&region_start.to_string()),
                            'e' => output.push_str(&region_end.to_string()),
                            'l' => output.push_str(&region_sequence.len().to_string()),
                            's' => output.push_str(&fasta_block_text(region_sequence.as_bytes())),
                            other => {
                                return Err(format!(
                                    "unknown ryo region token %{token}{field}{other}"
                                ));
                            }
                        }
                    }
                    other => return Err(format!("unknown ryo token %{token}{other}")),
                }
            }
            other => return Err(format!("unknown ryo token %{other}")),
        }
    }
    Ok(output)
}

fn transition_metadata(model: &str, transition_id: u16) -> (&'static str, &'static str) {
    if model == "cdna2genome" {
        return match transition_id {
            0 | 300 | 400 => ("start to match", "none"),
            1 | 301 | 401 => ("match", "match"),
            2 | 302 | 402 => ("match to insert", "gap"),
            3 | 303 | 403 => ("insert", "gap"),
            4 | 304 | 404 => ("insert to match", "none"),
            5 | 305 | 405 => ("match to delete", "gap"),
            6 | 306 | 406 => ("delete", "gap"),
            7 | 307 | 407 => ("delete to match", "none"),
            8 => ("frameshift open 1 query", "frameshift"),
            9 => ("frameshift open 2 query", "frameshift"),
            10 => ("frameshift close 0 query", "none"),
            11 => ("frameshift close 3 query", "frameshift"),
            12 => ("frameshift open 1 target", "frameshift"),
            13 => ("frameshift open 2 target", "frameshift"),
            14 => ("frameshift close 0 target", "none"),
            15 => ("frameshift close 3 target", "frameshift"),
            16 | 308 | 408 | 320 | 321 => ("match to end", "none"),
            310 | 410 => ("(START) to intron forward", "5'ss"),
            311 | 411 => ("target intron loop forward", "intron"),
            312 | 412 => ("intron forward to (END)", "3'ss"),
            101 => ("(START) to phase1pre phase target intron -T", "split codon"),
            102 => ("(START) to phase2pre phase target intron -T", "split codon"),
            110 => ("(START) to intron 0:0 phase target intron -T", "5'ss"),
            111 => ("(START) to intron 1:2 phase target intron -T", "5'ss"),
            112 => ("(START) to intron 2:1 phase target intron -T", "5'ss"),
            120 => ("target intron loop 0:0 phase target intron -T", "intron"),
            121 => ("target intron loop 1:2 phase target intron -T", "intron"),
            122 => ("target intron loop 2:1 phase target intron -T", "intron"),
            130 => ("intron 0:0 phase target intron -T to (END)", "3'ss"),
            131 => ("intron 1:2 phase target intron -T to (END)", "3'ss"),
            132 => ("intron 2:1 phase target intron -T to (END)", "3'ss"),
            140 => ("phase0post phase target intron -T to (END)", "match"),
            141 => ("phase1post phase target intron -T to (END)", "split codon"),
            142 => ("phase2post phase target intron -T to (END)", "split codon"),
            _ => ("unknown", "none"),
        };
    }
    if model == "genome2genome" {
        return match transition_id {
            29 => ("start to match", "none"),
            30 | 51 => ("match", "match"),
            31 | 52 => ("match to insert", "gap"),
            32 | 54 => ("insert", "gap"),
            33 | 53 => ("match to delete", "gap"),
            34 | 55 => ("delete", "gap"),
            35 | 59 => ("insert to match", "none"),
            36 | 79 => ("delete to match", "none"),
            39 | 89 | 90 => ("match to end", "none"),
            40 => ("(START) to intron forward", "5'ss"),
            41 => ("target intron loop forward", "intron"),
            42 => ("intron forward to (END)", "3'ss"),
            43 => ("(START) to intron query", "5'ss"),
            44 => ("query intron loop query", "intron"),
            45 => ("intron query to (END)", "3'ss"),
            46 => ("(START) to intron joint", "5'ss"),
            47 => ("query intron loop joint", "intron"),
            49 => ("target intron loop joint", "intron"),
            48 => ("intron joint to (END)", "3'ss"),
            50 | 58 => ("match to end", "none"),
            80 => ("frameshift open 1 query", "frameshift"),
            81 => ("frameshift open 2 query", "frameshift"),
            82 => ("frameshift close 0 query", "none"),
            83 => ("frameshift close 3 query", "frameshift"),
            84 => ("frameshift open 1 target", "frameshift"),
            85 => ("frameshift open 2 target", "frameshift"),
            86 => ("frameshift close 0 target", "none"),
            87 => ("frameshift close 3 target", "frameshift"),
            100 => ("(START) to phase0pre phase target intron -T", "match"),
            101 => ("(START) to phase1pre phase target intron -T", "split codon"),
            102 => ("(START) to phase2pre phase target intron -T", "split codon"),
            110 => ("(START) to intron 0:0 phase target intron -T", "5'ss"),
            111 => ("(START) to intron 1:2 phase target intron -T", "5'ss"),
            112 => ("(START) to intron 2:1 phase target intron -T", "5'ss"),
            120 => ("target intron loop 0:0 phase target intron -T", "intron"),
            121 => ("target intron loop 1:2 phase target intron -T", "intron"),
            122 => ("target intron loop 2:1 phase target intron -T", "intron"),
            130 => ("intron 0:0 phase target intron -T to (END)", "3'ss"),
            131 => ("intron 1:2 phase target intron -T to (END)", "3'ss"),
            132 => ("intron 2:1 phase target intron -T to (END)", "3'ss"),
            140 => ("phase0post phase target intron -T to (END)", "match"),
            141 => ("phase1post phase target intron -T to (END)", "split codon"),
            142 => ("phase2post phase target intron -T to (END)", "split codon"),
            150 => ("(START) to phase0pre phase query Q-", "match"),
            151 => ("(START) to phase1pre phase query Q-", "split codon"),
            152 => ("(START) to phase2pre phase query Q-", "split codon"),
            160 => ("(START) to intron 0:0 phase query Q-", "5'ss"),
            161 => ("(START) to intron 1:2 phase query Q-", "5'ss"),
            162 => ("(START) to intron 2:1 phase query Q-", "5'ss"),
            170 => ("query intron loop 0:0 phase query Q-", "intron"),
            171 => ("query intron loop 1:2 phase query Q-", "intron"),
            172 => ("query intron loop 2:1 phase query Q-", "intron"),
            180 => ("intron 0:0 phase query Q- to (END)", "3'ss"),
            181 => ("intron 1:2 phase query Q- to (END)", "3'ss"),
            182 => ("intron 2:1 phase query Q- to (END)", "3'ss"),
            190 => ("phase0post phase query Q- to (END)", "match"),
            191 => ("phase1post phase query Q- to (END)", "split codon"),
            192 => ("phase2post phase query Q- to (END)", "split codon"),
            200 => ("(START) to phase0pre phase joint QT", "match"),
            201 => ("(START) to phase1pre phase joint QT", "split codon"),
            202 => ("(START) to phase2pre phase joint QT", "split codon"),
            210 => ("(START) to intron 0:0 phase joint QT", "5'ss"),
            211 => ("(START) to intron 1:2 phase joint QT", "5'ss"),
            212 => ("(START) to intron 2:1 phase joint QT", "5'ss"),
            220 => ("query intron loop 0:0 phase joint QT", "intron"),
            221 => ("query intron loop 1:2 phase joint QT", "intron"),
            222 => ("query intron loop 2:1 phase joint QT", "intron"),
            230 => ("target intron loop 0:0 phase joint QT", "intron"),
            231 => ("target intron loop 1:2 phase joint QT", "intron"),
            232 => ("target intron loop 2:1 phase joint QT", "intron"),
            240 => ("intron 0:0 phase joint QT to (END)", "3'ss"),
            241 => ("intron 1:2 phase joint QT to (END)", "3'ss"),
            242 => ("intron 2:1 phase joint QT to (END)", "3'ss"),
            250 => ("phase0post phase joint QT to (END)", "match"),
            251 => ("phase1post phase joint QT to (END)", "split codon"),
            252 => ("phase2post phase joint QT to (END)", "split codon"),
            _ => ("unknown", "none"),
        };
    }
    if model == "coding2genome" {
        return match transition_id {
            0 => ("start to match", "none"),
            1 => ("match", "match"),
            2 => ("match to insert", "gap"),
            3 => ("insert", "gap"),
            4 => ("insert to match", "none"),
            5 => ("match to delete", "gap"),
            6 => ("delete", "gap"),
            7 => ("delete to match", "none"),
            8 => ("frameshift open 1 query", "frameshift"),
            9 => ("frameshift open 2 query", "frameshift"),
            10 => ("frameshift close 0 query", "none"),
            11 => ("frameshift close 3 query", "frameshift"),
            12 => ("frameshift open 1 target", "frameshift"),
            13 => ("frameshift open 2 target", "frameshift"),
            14 => ("frameshift close 0 target", "none"),
            15 => ("frameshift close 3 target", "frameshift"),
            16 => ("match to end", "none"),
            101 => ("(START) to phase1pre phase target intron -T", "split codon"),
            102 => ("(START) to phase2pre phase target intron -T", "split codon"),
            110 => ("(START) to intron 0:0 phase target intron -T", "5'ss"),
            111 => ("(START) to intron 1:2 phase target intron -T", "5'ss"),
            112 => ("(START) to intron 2:1 phase target intron -T", "5'ss"),
            120 => ("target intron loop 0:0 phase target intron -T", "intron"),
            121 => ("target intron loop 1:2 phase target intron -T", "intron"),
            122 => ("target intron loop 2:1 phase target intron -T", "intron"),
            130 => ("intron 0:0 phase target intron -T to (END)", "3'ss"),
            131 => ("intron 1:2 phase target intron -T to (END)", "3'ss"),
            132 => ("intron 2:1 phase target intron -T to (END)", "3'ss"),
            140 => ("phase0post phase target intron -T to (END)", "match"),
            141 => ("phase1post phase target intron -T to (END)", "split codon"),
            142 => ("phase2post phase target intron -T to (END)", "split codon"),
            _ => ("unknown", "none"),
        };
    }
    if model == "est2genome" {
        return match transition_id {
            0 => ("start to match forward", "none"),
            1 => ("match forward", "match"),
            2 => ("match to insert forward", "gap"),
            3 => ("match to delete forward", "gap"),
            4 => ("insert forward", "gap"),
            5 => ("delete forward", "gap"),
            6 => ("(START) to intron forward", "5'ss"),
            7 => ("target intron loop forward", "intron"),
            8 => ("intron forward to (END)", "3'ss"),
            9 => ("insert to match forward", "none"),
            10 => ("delete to match forward", "none"),
            11 => ("match to end forward", "none"),
            _ => ("unknown", "none"),
        };
    }
    if model.starts_with("protein2genome") {
        return match transition_id {
            0 => ("start to match", "none"),
            1 => ("match", "match"),
            2 => ("match to insert", "gap"),
            3 => ("match to delete", "gap"),
            4 => ("insert", "gap"),
            5 => ("delete", "gap"),
            6 => ("insert to match", "none"),
            7 => ("delete to match", "none"),
            8 => ("frameshift open 1 p2g", "frameshift"),
            9 => ("frameshift open 2 p2g", "frameshift"),
            10 => ("frameshift close 0 p2g", "none"),
            11 => ("frameshift close 3 p2g", "frameshift"),
            12 => ("match to end", "none"),
            101 => ("(START) to phase1pre phase-T", "split codon"),
            102 => ("(START) to phase2pre phase-T", "split codon"),
            110 => ("(START) to intron 0:0 phase-T", "5'ss"),
            111 => ("(START) to intron 1:2 phase-T", "5'ss"),
            112 => ("(START) to intron 2:1 phase-T", "5'ss"),
            120 => ("target intron loop 0:0 phase-T", "intron"),
            121 => ("target intron loop 1:2 phase-T", "intron"),
            122 => ("target intron loop 2:1 phase-T", "intron"),
            130 => ("intron 0:0 phase-T to (END)", "3'ss"),
            131 => ("intron 1:2 phase-T to (END)", "3'ss"),
            132 => ("intron 2:1 phase-T to (END)", "3'ss"),
            140 => ("phase0post phase-T to (END)", "match"),
            141 => ("phase1post phase-T to (END)", "split codon"),
            142 => ("phase2post phase-T to (END)", "split codon"),
            _ => ("unknown", "none"),
        };
    }
    if model.starts_with("coding2coding") {
        return match transition_id {
            0 => ("start to match", "none"),
            1 => ("match", "match"),
            2 => ("match to insert", "gap"),
            3 => ("insert", "gap"),
            4 => ("insert to match", "none"),
            5 => ("match to delete", "gap"),
            6 => ("delete", "gap"),
            7 => ("delete to match", "none"),
            8 => ("frameshift open 1 query", "frameshift"),
            9 => ("frameshift open 2 query", "frameshift"),
            10 => ("frameshift close 0 query", "none"),
            11 => ("frameshift close 3 query", "frameshift"),
            12 => ("frameshift open 1 target", "frameshift"),
            13 => ("frameshift open 2 target", "frameshift"),
            14 => ("frameshift close 0 target", "none"),
            15 => ("frameshift close 3 target", "frameshift"),
            16 => ("match to end", "none"),
            _ => ("unknown", "none"),
        };
    }
    if model.starts_with("protein2dna") {
        return match transition_id {
            0 => ("start to match", "none"),
            1 => ("match", "match"),
            2 => ("match to insert", "gap"),
            3 => ("match to delete", "gap"),
            4 => ("insert", "gap"),
            5 => ("delete", "gap"),
            6 => ("insert to match", "none"),
            7 => ("delete to match", "none"),
            8 => ("frameshift open 1 p2d", "frameshift"),
            9 => ("frameshift open 2 p2d", "frameshift"),
            10 => ("frameshift close 0 p2d", "none"),
            11 => ("frameshift close 3 p2d", "frameshift"),
            12 => ("match to end", "none"),
            _ => ("unknown", "none"),
        };
    }
    if model.starts_with("ungapped") {
        return match transition_id {
            0 => ("start to match", "none"),
            1 => ("match", "match"),
            2 => ("match to end", "none"),
            _ => ("unknown", "none"),
        };
    }
    match transition_id {
        0 => ("start to match", "none"),
        1 => ("match", "match"),
        2 => ("match to insert", "gap"),
        3 => ("match to delete", "gap"),
        4 => ("insert", "gap"),
        5 => ("delete", "gap"),
        6 => ("insert to match", "none"),
        7 => ("delete to match", "none"),
        8 => ("match to end", "none"),
        _ => ("unknown", "none"),
    }
}

#[allow(clippy::too_many_arguments)]
fn render_transition_ryo(
    format: &str,
    step: &RawStep,
    model: &str,
    query: &Sequence,
    target: &Sequence,
    query_begin: u64,
    target_begin: u64,
    query_strand: Strand,
    target_strand: Strand,
) -> Result<String, String> {
    let query_advance = u64::from(step.query_advance);
    let target_advance = u64::from(step.target_advance);
    let (transition_name, transition_label) = transition_metadata(model, step.transition_id);
    let query_end = if query_strand == Strand::Reverse {
        query_begin.saturating_sub(query_advance)
    } else {
        query_begin + query_advance
    };
    let target_end = if target_strand == Strand::Reverse {
        target_begin.saturating_sub(target_advance)
    } else {
        target_begin + target_advance
    };
    let chars = format.chars().collect::<Vec<_>>();
    let mut output = String::new();
    let mut index = 0;
    while index < chars.len() {
        if chars[index] == char::from(92) {
            index += 1;
            if index >= chars.len() {
                output.push(char::from(92));
                break;
            }
            output.push(match chars[index] {
                'n' => char::from(10),
                't' => char::from(9),
                other => other,
            });
            index += 1;
            continue;
        }
        if chars[index] != '%' {
            output.push(chars[index]);
            index += 1;
            continue;
        }
        if chars.get(index + 1) != Some(&'P') {
            return Err("only %P tokens are valid inside a transition block".into());
        }
        index += 2;
        let Some(&kind) = chars.get(index) else {
            return Err("incomplete %P token".into());
        };
        index += 1;
        match kind {
            'n' => output.push_str(transition_name),
            'l' => output.push_str(transition_label),
            's' => output.push_str(&step.score.to_string()),
            'q' | 't' => {
                let is_query = kind == 'q';
                let Some(&field) = chars.get(index) else {
                    return Err("incomplete %P query or target token".into());
                };
                index += 1;
                let (sequence, begin, end, advance, strand) = if is_query {
                    (query, query_begin, query_end, query_advance, query_strand)
                } else {
                    (
                        target,
                        target_begin,
                        target_end,
                        target_advance,
                        target_strand,
                    )
                };
                match field {
                    's' => output.push_str(&aligned_sequence_text(sequence, begin, end, strand)),
                    'a' => output.push_str(&advance.to_string()),
                    'b' => output.push_str(&begin.to_string()),
                    'e' => output.push_str(&end.to_string()),
                    other => return Err(format!("unknown transition token %P{kind}{other}")),
                }
            }
            other => return Err(format!("unknown transition token %P{other}")),
        }
    }
    Ok(output)
}

#[allow(clippy::too_many_arguments)]
fn render_ryo(
    format: &str,
    alignment: &Alignment,
    query: &Sequence,
    target: &Sequence,
    query_definition: &str,
    target_definition: &str,
    model: &str,
    rank: usize,
) -> Result<String, String> {
    let chars = format.chars().collect::<Vec<_>>();
    let Some(open) = chars.iter().position(|character| *character == '{') else {
        return render_ryo_plain(
            format,
            alignment,
            query,
            target,
            query_definition,
            target_definition,
            model,
            rank,
        );
    };
    let close = chars
        .iter()
        .enumerate()
        .skip(open + 1)
        .find_map(|(index, character)| (*character == '}').then_some(index))
        .ok_or("unclosed transition block in ryo format")?;
    let prefix = chars[..open].iter().collect::<String>();
    let body = chars[open + 1..close].iter().collect::<String>();
    let suffix = chars[close + 1..].iter().collect::<String>();
    let mut output = render_ryo_plain(
        &prefix,
        alignment,
        query,
        target,
        query_definition,
        target_definition,
        model,
        rank,
    )?;
    let mut query_position = alignment.query_start;
    let mut target_position = alignment.target_start;
    let raw_steps = if alignment.raw_trace.is_empty() {
        alignment
            .trace
            .iter()
            .flat_map(|run| {
                (0..run.repeats).map(move |_| RawStep {
                    transition_id: run.transition_id,
                    query_advance: run.query_advance,
                    target_advance: run.target_advance,
                    score: 0,
                })
            })
            .collect::<Vec<_>>()
    } else {
        alignment.raw_trace.clone()
    };
    for step in &raw_steps {
        output.push_str(&render_transition_ryo(
            &body,
            step,
            model,
            query,
            target,
            query_position,
            target_position,
            alignment.query_strand,
            alignment.target_strand,
        )?);
        let query_advance = u64::from(step.query_advance);
        let target_advance = u64::from(step.target_advance);
        query_position = if alignment.query_strand == Strand::Reverse {
            query_position.saturating_sub(query_advance)
        } else {
            query_position + query_advance
        };
        target_position = if alignment.target_strand == Strand::Reverse {
            target_position.saturating_sub(target_advance)
        } else {
            target_position + target_advance
        };
    }
    output.push_str(&render_ryo(
        &suffix,
        alignment,
        query,
        target,
        query_definition,
        target_definition,
        model,
        rank,
    )?);
    Ok(output)
}

fn write_stdout(text: &str, newline: bool) -> Result<(), String> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let result = if newline {
        writeln!(handle, "{text}")
    } else {
        write!(handle, "{text}")
    };
    match result {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::BrokenPipe => {
            std::process::exit(0);
        }
        Err(error) => Err(format!("failed writing output: {error}")),
    }
}

fn pretty_alignment(
    alignment: &Alignment,
    query: &Sequence,
    target: &Sequence,
    report_model: &str,
) -> String {
    let rendered = alignment.pretty(query, target);
    // Upstream's generic affine views omit a model line, while the composed
    // and translated model views identify their concrete C4 graph there.
    let needs_model_line = !matches!(
        report_model,
        "ungapped:dna2dna"
            | "affine:global:dna2dna"
            | "affine:bestfit:dna2dna"
            | "affine:local:dna2dna"
            | "affine:overlap:dna2dna"
    );
    if !needs_model_line {
        return rendered;
    }
    let target_header = format!("        Target: {}\n", target.id);
    rendered.replacen(
        &target_header,
        &format!("{target_header}         Model: {report_model}\n"),
        1,
    )
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("exonerate-rs: {e}\n{}", usage());
            ExitCode::from(2)
        }
    }
}
fn run() -> Result<(), String> {
    let mut args = std::env::args().skip(1).peekable();
    let mut model = Model::Ungapped;
    let mut model_name = "ungapped".to_owned();
    let mut model_selected = false;
    let mut scoring = Scoring::default();
    let mut intron = IntronScoring::default();
    let mut est2genome = false;
    let mut coding2coding = false;
    let mut coding2genome = false;
    let mut cdna2genome = false;
    let mut protein2genome = false;
    let mut protein2genome_bestfit = false;
    let mut genome2genome = false;
    let mut ner = false;
    let mut ungapped_translated = false;
    let (mut min_ner, mut max_ner, mut ner_open) = (10_u32, 50_000_u32, -20_i32);
    let mut both = true;
    let mut forward_coordinates = true;
    let mut query_type = "dna".to_owned();
    let mut target_type = "dna".to_owned();
    let mut query_type_explicit = false;
    let mut target_type_explicit = false;
    let (mut query_chunk_id, mut query_chunk_total) = (0_usize, 0_usize);
    let (mut target_chunk_id, mut target_chunk_total) = (0_usize, 0_usize);
    let (mut sugar, mut cigar, mut vulgar, mut gff) = (false, false, true, false);
    let mut query_gff = false;
    let mut min_score: Option<i32> = Some(100);
    let mut percent: Option<f64> = None;
    let mut best_n: Option<usize> = None;
    let mut ryo: Option<String> = None;
    let (mut query_file, mut target_file): (Option<String>, Option<String>) = (None, None);
    let mut show_alignment = false;
    let mut subopt = true;
    let mut exhaustive = false;
    let mut dp_memory_mb = 32_usize;
    let mut heuristic = HeuristicConfig::default();
    let mut files = Vec::new();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-m" | "--model" => {
                if model_selected {
                    return Err("model was already specified".into());
                }
                let selected = args.next().ok_or("missing model")?;
                model_selected = true;
                model_name = selected.clone();
                match selected.as_str() {
                    "ungapped:trans" | "u:t" => {
                        ungapped_translated = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "protein2dna" | "p2d" => {
                        model = Model::Local;
                        query_type = "protein".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "protein2dna:bestfit" | "p2d:b" => {
                        model = Model::BestFit;
                        query_type = "protein".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "protein2genome:bestfit" | "p2g:b" => {
                        protein2genome_bestfit = true;
                        query_type = "protein".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "protein2genome" | "p2g" => {
                        protein2genome = true;
                        query_type = "protein".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "genome2genome" | "g2g" => {
                        genome2genome = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "cdna2genome" | "cd2g" => {
                        cdna2genome = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "coding2genome" | "c2g" => {
                        coding2genome = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "coding2coding" | "c2c" => {
                        coding2coding = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "ner" => {
                        ner = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    "est2genome" | "e2g" => {
                        est2genome = true;
                        query_type = "dna".to_owned();
                        target_type = "dna".to_owned();
                    }
                    _ => {
                        model = selected
                            .parse()
                            .map_err(|e: exonerate_core::Error| e.to_string())?;
                    }
                }
            }
            "-o" | "--gapopen" => {
                scoring.gap_open = args
                    .next()
                    .ok_or("missing gap-open penalty")?
                    .parse()
                    .map_err(|_| "invalid gap-open penalty")?
            }
            "-e" | "--gapextend" => {
                scoring.gap_extend = args
                    .next()
                    .ok_or("missing gap-extension penalty")?
                    .parse()
                    .map_err(|_| "invalid gap-extension penalty")?
            }
            "--codongapopen" => {
                scoring.codon_gap_open = args
                    .next()
                    .ok_or("missing codon gap-open penalty")?
                    .parse()
                    .map_err(|_| "invalid codon gap-open penalty")?
            }
            "--codongapextend" => {
                scoring.codon_gap_extend = args
                    .next()
                    .ok_or("missing codon gap-extend penalty")?
                    .parse()
                    .map_err(|_| "invalid codon gap-extend penalty")?
            }
            "-f" | "--frameshift" => {
                scoring.frameshift = args
                    .next()
                    .ok_or("missing frameshift penalty")?
                    .parse()
                    .map_err(|_| "invalid frameshift penalty")?
            }
            "--minintron" => {
                intron.min_len = args
                    .next()
                    .ok_or("missing minimum intron length")?
                    .parse()
                    .map_err(|_| "invalid minimum intron length")?
            }
            "--maxintron" => {
                intron.max_len = args
                    .next()
                    .ok_or("missing maximum intron length")?
                    .parse()
                    .map_err(|_| "invalid maximum intron length")?
            }
            "-i" | "--intronpenalty" => {
                intron.open_penalty = args
                    .next()
                    .ok_or("missing intron penalty")?
                    .parse()
                    .map_err(|_| "invalid intron penalty")?
            }
            "--forcegtag" => {
                intron.force_gtag = yes_no(&args.next().ok_or("missing forcegtag value")?)?
            }
            "--minner" => {
                min_ner = args
                    .next()
                    .ok_or("missing minimum NER length")?
                    .parse()
                    .map_err(|_| "invalid minimum NER length")?
            }
            "--maxner" => {
                max_ner = args
                    .next()
                    .ok_or("missing maximum NER length")?
                    .parse()
                    .map_err(|_| "invalid maximum NER length")?
            }
            "--neropen" => {
                ner_open = args
                    .next()
                    .ok_or("missing NER open penalty")?
                    .parse()
                    .map_err(|_| "invalid NER open penalty")?
            }
            "--wordlen" => {
                heuristic.word_len = args
                    .next()
                    .ok_or("missing word length")?
                    .parse()
                    .map_err(|_| "invalid word length")?
            }
            "--seedpadding" => {
                heuristic.padding = args
                    .next()
                    .ok_or("missing seed padding")?
                    .parse()
                    .map_err(|_| "invalid seed padding")?
            }
            "--seedrepeat" => {
                heuristic.max_word_occurrences = args
                    .next()
                    .ok_or("missing seed repeat cap")?
                    .parse()
                    .map_err(|_| "invalid seed repeat cap")?
            }
            "-q" | "--query" => query_file = Some(args.next().ok_or("missing query path")?),
            "-t" | "--target" => target_file = Some(args.next().ok_or("missing target path")?),
            "-Q" | "--querytype" => {
                query_type = args.next().ok_or("missing query type")?;
                query_type_explicit = true;
            }
            "-T" | "--targettype" => {
                target_type = args.next().ok_or("missing target type")?;
                target_type_explicit = true;
            }
            "--querychunkid" => {
                query_chunk_id = args
                    .next()
                    .ok_or("missing query chunk id")?
                    .parse()
                    .map_err(|_| "invalid query chunk id")?
            }
            "--querychunktotal" => {
                query_chunk_total = args
                    .next()
                    .ok_or("missing query chunk total")?
                    .parse()
                    .map_err(|_| "invalid query chunk total")?
            }
            "--targetchunkid" => {
                target_chunk_id = args
                    .next()
                    .ok_or("missing target chunk id")?
                    .parse()
                    .map_err(|_| "invalid target chunk id")?
            }
            "--targetchunktotal" => {
                target_chunk_total = args
                    .next()
                    .ok_or("missing target chunk total")?
                    .parse()
                    .map_err(|_| "invalid target chunk total")?
            }
            "-E" => exhaustive = true,
            "--exhaustive" => {
                exhaustive = if args
                    .peek()
                    .is_some_and(|value| matches!(value.as_str(), "yes" | "no"))
                {
                    yes_no(&args.next().expect("checked exhaustive value"))?
                } else {
                    true
                }
            }
            "-D" | "--dpmemory" => {
                dp_memory_mb = args
                    .next()
                    .ok_or("missing DP memory limit")?
                    .parse()
                    .map_err(|_| "invalid DP memory limit")?
            }
            "--score" => {
                min_score = Some(
                    args.next()
                        .ok_or("missing score threshold")?
                        .parse()
                        .map_err(|_| "invalid score threshold")?,
                )
            }
            "--percent" => {
                percent = Some(
                    args.next()
                        .ok_or("missing percent threshold")?
                        .parse()
                        .map_err(|_| "invalid percent threshold")?,
                )
            }
            "--ryo" => ryo = Some(args.next().ok_or("missing ryo format")?),
            "-n" | "--bestn" => {
                best_n = Some(
                    args.next()
                        .ok_or("missing bestn value")?
                        .parse()
                        .map_err(|_| "invalid bestn value")?,
                )
            }
            "-S" | "--subopt" => subopt = yes_no(&args.next().ok_or("missing subopt value")?)?,
            "--showalignment" => {
                show_alignment = yes_no(&args.next().ok_or("missing showalignment value")?)?
            }
            "-r" | "--revcomp" => both = yes_no(&args.next().ok_or("missing revcomp value")?)?,
            "--forwardcoordinates" => {
                forward_coordinates =
                    yes_no(&args.next().ok_or("missing forwardcoordinates value")?)?
            }
            "--forwardonly" => both = false,
            "--showsugar" => sugar = yes_no(&args.next().ok_or("missing showsugar value")?)?,
            "--showcigar" => cigar = yes_no(&args.next().ok_or("missing showcigar value")?)?,
            "--showvulgar" => vulgar = yes_no(&args.next().ok_or("missing showvulgar value")?)?,
            "--showgff" | "--showtargetgff" => {
                gff = yes_no(&args.next().ok_or("missing target GFF value")?)?
            }
            "--showquerygff" => query_gff = yes_no(&args.next().ok_or("missing query GFF value")?)?,
            "-h" | "--shorthelp" | "--help" => {
                write_stdout(usage(), true)?;
                return Ok(());
            }
            "-v" | "--version" => {
                write_stdout(&version(), true)?;
                return Ok(());
            }
            x if x.starts_with('-') => return Err(format!("unknown option {x}")),
            x => files.push(x.to_owned()),
        }
    }
    if query_file.is_some() || target_file.is_some() {
        if !files.is_empty() || query_file.is_none() || target_file.is_none() {
            return Err(
                "use either two positional FASTA files or both --query and --target".into(),
            );
        }
        files.push(query_file.expect("checked query path"));
        files.push(target_file.expect("checked target path"));
    }
    if files.len() != 2 {
        return Err("exactly QUERY.fa and TARGET.fa are required".into());
    }
    if intron.min_len > intron.max_len {
        return Err("minimum intron length must not exceed maximum intron length".into());
    }
    if heuristic.word_len == 0 || heuristic.word_len > 31 {
        return Err("word length must be between 1 and 31".into());
    }
    if heuristic.max_word_occurrences == 0 {
        return Err("seed repeat cap must be positive".into());
    }
    if min_ner > max_ner {
        return Err("minimum NER length must not exceed maximum NER length".into());
    }
    if !matches!(query_type.as_str(), "dna" | "protein") {
        return Err(format!(
            "query type must be dna or protein, got {query_type:?}"
        ));
    }
    if !matches!(target_type.as_str(), "dna" | "protein") {
        return Err(format!(
            "target type must be dna or protein, got {target_type:?}"
        ));
    }
    let requires_dna_pair = ungapped_translated
        || ner
        || est2genome
        || coding2coding
        || coding2genome
        || cdna2genome
        || genome2genome;
    let requires_protein_to_dna = protein2genome
        || protein2genome_bestfit
        || matches!(
            model_name.as_str(),
            "protein2dna" | "p2d" | "protein2dna:bestfit" | "p2d:b"
        );
    let q = fasta_chunk(
        &files[0],
        read_fasta(&files[0]).map_err(|e| e.to_string())?,
        query_chunk_id,
        query_chunk_total,
    )?;
    let t = fasta_chunk(
        &files[1],
        read_fasta(&files[1]).map_err(|e| e.to_string())?,
        target_chunk_id,
        target_chunk_total,
    )?;
    if !query_type_explicit && !requires_dna_pair && !requires_protein_to_dna {
        query_type = if q
            .iter()
            .any(|sequence| sequence_type(sequence) == "protein")
        {
            "protein".to_owned()
        } else {
            "dna".to_owned()
        };
    }
    if !target_type_explicit && !requires_dna_pair && !requires_protein_to_dna {
        target_type = if t
            .iter()
            .any(|sequence| sequence_type(sequence) == "protein")
        {
            "protein".to_owned()
        } else {
            "dna".to_owned()
        };
    }
    if requires_dna_pair && (query_type != "dna" || target_type != "dna") {
        return Err(format!(
            "model {model_name} requires DNA query and target sequences"
        ));
    }
    if requires_protein_to_dna && (query_type != "protein" || target_type != "dna") {
        return Err(format!(
            "model {model_name} requires a protein query and DNA target"
        ));
    }
    let query_definitions = if ryo.is_some() {
        read_fasta_definitions(&files[0])?
    } else {
        HashMap::new()
    };
    let target_definitions = if ryo.is_some() {
        read_fasta_definitions(&files[1])?
    } else {
        HashMap::new()
    };
    // Upstream's `--bestn` requests an HSP set, which necessarily enumerates
    // suboptimal paths before it applies the per-query rank cutoff.  Keep the
    // pre-existing single-path fallback for the model families whose
    // suboptimal executor has not been implemented yet.
    let suboptimal_supported = ungapped_translated
        || ner
        || est2genome
        || coding2coding
        || genome2genome
        || cdna2genome
        || coding2genome
        || protein2genome
        || protein2genome_bestfit
        || (query_type == "dna"
            && target_type == "dna"
            && matches!(model, Model::Ungapped | Model::Local))
        || (query_type == "protein"
            && target_type == "protein"
            && matches!(model, Model::Ungapped | Model::Local))
        || (query_type == "protein"
            && target_type == "dna"
            && matches!(model, Model::Ungapped | Model::Local | Model::BestFit));
    let mut alignments = if suboptimal_supported && (subopt || best_n.is_some()) {
        let plain_dna = !ungapped_translated
            && !ner
            && !coding2coding
            && !genome2genome
            && !cdna2genome
            && !coding2genome
            && !protein2genome
            && !protein2genome_bestfit
            && !est2genome
            && query_type == "dna"
            && target_type == "dna"
            && matches!(model, Model::Ungapped | Model::Local);
        let protein_dna = !ungapped_translated
            && !ner
            && !coding2coding
            && !genome2genome
            && !cdna2genome
            && !coding2genome
            && !protein2genome
            && !protein2genome_bestfit
            && !est2genome
            && query_type == "protein"
            && target_type == "dna"
            && matches!(model, Model::Ungapped | Model::Local | Model::BestFit);
        let protein_affine = !ungapped_translated
            && !ner
            && !coding2coding
            && !genome2genome
            && !cdna2genome
            && !coding2genome
            && !protein2genome
            && !protein2genome_bestfit
            && !est2genome
            && query_type == "protein"
            && target_type == "protein"
            && matches!(model, Model::Ungapped | Model::Local);
        if ungapped_translated {
            align_ungapped_translated_database_suboptimal(
                &q,
                &t,
                scoring,
                min_score.unwrap_or(100),
                both,
            )
        } else if ner {
            align_ner_database_suboptimal(
                &q,
                &t,
                scoring,
                min_ner,
                max_ner,
                ner_open,
                min_score.unwrap_or(100),
                both,
            )
        } else if est2genome {
            align_est2genome_database_suboptimal(&q, &t, intron, min_score.unwrap_or(100), both)
        } else if coding2coding {
            align_coding2coding_database_suboptimal(&q, &t, scoring, min_score.unwrap_or(100))
        } else if genome2genome {
            align_genome_to_genome_database_suboptimal_stranded(
                &q,
                &t,
                scoring,
                intron,
                min_score.unwrap_or(100),
                both,
            )
        } else if cdna2genome {
            align_cdna_to_genome_database_suboptimal(
                &q,
                &t,
                scoring,
                intron,
                min_score.unwrap_or(100),
            )
        } else if coding2genome {
            align_coding_to_genome_database_suboptimal(
                &q,
                &t,
                scoring,
                intron,
                min_score.unwrap_or(100),
                both,
            )
        } else if protein2genome || protein2genome_bestfit {
            align_protein_to_genome_database_suboptimal(
                &q,
                &t,
                scoring,
                intron,
                min_score.unwrap_or(100),
                both,
                protein2genome_bestfit,
            )
        } else if protein_dna {
            align_protein_to_dna_database_suboptimal(
                &q,
                &t,
                model,
                scoring,
                min_score.unwrap_or(100),
                both,
            )
        } else if protein_affine {
            if model == Model::Ungapped {
                align_protein_ungapped_database_suboptimal(
                    &q,
                    &t,
                    scoring,
                    min_score.unwrap_or(100),
                )
            } else {
                align_protein_database_suboptimal(&q, &t, scoring, min_score.unwrap_or(100))
            }
        } else if plain_dna {
            if model == Model::Ungapped {
                align_ungapped_database_suboptimal(&q, &t, scoring, min_score.unwrap_or(100), both)
            } else {
                align_database_suboptimal(&q, &t, scoring, min_score.unwrap_or(100), both)
            }
        } else {
            return Err("suboptimal enumeration is not implemented for this model".into());
        }
    } else if ungapped_translated {
        align_ungapped_translated_database(&q, &t, scoring, both)
    } else if ner {
        align_ner_database(&q, &t, scoring, min_ner, max_ner, ner_open, both)
    } else if coding2coding {
        align_coding2coding_database(&q, &t, scoring)
    } else if genome2genome {
        if exhaustive {
            align_genome_to_genome_database_stranded(&q, &t, scoring, intron, both, dp_memory_mb)
        } else if dp_memory_mb == 0 {
            align_genome_to_genome_database_heuristic_stranded_with_rolling_scores(
                &q, &t, scoring, intron, both, heuristic, true,
            )
        } else {
            align_genome_to_genome_database_heuristic_stranded(
                &q, &t, scoring, intron, both, heuristic,
            )
        }
    } else if cdna2genome {
        if exhaustive {
            align_cdna_to_genome_database_with_dp_memory_stranded(
                &q,
                &t,
                scoring,
                intron,
                both,
                dp_memory_mb,
            )
        } else {
            align_cdna_to_genome_database_heuristic(&q, &t, scoring, intron, heuristic)
        }
    } else if coding2genome {
        if exhaustive {
            align_coding_to_genome_database_with_dp_memory(
                &q,
                &t,
                scoring,
                intron,
                both,
                dp_memory_mb,
            )
        } else {
            align_coding_to_genome_database_heuristic(&q, &t, scoring, intron, both, heuristic)
        }
    } else if protein2genome_bestfit {
        if exhaustive {
            align_protein_to_genome_database_with_dp_memory(
                &q,
                &t,
                scoring,
                intron,
                both,
                true,
                dp_memory_mb,
            )
        } else {
            align_protein_to_genome_bestfit_database_heuristic(
                &q, &t, scoring, intron, both, heuristic,
            )
        }
    } else if protein2genome {
        if exhaustive {
            align_protein_to_genome_database_with_dp_memory(
                &q,
                &t,
                scoring,
                intron,
                both,
                false,
                dp_memory_mb,
            )
        } else {
            align_protein_to_genome_database_heuristic(&q, &t, scoring, intron, both, heuristic)
        }
    } else if est2genome {
        if exhaustive {
            align_est2genome_database_with_dp_memory(&q, &t, intron, both, dp_memory_mb)
        } else {
            align_est2genome_database_heuristic(&q, &t, intron, both, heuristic)
        }
    } else {
        match (query_type.as_str(), target_type.as_str()) {
            ("dna", "dna") if model == Model::Local && !exhaustive => {
                align_database_heuristic(&q, &t, model, scoring, both, heuristic)
            }
            ("dna", "dna") => {
                align_database_with_dp_memory(&q, &t, model, scoring, both, dp_memory_mb)
            }
            ("protein", "protein") => {
                align_protein_database_with_dp_memory(&q, &t, model, scoring, dp_memory_mb)
            }
            ("protein", "dna") => align_protein_to_dna_database(&q, &t, model, scoring, both),
            _ => {
                return Err(
                    "implemented type pairs are dna:dna, protein:protein, and protein:dna".into(),
                );
            }
        }
    };
    if let Some(threshold) = min_score {
        alignments.retain(|alignment| alignment.score >= threshold);
    }
    if let Some(percent) = percent {
        if !(0.0..=100.0).contains(&percent) {
            return Err("percent threshold must be between 0 and 100".into());
        }
        let self_scores: HashMap<_, _> = q
            .iter()
            .map(|sequence| {
                let score = if query_type == "protein" {
                    protein_self_score(sequence)
                } else if ungapped_translated || coding2coding || coding2genome {
                    translated_self_score(sequence)
                } else if cdna2genome || genome2genome {
                    dna_self_score(sequence).max(translated_self_score(sequence))
                } else {
                    dna_self_score(sequence)
                };
                (sequence.id.as_str(), score)
            })
            .collect();
        alignments.retain(|alignment| {
            let self_score = self_scores
                .get(alignment.query_id.as_str())
                .copied()
                .unwrap_or(0);
            f64::from(alignment.score) * 100.0 >= f64::from(self_score) * percent
        });
    }
    let best_n = best_n.filter(|&count| count > 0);
    if let Some(best_n) = best_n {
        alignments.sort_by(|left, right| {
            left.query_id
                .cmp(&right.query_id)
                .then_with(|| right.score.cmp(&left.score))
                .then_with(|| left.target_id.cmp(&right.target_id))
                .then_with(|| left.target_start.cmp(&right.target_start))
                .then_with(|| left.query_start.cmp(&right.query_start))
                .then_with(|| left.target_end.cmp(&right.target_end))
                .then_with(|| left.query_end.cmp(&right.query_end))
        });
        let mut seen: HashMap<String, (usize, Option<i32>)> = HashMap::new();
        alignments.retain(|alignment| {
            let (count, cutoff) = seen.entry(alignment.query_id.clone()).or_default();
            *count += 1;
            if *count <= best_n {
                if *count == best_n {
                    *cutoff = Some(alignment.score);
                }
                true
            } else {
                cutoff.is_some_and(|score| alignment.score == score)
            }
        });
    }
    if !forward_coordinates {
        let query_lengths: HashMap<_, _> = q
            .iter()
            .map(|sequence| (sequence.id.as_str(), sequence.bases.len() as u64))
            .collect();
        for alignment in &mut alignments {
            if alignment.query_strand == Strand::Reverse {
                let length = query_lengths
                    .get(alignment.query_id.as_str())
                    .copied()
                    .expect("alignment query comes from the input database");
                alignment.query_start = length - alignment.query_start;
                alignment.query_end = length - alignment.query_end;
            }
            if alignment.target_strand == Strand::Reverse {
                alignment.target_start = alignment.target_len - alignment.target_start;
                alignment.target_end = alignment.target_len - alignment.target_end;
            }
        }
    }
    let report_model = if ungapped_translated {
        "ungapped:trans".to_owned()
    } else if ner {
        "NER:affine:local:dna2dna".to_owned()
    } else if est2genome {
        "est2genome".to_owned()
    } else if coding2coding {
        "coding2coding".to_owned()
    } else if coding2genome {
        "coding2genome".to_owned()
    } else if cdna2genome {
        "cdna2genome".to_owned()
    } else if genome2genome {
        "genome2genome".to_owned()
    } else if protein2genome_bestfit {
        "protein2genome:bestfit".to_owned()
    } else if protein2genome {
        "protein2genome:local".to_owned()
    } else if matches!(model_name.as_str(), "protein2dna" | "p2d") {
        "protein2dna".to_owned()
    } else if matches!(model_name.as_str(), "protein2dna:bestfit" | "p2d:b") {
        "protein2dna:bestfit".to_owned()
    } else {
        let base = match model {
            Model::Ungapped => "ungapped",
            Model::Global => "affine:global",
            Model::BestFit => "affine:bestfit",
            Model::Local => "affine:local",
            Model::Overlap => "affine:overlap",
        };
        format!("{base}:{query_type}2{target_type}")
    };
    let mut output_ranks: HashMap<String, usize> = HashMap::new();
    for a in alignments {
        let rank = output_ranks.entry(a.query_id.clone()).or_default();
        *rank += 1;
        if show_alignment {
            if let (Some(query), Some(target)) = (
                q.iter().find(|sequence| sequence.id == a.query_id),
                t.iter().find(|sequence| sequence.id == a.target_id),
            ) {
                write_stdout(&pretty_alignment(&a, query, target, &report_model), false)?;
            }
        }
        if sugar {
            write_stdout(&a.sugar(), true)?
        }
        if cigar {
            let cigar_output = if cdna2genome || genome2genome {
                a.cdna_cigar()
            } else {
                a.cigar()
            };
            write_stdout(&cigar_output, true)?
        }
        if vulgar {
            write_stdout(&a.vulgar(), true)?
        }
        if gff {
            write_stdout(&a.gff3(), true)?
        }
        if query_gff {
            write_stdout(&a.query_gff3(), true)?
        }
        if let Some(format) = &ryo {
            let query = q
                .iter()
                .find(|sequence| sequence.id == a.query_id)
                .ok_or_else(|| format!("query record not found for {}", a.query_id))?;
            let target = t
                .iter()
                .find(|sequence| sequence.id == a.target_id)
                .ok_or_else(|| format!("target record not found for {}", a.target_id))?;
            let query_definition = query_definitions
                .get(&query.id)
                .map(String::as_str)
                .unwrap_or(&query.id);
            let target_definition = target_definitions
                .get(&target.id)
                .map(String::as_str)
                .unwrap_or(&target.id);
            let report_rank = if best_n.is_some() { *rank } else { 0 };
            let rendered = render_ryo(
                format,
                &a,
                query,
                target,
                query_definition,
                target_definition,
                &report_model,
                report_rank,
            )?;
            write_stdout(&rendered, false)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use exonerate_core::{Strand, align};

    #[test]
    fn common_ryo_fields_match_upstream_layout() {
        let query = Sequence {
            id: "q".into(),
            bases: b"ACGT".to_vec(),
        };
        let target = Sequence {
            id: "t".into(),
            bases: b"TTACGTGG".to_vec(),
        };
        let alignment = align(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        let rendered = render_ryo(
            "%qi %ti %qab %qae %tab %tae %s %r %m %S %V\n",
            &alignment,
            &query,
            &target,
            "q",
            "t",
            "affine:local:dna2dna",
            0,
        )
        .unwrap();
        assert_eq!(
            rendered,
            "q t 0 4 2 6 20 0 affine:local:dna2dna q 0 4 + t 2 6 + 20 M 4 4\n"
        );
    }

    #[test]
    fn ryo_rejects_unknown_tokens() {
        let sequence = Sequence {
            id: "s".into(),
            bases: b"A".to_vec(),
        };
        let alignment = align(
            &sequence,
            &sequence,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        assert!(render_ryo("%x", &alignment, &sequence, &sequence, "s", "s", "m", 0).is_err());
    }

    #[test]
    fn transition_ryo_preserves_atomic_scores_and_epsilon_edges() {
        let query = Sequence {
            id: "q".into(),
            bases: b"AC".to_vec(),
        };
        let target = Sequence {
            id: "t".into(),
            bases: b"ATC".to_vec(),
        };
        let alignment = align(
            &query,
            &target,
            Model::Global,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<i32>(),
            alignment.score
        );
        let rendered = render_ryo(
            "{%Pn|%Pl|%Ps|%Pqa|%Pta\n}",
            &alignment,
            &query,
            &target,
            "q",
            "t",
            "affine:global:dna2dna",
            0,
        )
        .unwrap();
        assert_eq!(
            rendered,
            concat!(
                "start to match|none|0|0|0\n",
                "match|match|5|1|1\n",
                "match to delete|gap|-12|0|1\n",
                "delete to match|none|0|0|0\n",
                "match|match|5|1|1\n",
                "match to end|none|0|0|0\n",
            )
        );
    }
}
