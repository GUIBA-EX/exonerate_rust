use exonerate_core::{
    Alignment, HeuristicConfig, IntronScoring, Model, RawStep, Scoring, Sequence, Strand,
    align_cdna_to_genome_database, align_cdna_to_genome_database_heuristic,
    align_cdna_to_genome_database_suboptimal, align_coding_to_genome_database,
    align_coding_to_genome_database_heuristic, align_coding_to_genome_database_suboptimal,
    align_coding2coding_database, align_coding2coding_database_suboptimal,
    align_database_heuristic, align_database_suboptimal, align_database_with_dp_memory,
    align_est2genome_database, align_est2genome_database_heuristic,
    align_est2genome_database_suboptimal, align_genome_to_genome_database,
    align_genome_to_genome_database_heuristic, align_genome_to_genome_database_suboptimal,
    align_ner_database, align_ner_database_suboptimal, align_protein_database_with_dp_memory,
    align_protein_to_dna_database, align_protein_to_dna_database_suboptimal,
    align_protein_to_genome_bestfit_database, align_protein_to_genome_bestfit_database_heuristic,
    align_protein_to_genome_database, align_protein_to_genome_database_heuristic,
    align_protein_to_genome_database_suboptimal, align_ungapped_translated_database,
    align_ungapped_translated_database_suboptimal, dna_self_score, protein_self_score, read_fasta,
    reverse_complement, translated_self_score,
};
use std::collections::HashMap;
use std::io::{self, Write};
use std::process::ExitCode;

fn usage() -> &'static str {
    "Usage: exonerate-rs [--model MODEL] [--querytype dna|protein] [--targettype dna|protein] [--gapopen N] [--gapextend N] [--codongapopen N] [--codongapextend N] [--frameshift N] [--minintron N] [--maxintron N] [--intronpenalty N] [--forcegtag yes|no] [--minner N] [--maxner N] [--neropen N] [--wordlen N] [--seedpadding N] [--seedrepeat N] [-D N|--dpmemory N] [--score N] [--percent N] [--bestn N] [--ryo FORMAT] [-q QUERY.fa] [-t TARGET.fa] [--subopt yes|no] [--exhaustive] [--revcomp yes|no] [--forwardonly] [--showsugar yes|no] [--showcigar yes|no] [--showvulgar yes|no] [--showgff yes|no] [--showquerygff yes|no] [--showtargetgff yes|no] QUERY.fa TARGET.fa\n\nImplemented models: ungapped, ungapped:trans, affine:global, affine:bestfit, affine:local, affine:overlap, coding2coding, coding2genome, cdna2genome, protein2dna, protein2dna:bestfit, protein2genome, protein2genome:bestfit, est2genome, genome2genome, ner"
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

fn sequence_text(sequence: &Sequence) -> String {
    String::from_utf8_lossy(&sequence.bases).into_owned()
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

fn report_tail(report: &str) -> String {
    report
        .split_whitespace()
        .skip(10)
        .collect::<Vec<_>>()
        .join(" ")
}

fn render_ryo_plain(
    format: &str,
    alignment: &Alignment,
    query: &Sequence,
    target: &Sequence,
    model: &str,
    rank: usize,
) -> Result<String, String> {
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
                    'i' | 'd' => output.push_str(&sequence.id),
                    'l' => output.push_str(&sequence.bases.len().to_string()),
                    's' => output.push_str(&sequence_text(sequence)),
                    'S' => output.push(strand_symbol(strand)),
                    't' => output.push_str(sequence_type(sequence)),
                    'a' | 'c' => {
                        let Some(&region_field) = chars.get(index) else {
                            return Err("incomplete aligned-region ryo token".into());
                        };
                        index += 1;
                        match region_field {
                            'b' => output.push_str(&start.to_string()),
                            'e' => output.push_str(&end.to_string()),
                            'l' => output.push_str(&start.abs_diff(end).to_string()),
                            's' => output
                                .push_str(&aligned_sequence_text(sequence, start, end, strand)),
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

fn render_ryo(
    format: &str,
    alignment: &Alignment,
    query: &Sequence,
    target: &Sequence,
    model: &str,
    rank: usize,
) -> Result<String, String> {
    let chars = format.chars().collect::<Vec<_>>();
    let Some(open) = chars.iter().position(|character| *character == '{') else {
        return render_ryo_plain(format, alignment, query, target, model, rank);
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
    let mut output = render_ryo_plain(&prefix, alignment, query, target, model, rank)?;
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
    output.push_str(&render_ryo(&suffix, alignment, query, target, model, rank)?);
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
    let mut args = std::env::args().skip(1);
    let mut model = Model::Ungapped;
    let mut model_name = "ungapped".to_owned();
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
    let mut query_type = "dna".to_owned();
    let mut target_type = "dna".to_owned();
    let (mut sugar, mut cigar, mut vulgar, mut gff) = (false, false, true, false);
    let mut query_gff = false;
    let mut min_score: Option<i32> = None;
    let mut percent: Option<f64> = None;
    let mut best_n: Option<usize> = None;
    let mut ryo: Option<String> = None;
    let (mut query_file, mut target_file): (Option<String>, Option<String>) = (None, None);
    let mut show_alignment = false;
    let mut subopt = false;
    let mut exhaustive = false;
    let mut dp_memory_mb = 32_usize;
    let mut heuristic = HeuristicConfig::default();
    let mut files = Vec::new();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-m" | "--model" => {
                let selected = args.next().ok_or("missing model")?;
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
            "-Q" | "--querytype" => query_type = args.next().ok_or("missing query type")?,
            "-T" | "--targettype" => target_type = args.next().ok_or("missing target type")?,
            "-E" | "--exhaustive" => exhaustive = true,
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
            "--forwardonly" => both = false,
            "--showsugar" => sugar = yes_no(&args.next().ok_or("missing showsugar value")?)?,
            "--showcigar" => cigar = yes_no(&args.next().ok_or("missing showcigar value")?)?,
            "--showvulgar" => vulgar = yes_no(&args.next().ok_or("missing showvulgar value")?)?,
            "--showgff" | "--showtargetgff" => {
                gff = yes_no(&args.next().ok_or("missing target GFF value")?)?
            }
            "--showquerygff" => query_gff = yes_no(&args.next().ok_or("missing query GFF value")?)?,
            "-h" | "--help" => {
                write_stdout(usage(), true)?;
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
    let q = read_fasta(&files[0]).map_err(|e| e.to_string())?;
    let t = read_fasta(&files[1]).map_err(|e| e.to_string())?;
    let mut alignments = if subopt {
        let plain_dna_local = !ungapped_translated
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
            && model == Model::Local;
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
            && matches!(model, Model::Local | Model::BestFit);
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
            align_genome_to_genome_database_suboptimal(
                &q,
                &t,
                scoring,
                intron,
                min_score.unwrap_or(100),
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
        } else if plain_dna_local {
            align_database_suboptimal(&q, &t, scoring, min_score.unwrap_or(100), both)
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
            align_genome_to_genome_database(&q, &t, scoring, intron)
        } else {
            align_genome_to_genome_database_heuristic(&q, &t, scoring, intron, heuristic)
        }
    } else if cdna2genome {
        if exhaustive {
            align_cdna_to_genome_database(&q, &t, scoring, intron)
        } else {
            align_cdna_to_genome_database_heuristic(&q, &t, scoring, intron, heuristic)
        }
    } else if coding2genome {
        if exhaustive {
            align_coding_to_genome_database(&q, &t, scoring, intron, both)
        } else {
            align_coding_to_genome_database_heuristic(&q, &t, scoring, intron, both, heuristic)
        }
    } else if protein2genome_bestfit {
        if exhaustive {
            align_protein_to_genome_bestfit_database(&q, &t, scoring, intron, both)
        } else {
            align_protein_to_genome_bestfit_database_heuristic(
                &q, &t, scoring, intron, both, heuristic,
            )
        }
    } else if protein2genome {
        if exhaustive {
            align_protein_to_genome_database(&q, &t, scoring, intron, both)
        } else {
            align_protein_to_genome_database_heuristic(&q, &t, scoring, intron, both, heuristic)
        }
    } else if est2genome {
        if exhaustive {
            align_est2genome_database(&q, &t, intron, both)
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
    if let Some(best_n) = best_n {
        alignments.sort_by(|left, right| {
            left.query_id
                .cmp(&right.query_id)
                .then_with(|| right.score.cmp(&left.score))
                .then_with(|| left.target_id.cmp(&right.target_id))
                .then_with(|| left.target_start.cmp(&right.target_start))
        });
        let mut seen: HashMap<String, usize> = HashMap::new();
        alignments.retain(|alignment| {
            let count = seen.entry(alignment.query_id.clone()).or_default();
            let keep = *count < best_n;
            *count += 1;
            keep
        });
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
                write_stdout(&a.pretty(query, target), false)?;
            }
        }
        if sugar {
            write_stdout(&a.sugar(), true)?
        }
        if cigar {
            write_stdout(&a.cigar(), true)?
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
            let report_rank = if best_n.is_some() { *rank } else { 0 };
            let rendered = render_ryo(format, &a, query, target, &report_model, report_rank)?;
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
        assert!(render_ryo("%x", &alignment, &sequence, &sequence, "m", 0).is_err());
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
