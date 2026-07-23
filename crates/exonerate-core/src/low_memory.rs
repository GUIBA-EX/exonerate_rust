//! Exact checkpointed traceback for affine dynamic programming.
//!
//! The forward pass keeps rolling score rows and periodic score checkpoints.
//! Traceback recomputes parent pointers for one checkpoint section at a time.

use super::*;

#[derive(Clone)]
struct ScoreRow {
    row: usize,
    m: Vec<Score>,
    i: Vec<Score>,
    d: Vec<Score>,
}

#[derive(Clone, Copy)]
struct EndPoint {
    i: usize,
    j: usize,
    state: State,
    score: Score,
}

fn valid_end(model: Model, i: usize, j: usize, n: usize, m: usize) -> bool {
    match model {
        Model::Global => i == n && j == m,
        Model::BestFit => i == n,
        Model::Local => true,
        Model::Overlap => i == n || j == m,
        Model::Ungapped => false,
    }
}

fn consider_end(end: &mut EndPoint, i: usize, j: usize, state: State, score: Score) {
    if score > end.score
        || (score == end.score && (i, j, state.rank()) > (end.i, end.j, end.state.rank()))
    {
        *end = EndPoint { i, j, state, score };
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_row(
    query: &Sequence,
    target: &[u8],
    row: usize,
    model: Model,
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
    previous: Option<(&[Score], &[Score], &[Score])>,
    current_m: &mut [Score],
    current_i: &mut [Score],
    current_d: &mut [Score],
    mut parents: Option<(&mut [State], &mut [State], &mut [State])>,
) {
    let local = model == Model::Local;
    let bestfit = model == Model::BestFit;
    let overlap = model == Model::Overlap;
    current_m.fill(NEG_INF);
    current_i.fill(NEG_INF);
    current_d.fill(NEG_INF);
    if let Some((pm, pi, pd)) = parents.as_mut() {
        pm.fill(State::Stop);
        pi.fill(State::Stop);
        pd.fill(State::Stop);
    }
    for j in 0..=target.len() {
        if row == 0 && j == 0 {
            current_m[j] = 0;
            continue;
        }
        if (local && (row == 0 || j == 0))
            || (bestfit && row == 0)
            || (overlap && (row == 0 || j == 0))
        {
            current_m[j] = 0;
        }
        if model == Model::Global && row == 0 && j > 0 {
            let (value, state) = best(&[
                (add(current_m[j - 1], scoring.gap_open), State::M),
                (add(current_d[j - 1], scoring.gap_extend), State::D),
            ]);
            current_d[j] = value;
            if let Some((_, _, pd)) = parents.as_mut() {
                pd[j] = state;
            }
        }
        if row > 0 {
            let (previous_m, previous_i, previous_d) = previous.expect("previous score row");
            if model == Model::Global && j == 0 {
                let (value, state) = best(&[
                    (add(previous_m[j], scoring.gap_open), State::M),
                    (add(previous_i[j], scoring.gap_extend), State::I),
                ]);
                current_i[j] = value;
                if let Some((_, pi, _)) = parents.as_mut() {
                    pi[j] = state;
                }
            }
            if j > 0 {
                let substitution = scorer(query.bases[row - 1], target[j - 1], scoring);
                let (mut value, mut state) = best(&[
                    (add(previous_m[j - 1], substitution), State::M),
                    (add(previous_i[j - 1], substitution), State::I),
                    (add(previous_d[j - 1], substitution), State::D),
                ]);
                if local && value < 0 {
                    value = 0;
                    state = State::Stop;
                }
                current_m[j] = value;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = state;
                }
            }
            if !(model == Model::Global && j == 0) {
                let (mut value, mut state) = best(&[
                    (add(previous_m[j], scoring.gap_open), State::M),
                    (add(previous_i[j], scoring.gap_extend), State::I),
                ]);
                if local && value < 0 {
                    value = NEG_INF;
                    state = State::Stop;
                }
                current_i[j] = value;
                if let Some((_, pi, _)) = parents.as_mut() {
                    pi[j] = state;
                }
            }
        }
        if j > 0 && !(model == Model::Global && row == 0) {
            let (mut value, mut state) = best(&[
                (add(current_m[j - 1], scoring.gap_open), State::M),
                (add(current_d[j - 1], scoring.gap_extend), State::D),
            ]);
            if local && value < 0 {
                value = NEG_INF;
                state = State::Stop;
            }
            current_d[j] = value;
            if let Some((_, _, pd)) = parents.as_mut() {
                pd[j] = state;
            }
        }
    }
}

fn checkpoint_stride(query_len: usize) -> usize {
    (query_len as f64).sqrt().ceil().max(1.0) as usize
}

fn score_checkpoints(
    query: &Sequence,
    target: &[u8],
    model: Model,
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
) -> (Vec<ScoreRow>, EndPoint, usize) {
    let (n, m) = (query.bases.len(), target.len());
    let stride = checkpoint_stride(n);
    let mut previous_m = vec![NEG_INF; m + 1];
    let mut previous_i = vec![NEG_INF; m + 1];
    let mut previous_d = vec![NEG_INF; m + 1];
    let mut current_m = vec![NEG_INF; m + 1];
    let mut current_i = vec![NEG_INF; m + 1];
    let mut current_d = vec![NEG_INF; m + 1];
    let mut checkpoints = Vec::new();
    let mut end = EndPoint {
        i: 0,
        j: 0,
        state: State::M,
        score: NEG_INF,
    };
    for row in 0..=n {
        fill_row(
            query,
            target,
            row,
            model,
            scoring,
            scorer,
            (row > 0).then_some((&previous_m, &previous_i, &previous_d)),
            &mut current_m,
            &mut current_i,
            &mut current_d,
            None,
        );
        for j in 0..=m {
            if valid_end(model, row, j, n, m) {
                consider_end(&mut end, row, j, State::M, current_m[j]);
                consider_end(&mut end, row, j, State::I, current_i[j]);
                consider_end(&mut end, row, j, State::D, current_d[j]);
            }
        }
        if row % stride == 0 || row == n {
            checkpoints.push(ScoreRow {
                row,
                m: current_m.clone(),
                i: current_i.clone(),
                d: current_d.clone(),
            });
        }
        std::mem::swap(&mut previous_m, &mut current_m);
        std::mem::swap(&mut previous_i, &mut current_i);
        std::mem::swap(&mut previous_d, &mut current_d);
    }
    (checkpoints, end, stride)
}

#[allow(clippy::too_many_arguments)]
fn rebuild_block(
    query: &Sequence,
    target: &[u8],
    model: Model,
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
    checkpoint: &ScoreRow,
    end_row: usize,
) -> (Vec<State>, Vec<State>, Vec<State>) {
    let cols = target.len() + 1;
    let rows = end_row - checkpoint.row;
    let mut pm = vec![State::Stop; (rows + 1) * cols];
    let mut pi = vec![State::Stop; (rows + 1) * cols];
    let mut pd = vec![State::Stop; (rows + 1) * cols];
    let (mut previous_m, mut previous_i, mut previous_d) = (
        checkpoint.m.clone(),
        checkpoint.i.clone(),
        checkpoint.d.clone(),
    );
    let mut current_m = vec![NEG_INF; cols];
    let mut current_i = vec![NEG_INF; cols];
    let mut current_d = vec![NEG_INF; cols];
    for row in checkpoint.row + 1..=end_row {
        let offset = (row - checkpoint.row) * cols;
        fill_row(
            query,
            target,
            row,
            model,
            scoring,
            scorer,
            Some((&previous_m, &previous_i, &previous_d)),
            &mut current_m,
            &mut current_i,
            &mut current_d,
            Some((
                &mut pm[offset..offset + cols],
                &mut pi[offset..offset + cols],
                &mut pd[offset..offset + cols],
            )),
        );
        std::mem::swap(&mut previous_m, &mut current_m);
        std::mem::swap(&mut previous_i, &mut current_i);
        std::mem::swap(&mut previous_d, &mut current_d);
    }
    (pm, pi, pd)
}

#[allow(clippy::too_many_arguments)]
fn finish_alignment(
    query: &Sequence,
    target: &Sequence,
    oriented_target: &[u8],
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
    strand: Strand,
    end: EndPoint,
    query_start: usize,
    target_start: usize,
    mut ops: Vec<(Op, u16)>,
) -> Alignment {
    ops.reverse();
    let mut raw_trace = Vec::new();
    if !ops.is_empty() {
        raw_trace.push(RawStep {
            transition_id: 0,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
        let (mut qi, mut tj) = (query_start, target_start);
        let mut previous_op = None;
        for &(op, transition_id) in &ops {
            if let Some(epsilon_id) = match previous_op {
                Some(Op::Insert) if op != Op::Insert => Some(6),
                Some(Op::Delete) if op != Op::Delete => Some(7),
                _ => None,
            } {
                raw_trace.push(RawStep {
                    transition_id: epsilon_id,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                });
            }
            let (query_advance, target_advance) = op.advances();
            let step_score = match op {
                Op::Match => scorer(query.bases[qi], oriented_target[tj], scoring),
                Op::Insert | Op::Delete if matches!(transition_id, 2 | 3) => scoring.gap_open,
                Op::Insert | Op::Delete => scoring.gap_extend,
                _ => unreachable!("affine traceback operation"),
            };
            raw_trace.push(RawStep {
                transition_id,
                query_advance,
                target_advance,
                score: step_score,
            });
            qi += query_advance as usize;
            tj += target_advance as usize;
            previous_op = Some(op);
        }
        if let Some(transition_id) = match previous_op {
            Some(Op::Insert) => Some(6),
            Some(Op::Delete) => Some(7),
            _ => None,
        } {
            raw_trace.push(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
        }
        raw_trace.push(RawStep {
            transition_id: 8,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    let mut trace: Vec<TraceRun> = Vec::new();
    for (op, transition_id) in ops {
        let (query_advance, target_advance) = op.advances();
        if let Some(last) = trace.last_mut()
            && last.op == op
            && last.transition_id == transition_id
            && last.query_advance == query_advance
            && last.target_advance == target_advance
        {
            last.repeats += 1;
        } else {
            trace.push(TraceRun {
                transition_id,
                op,
                query_advance,
                target_advance,
                repeats: 1,
            });
        }
    }
    let (target_start, target_end) = if strand == Strand::Forward {
        (target_start as u64, end.j as u64)
    } else {
        (
            target.bases.len() as u64 - target_start as u64,
            target.bases.len() as u64 - end.j as u64,
        )
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: query_start as u64,
        query_end: end.i as u64,
        query_strand: Strand::Forward,
        target_start,
        target_end,
        target_len: target.bases.len() as u64,
        target_strand: strand,
        score: end.score,
        raw_trace,
        trace,
    }
}

fn align_checkpointed_with_scorer(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
    scorer: fn(u8, u8, Scoring) -> Score,
) -> Alignment {
    if model == Model::Ungapped {
        return align_with_scorer(query, target, model, scoring, strand, scorer);
    }
    let oriented_target = if strand == Strand::Forward {
        target.bases.clone()
    } else {
        reverse_complement(&target.bases)
    };
    let (checkpoints, end, stride) =
        score_checkpoints(query, &oriented_target, model, scoring, scorer);
    let cols = oriented_target.len() + 1;
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let mut reversed_ops = Vec::new();
    while state != State::Stop {
        let block_start = if i == 0 {
            0
        } else {
            ((i - 1) / stride) * stride
        };
        let checkpoint = checkpoints
            .iter()
            .find(|checkpoint| checkpoint.row == block_start)
            .expect("checkpoint row");
        let (pm, pi, pd) = rebuild_block(
            query,
            &oriented_target,
            model,
            scoring,
            scorer,
            checkpoint,
            i,
        );
        loop {
            if state == State::Stop {
                break;
            }
            if i == block_start {
                if state == State::D && j > 0 {
                    let (_, previous) = best(&[
                        (add(checkpoint.m[j - 1], scoring.gap_open), State::M),
                        (add(checkpoint.d[j - 1], scoring.gap_extend), State::D),
                    ]);
                    reversed_ops.push((Op::Delete, if previous == State::M { 3 } else { 5 }));
                    j -= 1;
                    state = previous;
                    continue;
                }
                if i == 0 {
                    state = State::Stop;
                }
                break;
            }
            let k = (i - block_start) * cols + j;
            let current = state;
            let previous = match current {
                State::M => {
                    if j == 0 || pm[k] == State::Stop {
                        state = State::Stop;
                        break;
                    }
                    i -= 1;
                    j -= 1;
                    pm[k]
                }
                State::I => {
                    i -= 1;
                    pi[k]
                }
                State::D => {
                    if j == 0 {
                        state = State::Stop;
                        break;
                    }
                    j -= 1;
                    pd[k]
                }
                State::Stop => break,
            };
            let op = match current {
                State::M => Op::Match,
                State::I => Op::Insert,
                State::D => Op::Delete,
                State::Stop => unreachable!(),
            };
            let transition_id = match (current, previous) {
                (State::M, _) => 1,
                (State::I, State::M) => 2,
                (State::I, _) => 4,
                (State::D, State::M) => 3,
                (State::D, _) => 5,
                (State::Stop, _) => unreachable!(),
            };
            reversed_ops.push((op, transition_id));
            state = previous;
        }
    }
    finish_alignment(
        query,
        target,
        &oriented_target,
        scoring,
        scorer,
        strand,
        end,
        i,
        j,
        reversed_ops,
    )
}

/// Exact affine alignment using rolling score rows and checkpointed traceback.
pub fn align_low_memory(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
) -> Alignment {
    align_checkpointed_with_scorer(query, target, model, scoring, strand, dna_score)
}

/// Protein counterpart of [`align_low_memory`].
pub fn align_protein_low_memory(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
) -> Alignment {
    let mut alignment = align_checkpointed_with_scorer(
        query,
        target,
        model,
        scoring,
        Strand::Forward,
        protein_score,
    );
    alignment.query_strand = Strand::Unknown;
    alignment.target_strand = Strand::Unknown;
    alignment
}

fn full_affine_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * (3 * size_of::<Score>() as u128 + 3 * size_of::<State>() as u128)
}

fn use_checkpoints(query_len: usize, target_len: usize, memory_mb: usize) -> bool {
    full_affine_dp_bytes(query_len, target_len) > ((memory_mb as u128) << 20)
}

pub fn align_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
    memory_mb: usize,
) -> Alignment {
    if use_checkpoints(query.bases.len(), target.bases.len(), memory_mb) {
        align_low_memory(query, target, model, scoring, strand)
    } else {
        align(query, target, model, scoring, strand)
    }
}

pub fn align_protein_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    memory_mb: usize,
) -> Alignment {
    if use_checkpoints(query.bases.len(), target.bases.len(), memory_mb) {
        align_protein_low_memory(query, target, model, scoring)
    } else {
        align_protein(query, target, model, scoring)
    }
}

pub fn align_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    both_strands: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_with_dp_memory(
                query,
                target,
                model,
                scoring,
                Strand::Forward,
                memory_mb,
            ));
            if both_strands {
                let reverse_query = Sequence {
                    id: query.id.clone(),
                    bases: reverse_complement(&query.bases),
                };
                let mut alignment = align_with_dp_memory(
                    &reverse_query,
                    target,
                    model,
                    scoring,
                    Strand::Forward,
                    memory_mb,
                );
                alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                alignment.query_strand = Strand::Reverse;
                out.push(alignment);
            }
        }
    }
    out
}

pub fn align_protein_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    memory_mb: usize,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_protein_with_dp_memory(
                query, target, model, scoring, memory_mb,
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sequence(id: &str, bases: &str) -> Sequence {
        Sequence {
            id: id.into(),
            bases: bases.as_bytes().to_vec(),
        }
    }

    fn assert_same(left: &Alignment, right: &Alignment) {
        assert_eq!(left.query_id, right.query_id);
        assert_eq!(left.target_id, right.target_id);
        assert_eq!(
            (left.query_start, left.query_end),
            (right.query_start, right.query_end)
        );
        assert_eq!(
            (left.target_start, left.target_end),
            (right.target_start, right.target_end)
        );
        assert_eq!(
            (left.query_strand, left.target_strand),
            (right.query_strand, right.target_strand)
        );
        assert_eq!(left.score, right.score);
        assert_eq!(left.raw_trace, right.raw_trace);
        assert_eq!(left.trace, right.trace);
        assert_eq!(left.sugar(), right.sugar());
        assert_eq!(left.cigar(), right.cigar());
        assert_eq!(left.vulgar(), right.vulgar());
    }

    #[test]
    fn checkpointed_traceback_is_identical_for_all_affine_scopes() {
        let cases = [
            ("ACGTACGT", "TTACGTCGTAA"),
            ("AAAACCCCGGGG", "AAAAGGGG"),
            ("ACGT", "ACGT"),
            ("GATTACA", "GCATGCU"),
            ("", "ACGT"),
            ("ACGT", ""),
        ];
        for model in [Model::Global, Model::BestFit, Model::Local, Model::Overlap] {
            for strand in [Strand::Forward, Strand::Reverse] {
                for (query_bases, target_bases) in cases {
                    let query = sequence("q", query_bases);
                    let target = sequence("t", target_bases);
                    let full = align(&query, &target, model, Scoring::default(), strand);
                    let low = align_low_memory(&query, &target, model, Scoring::default(), strand);
                    assert_same(&low, &full);
                }
            }
        }
    }

    #[test]
    fn checkpointed_traceback_crosses_multiple_sections() {
        let query = sequence("q", &"ACGT".repeat(80));
        let mut target_bases = "ACGT".repeat(40);
        target_bases.push_str(&"TT".repeat(17));
        target_bases.push_str(&"ACGT".repeat(40));
        let target = sequence("t", &target_bases);
        let full = align(
            &query,
            &target,
            Model::Global,
            Scoring::default(),
            Strand::Forward,
        );
        let low = align_low_memory(
            &query,
            &target,
            Model::Global,
            Scoring::default(),
            Strand::Forward,
        );
        assert_same(&low, &full);
    }

    #[test]
    fn zero_memory_limit_forces_checkpointed_equivalent_path() {
        let query = sequence("q", "ACGTACGT");
        let target = sequence("t", "TTACGTACGTAA");
        assert!(use_checkpoints(query.bases.len(), target.bases.len(), 0));
        assert_same(
            &align_with_dp_memory(
                &query,
                &target,
                Model::Local,
                Scoring::default(),
                Strand::Forward,
                0,
            ),
            &align(
                &query,
                &target,
                Model::Local,
                Scoring::default(),
                Strand::Forward,
            ),
        );
    }

    #[test]
    fn checkpointed_protein_traceback_is_identical() {
        let query = sequence("q", "MTEYKLVVVGAGGVGKSALTIQLIQ");
        let target = sequence("t", "MTEYRLVVVGAGGVGKSALTIQLIQ");
        for model in [Model::Global, Model::BestFit, Model::Local, Model::Overlap] {
            assert_same(
                &align_protein_low_memory(&query, &target, model, Scoring::default()),
                &align_protein(&query, &target, model, Scoring::default()),
            );
        }
    }
}
