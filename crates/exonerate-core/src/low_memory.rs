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

type EstParentSlices<'a> = (
    &'a mut [Option<EstParent>],
    &'a mut [Option<EstParent>],
    &'a mut [Option<EstParent>],
);
type EstParentRows = (
    Vec<Option<EstParent>>,
    Vec<Option<EstParent>>,
    Vec<Option<EstParent>>,
);
type P2ParentSlices<'a> = (
    &'a mut [Option<P2GenomeParent>],
    &'a mut [Option<P2GenomeParent>],
    &'a mut [Option<P2GenomeParent>],
);
type P2ParentRows = (
    Vec<Option<P2GenomeParent>>,
    Vec<Option<P2GenomeParent>>,
    Vec<Option<P2GenomeParent>>,
);

#[derive(Clone)]
struct ScoreCheckpoint {
    row: usize,
    /// The checkpoint row together with its five immediate predecessors.
    /// Coding frameshifts can consume up to five query bases in one edge.
    history: Vec<ScoreRow>,
}

#[derive(Clone)]
struct CdnaScoreRow {
    row: usize,
    u5m: Vec<Score>,
    u5i: Vec<Score>,
    u5d: Vec<Score>,
    cm: Vec<Score>,
    ci: Vec<Score>,
    cd: Vec<Score>,
    u3m: Vec<Score>,
    u3i: Vec<Score>,
    u3d: Vec<Score>,
}

#[derive(Clone)]
struct CdnaCheckpoint {
    row: usize,
    history: Vec<CdnaScoreRow>,
}

/// One rolling row of the full genome-to-genome composed graph.  The two UTR
/// triples are separate because the CDS epsilon exit is one-way.
#[allow(dead_code)] // consumed by the genome2genome checkpoint forward pass
#[derive(Clone)]
pub(crate) struct GenomeScoreRow {
    pub(crate) row: usize,
    pub(crate) u5m: Vec<Score>,
    pub(crate) u5i: Vec<Score>,
    pub(crate) u5d: Vec<Score>,
    pub(crate) cm: Vec<Score>,
    pub(crate) ci: Vec<Score>,
    pub(crate) cd: Vec<Score>,
    pub(crate) u3m: Vec<Score>,
    pub(crate) u3i: Vec<Score>,
    pub(crate) u3d: Vec<Score>,
}

/// A self-contained genome-to-genome reconstruction boundary.  Unlike the
/// other intron models, the live candidate deques carry information from an
/// unbounded number of earlier rows and must be saved alongside scores.
#[allow(dead_code)] // consumed by the genome2genome checkpoint forward pass
#[derive(Clone)]
pub(crate) struct GenomeCheckpoint {
    pub(crate) row: usize,
    pub(crate) history: Vec<GenomeScoreRow>,
    pub(crate) queues: GenomeIntronQueues,
}

#[allow(dead_code)] // consumed by genome2genome block reconstruction
struct GenomeParentRows {
    u5m: Vec<Option<GenomeParent>>,
    u5i: Vec<Option<GenomeParent>>,
    u5d: Vec<Option<GenomeParent>>,
    cm: Vec<Option<GenomeParent>>,
    ci: Vec<Option<GenomeParent>>,
    cd: Vec<Option<GenomeParent>>,
    u3m: Vec<Option<GenomeParent>>,
    u3i: Vec<Option<GenomeParent>>,
    u3d: Vec<Option<GenomeParent>>,
}

#[allow(dead_code)] // consumed by genome2genome block reconstruction
impl GenomeParentRows {
    fn new(cells: usize) -> Self {
        Self {
            u5m: vec![None; cells],
            u5i: vec![None; cells],
            u5d: vec![None; cells],
            cm: vec![None; cells],
            ci: vec![None; cells],
            cd: vec![None; cells],
            u3m: vec![None; cells],
            u3i: vec![None; cells],
            u3d: vec![None; cells],
        }
    }

    fn get(&self, state: GenomeState, index: usize) -> Option<GenomeParent> {
        match state {
            GenomeState::UtrM => self.u5m[index],
            GenomeState::UtrI => self.u5i[index],
            GenomeState::UtrD => self.u5d[index],
            GenomeState::CdsM => self.cm[index],
            GenomeState::CdsI => self.ci[index],
            GenomeState::CdsD => self.cd[index],
            GenomeState::U3M => self.u3m[index],
            GenomeState::U3I => self.u3i[index],
            GenomeState::U3D => self.u3d[index],
        }
    }
}

#[allow(dead_code)] // consumed by the genome2genome checkpoint forward pass
fn empty_genome_score_row(row: usize, cols: usize) -> GenomeScoreRow {
    GenomeScoreRow {
        row,
        // The pre-CDS UTR is a local graph and therefore may start at every
        // query/target coordinate; all later graph portions remain unreachable
        // until an explicit transition enters them.
        u5m: vec![0; cols],
        u5i: vec![NEG_INF; cols],
        u5d: vec![NEG_INF; cols],
        cm: vec![NEG_INF; cols],
        ci: vec![NEG_INF; cols],
        cd: vec![NEG_INF; cols],
        u3m: vec![NEG_INF; cols],
        u3i: vec![NEG_INF; cols],
        u3d: vec![NEG_INF; cols],
    }
}

#[derive(Clone, Copy)]
struct CdnaEndPoint {
    i: usize,
    j: usize,
    state: CdnaState,
    score: Score,
}

struct CdnaParentRows {
    u5m: Vec<Option<CdnaParent>>,
    u5i: Vec<Option<CdnaParent>>,
    u5d: Vec<Option<CdnaParent>>,
    cm: Vec<Option<CdnaParent>>,
    ci: Vec<Option<CdnaParent>>,
    cd: Vec<Option<CdnaParent>>,
    u3m: Vec<Option<CdnaParent>>,
    u3i: Vec<Option<CdnaParent>>,
    u3d: Vec<Option<CdnaParent>>,
}

impl CdnaParentRows {
    fn new(cells: usize) -> Self {
        Self {
            u5m: vec![None; cells],
            u5i: vec![None; cells],
            u5d: vec![None; cells],
            cm: vec![None; cells],
            ci: vec![None; cells],
            cd: vec![None; cells],
            u3m: vec![None; cells],
            u3i: vec![None; cells],
            u3d: vec![None; cells],
        }
    }

    fn fill_none(&mut self) {
        self.u5m.fill(None);
        self.u5i.fill(None);
        self.u5d.fill(None);
        self.cm.fill(None);
        self.ci.fill(None);
        self.cd.fill(None);
        self.u3m.fill(None);
        self.u3i.fill(None);
        self.u3d.fill(None);
    }

    fn get(&self, state: CdnaState, index: usize) -> Option<CdnaParent> {
        match state {
            CdnaState::U5M => self.u5m[index],
            CdnaState::U5I => self.u5i[index],
            CdnaState::U5D => self.u5d[index],
            CdnaState::CodingM => self.cm[index],
            CdnaState::CodingI => self.ci[index],
            CdnaState::CodingD => self.cd[index],
            CdnaState::U3M => self.u3m[index],
            CdnaState::U3I => self.u3i[index],
            CdnaState::U3D => self.u3d[index],
        }
    }
}

fn cdna_row_relax(
    row: &mut CdnaScoreRow,
    parents: &mut Option<&mut CdnaParentRows>,
    destination: CdnaState,
    index: usize,
    candidate: Score,
    source: CdnaState,
    fragment: TraceFragment,
) {
    let score = match destination {
        CdnaState::U5M => &mut row.u5m[index],
        CdnaState::U5I => &mut row.u5i[index],
        CdnaState::U5D => &mut row.u5d[index],
        CdnaState::CodingM => &mut row.cm[index],
        CdnaState::CodingI => &mut row.ci[index],
        CdnaState::CodingD => &mut row.cd[index],
        CdnaState::U3M => &mut row.u3m[index],
        CdnaState::U3I => &mut row.u3i[index],
        CdnaState::U3D => &mut row.u3d[index],
    };
    if candidate < *score || (candidate == *score && fragment.atoms.iter().all(Option::is_none)) {
        return;
    }
    *score = candidate;
    if let Some(parents) = parents.as_deref_mut() {
        let parent = Some(CdnaParent {
            state: source,
            fragment,
        });
        match destination {
            CdnaState::U5M => parents.u5m[index] = parent,
            CdnaState::U5I => parents.u5i[index] = parent,
            CdnaState::U5D => parents.u5d[index] = parent,
            CdnaState::CodingM => parents.cm[index] = parent,
            CdnaState::CodingI => parents.ci[index] = parent,
            CdnaState::CodingD => parents.cd[index] = parent,
            CdnaState::U3M => parents.u3m[index] = parent,
            CdnaState::U3I => parents.u3i[index] = parent,
            CdnaState::U3D => parents.u3d[index] = parent,
        }
    }
}

fn checkpoint_count(query_len: usize, stride: usize) -> usize {
    query_len.div_ceil(stride.max(1)) + 1
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

#[allow(clippy::too_many_arguments)]
fn fill_est_row(
    query: &Sequence,
    bases: &[u8],
    row: usize,
    intron: IntronScoring,
    previous: Option<(&[Score], &[Score], &[Score])>,
    current_m: &mut [Score],
    current_i: &mut [Score],
    current_d: &mut [Score],
    mut parents: Option<EstParentSlices<'_>>,
) {
    current_m.fill(0);
    current_i.fill(NEG_INF);
    current_d.fill(NEG_INF);
    if let Some((pm, pi, pd)) = parents.as_mut() {
        pm.fill(None);
        pi.fill(None);
        pd.fill(None);
    }
    let mut window = IntronCandidateWindow::default();
    for j in 0..=bases.len() {
        if j >= intron.min_len as usize {
            let start = j - intron.min_len as usize;
            let (source, source_state) = best(&[
                (current_m[start], State::M),
                (current_i[start], State::I),
                (current_d[start], State::D),
            ]);
            if source > 0 {
                if let Some(donor) =
                    splice_score(bases, start, SpliceType::DonorForward, intron.force_gtag)
                {
                    window.insert(IntronCandidate {
                        start,
                        score: add(source, intron.open_penalty.saturating_add(donor)),
                        state_rank: source_state.rank(),
                    });
                }
            }
        }
        if j > intron.max_len as usize {
            window.expire_before(j - intron.max_len as usize);
        }
        if j >= 2 {
            if let (Some(candidate), Some(acceptor)) = (
                window.best(),
                splice_score(bases, j - 2, SpliceType::AcceptorForward, intron.force_gtag),
            ) {
                let score = add(candidate.score, acceptor);
                if score > current_m[j] {
                    let length = (j - candidate.start) as u32;
                    current_m[j] = score;
                    if let Some((pm, _, _)) = parents.as_mut() {
                        pm[j] = Some(EstParent {
                            state: state_from_rank(candidate.state_rank),
                            fragment: TraceFragment::intron(
                                TraceRun {
                                    transition_id: 6,
                                    op: Op::Splice5,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 7,
                                    op: Op::Intron,
                                    query_advance: 0,
                                    target_advance: length - 4,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 8,
                                    op: Op::Splice3,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                            ),
                            raw_fragment: RawFragment::with_prefix(
                                match state_from_rank(candidate.state_rank) {
                                    State::I => Some(9),
                                    State::D => Some(10),
                                    _ => None,
                                },
                                &[
                                    RawStep {
                                        transition_id: 6,
                                        query_advance: 0,
                                        target_advance: 2,
                                        score: intron.open_penalty.saturating_add(
                                            splice_score(
                                                bases,
                                                candidate.start,
                                                SpliceType::DonorForward,
                                                intron.force_gtag,
                                            )
                                            .unwrap_or(0),
                                        ),
                                    },
                                    RawStep {
                                        transition_id: 7,
                                        query_advance: 0,
                                        target_advance: length - 4,
                                        score: 0,
                                    },
                                    RawStep {
                                        transition_id: 8,
                                        query_advance: 0,
                                        target_advance: 2,
                                        score: acceptor,
                                    },
                                ],
                            ),
                        });
                    }
                }
            }
        }
        if row > 0 && j > 0 {
            let (previous_m, previous_i, previous_d) = previous.expect("previous EST row");
            let (score, state) = best(&[
                (previous_m[j - 1], State::M),
                (previous_i[j - 1], State::I),
                (previous_d[j - 1], State::D),
            ]);
            let score = add(
                score,
                dna_score(query.bases[row - 1], bases[j - 1], Scoring::default()),
            );
            if score > current_m[j] {
                current_m[j] = score;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = Some(EstParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 1,
                            op: Op::Match,
                            query_advance: 1,
                            target_advance: 1,
                            repeats: 1,
                        }),
                        raw_fragment: RawFragment::with_prefix(
                            match state {
                                State::I => Some(9),
                                State::D => Some(10),
                                _ => None,
                            },
                            &[RawStep {
                                transition_id: 1,
                                query_advance: 1,
                                target_advance: 1,
                                score: dna_score(
                                    query.bases[row - 1],
                                    bases[j - 1],
                                    Scoring::default(),
                                ),
                            }],
                        ),
                    });
                }
            }
        }
        if row > 0 {
            let (previous_m, previous_i, _) = previous.expect("previous EST row");
            let (score, state, transition_id) = if add(previous_m[j], -12) >= add(previous_i[j], -4)
            {
                (add(previous_m[j], -12), State::M, 2)
            } else {
                (add(previous_i[j], -4), State::I, 4)
            };
            if score >= 0 {
                current_i[j] = score;
                if let Some((_, pi, _)) = parents.as_mut() {
                    pi[j] = Some(EstParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id,
                            op: Op::Insert,
                            query_advance: 1,
                            target_advance: 0,
                            repeats: 1,
                        }),
                        raw_fragment: RawFragment::one(RawStep {
                            transition_id,
                            query_advance: 1,
                            target_advance: 0,
                            score: if transition_id == 2 { -12 } else { -4 },
                        }),
                    });
                }
            }
        }
        if j > 0 {
            let (score, state, transition_id) =
                if add(current_m[j - 1], -12) >= add(current_d[j - 1], -4) {
                    (add(current_m[j - 1], -12), State::M, 3)
                } else {
                    (add(current_d[j - 1], -4), State::D, 5)
                };
            if score >= 0 {
                current_d[j] = score;
                if let Some((_, _, pd)) = parents.as_mut() {
                    pd[j] = Some(EstParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 1,
                            repeats: 1,
                        }),
                        raw_fragment: RawFragment::one(RawStep {
                            transition_id,
                            query_advance: 0,
                            target_advance: 1,
                            score: if transition_id == 3 { -12 } else { -4 },
                        }),
                    });
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_p2g_row(
    query: &Sequence,
    bases: &[u8],
    row: usize,
    scoring: Scoring,
    intron: IntronScoring,
    bestfit: bool,
    previous: (&[Score], &[Score], &[Score]),
    current_m: &mut [Score],
    current_i: &mut [Score],
    current_d: &mut [Score],
    mut parents: Option<P2ParentSlices<'_>>,
) {
    current_m.fill(if bestfit { NEG_INF } else { 0 });
    current_i.fill(NEG_INF);
    current_d.fill(NEG_INF);
    if let Some((pm, pi, pd)) = parents.as_mut() {
        pm.fill(None);
        pi.fill(None);
        pd.fill(None);
    }
    let (previous_m, previous_i, previous_d) = previous;
    let mut phase_windows = [
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
    ];
    for j in 0..=bases.len() {
        for (phase, window) in phase_windows.iter_mut().enumerate() {
            let post_len = 3 - phase;
            if j < post_len + intron.min_len as usize + phase {
                continue;
            }
            let post_start = j - post_len;
            let donor = post_start - intron.min_len as usize;
            let source_j = donor - phase;
            let (source, state) = best(&[
                (previous_m[source_j], State::M),
                (previous_i[source_j], State::I),
                (previous_d[source_j], State::D),
            ]);
            if source > 0 || bestfit {
                if let Some(donor_score) =
                    splice_score(bases, donor, SpliceType::DonorForward, intron.force_gtag)
                {
                    window.insert(IntronCandidate {
                        start: donor,
                        score: add(source, intron.open_penalty.saturating_add(donor_score)),
                        state_rank: state.rank(),
                    });
                }
            }
            if post_start > intron.max_len as usize {
                window.expire_before(post_start - intron.max_len as usize);
            }
            let Some(candidate) = window.best() else {
                continue;
            };
            if post_start < 2 {
                continue;
            }
            let Some(acceptor) = splice_score(
                bases,
                post_start - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            ) else {
                continue;
            };
            let pre_start = candidate.start - phase;
            let mut codon = Vec::with_capacity(3);
            codon.extend_from_slice(&bases[pre_start..candidate.start]);
            codon.extend_from_slice(&bases[post_start..j]);
            if codon.len() != 3 {
                continue;
            }
            let score = add(
                add(candidate.score, acceptor),
                protein_score(query.bases[row - 1], translate_dna(&codon, 0)[0], scoring),
            );
            if score > current_m[j] {
                let intron_len = (post_start - candidate.start) as u32;
                let splice5 = TraceRun {
                    transition_id: 20,
                    op: Op::Splice5,
                    query_advance: 0,
                    target_advance: 2,
                    repeats: 1,
                };
                let loop_run = TraceRun {
                    transition_id: 21,
                    op: Op::Intron,
                    query_advance: 0,
                    target_advance: intron_len - 4,
                    repeats: 1,
                };
                let splice3 = TraceRun {
                    transition_id: 22,
                    op: Op::Splice3,
                    query_advance: 0,
                    target_advance: 2,
                    repeats: 1,
                };
                let fragment = if phase == 0 {
                    TraceFragment::intron_match(
                        splice5,
                        loop_run,
                        splice3,
                        TraceRun {
                            transition_id: 1,
                            op: Op::Match,
                            query_advance: 1,
                            target_advance: 3,
                            repeats: 1,
                        },
                    )
                } else {
                    TraceFragment::phase_intron(
                        TraceRun {
                            transition_id: 23,
                            op: Op::SplitCodon,
                            query_advance: 0,
                            target_advance: phase as u32,
                            repeats: 1,
                        },
                        splice5,
                        loop_run,
                        splice3,
                        TraceRun {
                            transition_id: 24,
                            op: Op::SplitCodon,
                            query_advance: 1,
                            target_advance: post_len as u32,
                            repeats: 1,
                        },
                    )
                };
                current_m[j] = score;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = Some(P2GenomeParent {
                        state: state_from_rank(candidate.state_rank),
                        fragment,
                    });
                }
            }
        }
        if j >= 3 {
            let from_m = add(current_m[j - 3], scoring.codon_gap_open);
            let from_d = add(current_d[j - 3], scoring.codon_gap_extend);
            if from_m >= from_d && (from_m >= 0 || bestfit) {
                current_d[j] = from_m;
                if let Some((_, _, pd)) = parents.as_mut() {
                    pd[j] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 3,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            } else if from_d >= 0 || bestfit {
                current_d[j] = from_d;
                if let Some((_, _, pd)) = parents.as_mut() {
                    pd[j] = Some(P2GenomeParent {
                        state: State::D,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 5,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
        for advance in [1_usize, 2, 4, 5] {
            if j < advance {
                continue;
            }
            let (source, state) = best(&[
                (current_m[j - advance], State::M),
                (current_i[j - advance], State::I),
                (current_d[j - advance], State::D),
            ]);
            let value = add(source, scoring.frameshift);
            if value > current_m[j] {
                current_m[j] = value;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = Some(P2GenomeParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 6,
                            op: Op::Frameshift,
                            query_advance: 0,
                            target_advance: advance as u32,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
        if j >= 3 {
            let (source, state) = best(&[
                (previous_m[j - 3], State::M),
                (previous_i[j - 3], State::I),
                (previous_d[j - 3], State::D),
            ]);
            let value = add(
                source,
                protein_score(
                    query.bases[row - 1],
                    translate_dna(&bases[j - 3..j], 0)[0],
                    scoring,
                ),
            );
            if value >= current_m[j] {
                current_m[j] = value;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = Some(P2GenomeParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 1,
                            op: Op::Match,
                            query_advance: 1,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
        let from_m = add(previous_m[j], scoring.codon_gap_open);
        let from_i = add(previous_i[j], scoring.codon_gap_extend);
        if from_m >= from_i && (from_m >= 0 || bestfit) {
            current_i[j] = from_m;
            if let Some((_, pi, _)) = parents.as_mut() {
                pi[j] = Some(P2GenomeParent {
                    state: State::M,
                    fragment: TraceFragment::one(TraceRun {
                        transition_id: 2,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                });
            }
        } else if from_i >= 0 || bestfit {
            current_i[j] = from_i;
            if let Some((_, pi, _)) = parents.as_mut() {
                pi[j] = Some(P2GenomeParent {
                    state: State::I,
                    fragment: TraceFragment::one(TraceRun {
                        transition_id: 4,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                });
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_c2g_row(
    query: &Sequence,
    bases: &[u8],
    row: usize,
    scoring: Scoring,
    intron: IntronScoring,
    history: &[ScoreRow],
    current_m: &mut [Score],
    current_i: &mut [Score],
    current_d: &mut [Score],
    mut parents: Option<P2ParentSlices<'_>>,
) {
    current_m.fill(0);
    current_i.fill(NEG_INF);
    current_d.fill(NEG_INF);
    if let Some((pm, pi, pd)) = parents.as_mut() {
        pm.fill(None);
        pi.fill(None);
        pd.fill(None);
    }
    let prior = |advance| {
        history
            .iter()
            .find(|scores| scores.row == row - advance)
            .expect("coding2genome rolling history")
    };
    let mut phase_windows = [
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
    ];
    for j in 0..=bases.len() {
        if row >= 3 {
            let previous = prior(3);
            for (phase, window) in phase_windows.iter_mut().enumerate() {
                let post_len = 3 - phase;
                if j < post_len + intron.min_len as usize + phase {
                    continue;
                }
                let post_start = j - post_len;
                let donor = post_start - intron.min_len as usize;
                let source_j = donor - phase;
                let (source, state) = best(&[
                    (previous.m[source_j], State::M),
                    (previous.i[source_j], State::I),
                    (previous.d[source_j], State::D),
                ]);
                if source > 0 {
                    if let Some(donor_score) =
                        splice_score(bases, donor, SpliceType::DonorForward, intron.force_gtag)
                    {
                        window.insert(IntronCandidate {
                            start: donor,
                            score: add(source, intron.open_penalty.saturating_add(donor_score)),
                            state_rank: state.rank(),
                        });
                    }
                }
                if post_start > intron.max_len as usize {
                    window.expire_before(post_start - intron.max_len as usize);
                }
                let Some(candidate) = window.best() else {
                    continue;
                };
                if post_start < 2 {
                    continue;
                }
                let Some(acceptor) = splice_score(
                    bases,
                    post_start - 2,
                    SpliceType::AcceptorForward,
                    intron.force_gtag,
                ) else {
                    continue;
                };
                let pre_start = candidate.start - phase;
                let mut codon = Vec::with_capacity(3);
                codon.extend_from_slice(&bases[pre_start..candidate.start]);
                codon.extend_from_slice(&bases[post_start..j]);
                if codon.len() != 3 {
                    continue;
                }
                let score = add(
                    add(candidate.score, acceptor),
                    protein_score(
                        translate_dna(&query.bases[row - 3..row], 0)[0],
                        translate_dna(&codon, 0)[0],
                        scoring,
                    ),
                );
                if score > current_m[j] {
                    let intron_len = (post_start - candidate.start) as u32;
                    let splice5 = TraceRun {
                        transition_id: 20,
                        op: Op::Splice5,
                        query_advance: 0,
                        target_advance: 2,
                        repeats: 1,
                    };
                    let loop_run = TraceRun {
                        transition_id: 21,
                        op: Op::Intron,
                        query_advance: 0,
                        target_advance: intron_len - 4,
                        repeats: 1,
                    };
                    let splice3 = TraceRun {
                        transition_id: 22,
                        op: Op::Splice3,
                        query_advance: 0,
                        target_advance: 2,
                        repeats: 1,
                    };
                    let fragment = if phase == 0 {
                        TraceFragment::intron_match(
                            splice5,
                            loop_run,
                            splice3,
                            TraceRun {
                                transition_id: 1,
                                op: Op::Match,
                                query_advance: 3,
                                target_advance: 3,
                                repeats: 1,
                            },
                        )
                    } else {
                        TraceFragment::phase_intron(
                            TraceRun {
                                transition_id: 23,
                                op: Op::SplitCodon,
                                query_advance: phase as u32,
                                target_advance: phase as u32,
                                repeats: 1,
                            },
                            splice5,
                            loop_run,
                            splice3,
                            TraceRun {
                                transition_id: 24,
                                op: Op::SplitCodon,
                                query_advance: post_len as u32,
                                target_advance: post_len as u32,
                                repeats: 1,
                            },
                        )
                    };
                    current_m[j] = score;
                    if let Some((pm, _, _)) = parents.as_mut() {
                        pm[j] = Some(P2GenomeParent {
                            state: state_from_rank(candidate.state_rank),
                            fragment,
                        });
                    }
                }
            }
        }
        if j >= 3 {
            let from_m = add(current_m[j - 3], scoring.codon_gap_open);
            let from_d = add(current_d[j - 3], scoring.codon_gap_extend);
            if from_m >= from_d && from_m >= 0 {
                current_d[j] = from_m;
                if let Some((_, _, pd)) = parents.as_mut() {
                    pd[j] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 3,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            } else if from_d >= 0 {
                current_d[j] = from_d;
                if let Some((_, _, pd)) = parents.as_mut() {
                    pd[j] = Some(P2GenomeParent {
                        state: State::D,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 5,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
        for advance in [1_usize, 2, 4, 5] {
            if row >= advance {
                let source = prior(advance);
                let (value, state) = best(&[
                    (source.m[j], State::M),
                    (source.i[j], State::I),
                    (source.d[j], State::D),
                ]);
                let value = add(value, scoring.frameshift);
                if value > current_m[j] {
                    current_m[j] = value;
                    if let Some((pm, _, _)) = parents.as_mut() {
                        pm[j] = Some(P2GenomeParent {
                            state,
                            fragment: TraceFragment::one(TraceRun {
                                transition_id: 6,
                                op: Op::Frameshift,
                                query_advance: advance as u32,
                                target_advance: 0,
                                repeats: 1,
                            }),
                        });
                    }
                }
            }
            if j >= advance {
                let (value, state) = best(&[
                    (current_m[j - advance], State::M),
                    (current_i[j - advance], State::I),
                    (current_d[j - advance], State::D),
                ]);
                let value = add(value, scoring.frameshift);
                if value > current_m[j] {
                    current_m[j] = value;
                    if let Some((pm, _, _)) = parents.as_mut() {
                        pm[j] = Some(P2GenomeParent {
                            state,
                            fragment: TraceFragment::one(TraceRun {
                                transition_id: 7,
                                op: Op::Frameshift,
                                query_advance: 0,
                                target_advance: advance as u32,
                                repeats: 1,
                            }),
                        });
                    }
                }
            }
        }
        if row >= 3 && j >= 3 {
            let source = prior(3);
            let (value, state) = best(&[
                (source.m[j - 3], State::M),
                (source.i[j - 3], State::I),
                (source.d[j - 3], State::D),
            ]);
            let value = add(
                value,
                protein_score(
                    translate_dna(&query.bases[row - 3..row], 0)[0],
                    translate_dna(&bases[j - 3..j], 0)[0],
                    scoring,
                ),
            );
            if value >= current_m[j] {
                current_m[j] = value;
                if let Some((pm, _, _)) = parents.as_mut() {
                    pm[j] = Some(P2GenomeParent {
                        state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 1,
                            op: Op::Match,
                            query_advance: 3,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
        if row >= 3 {
            let source = prior(3);
            let from_m = add(source.m[j], scoring.codon_gap_open);
            let from_i = add(source.i[j], scoring.codon_gap_extend);
            if from_m >= from_i && from_m >= 0 {
                current_i[j] = from_m;
                if let Some((_, pi, _)) = parents.as_mut() {
                    pi[j] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 2,
                            op: Op::Insert,
                            query_advance: 3,
                            target_advance: 0,
                            repeats: 1,
                        }),
                    });
                }
            } else if from_i >= 0 {
                current_i[j] = from_i;
                if let Some((_, pi, _)) = parents.as_mut() {
                    pi[j] = Some(P2GenomeParent {
                        state: State::I,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 4,
                            op: Op::Insert,
                            query_advance: 3,
                            target_advance: 0,
                            repeats: 1,
                        }),
                    });
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_cdna_row(
    query: &Sequence,
    target: &Sequence,
    row_index: usize,
    scoring: Scoring,
    intron: IntronScoring,
    history: &[CdnaScoreRow],
    current: &mut CdnaScoreRow,
    mut parents: Option<&mut CdnaParentRows>,
) {
    current.u5m.fill(0);
    current.u5i.fill(NEG_INF);
    current.u5d.fill(NEG_INF);
    current.cm.fill(0);
    current.ci.fill(NEG_INF);
    current.cd.fill(NEG_INF);
    current.u3m.fill(NEG_INF);
    current.u3i.fill(NEG_INF);
    current.u3d.fill(NEG_INF);
    if let Some(parents) = parents.as_deref_mut() {
        parents.fill_none();
    }
    let prior = |advance| {
        history
            .iter()
            .find(|scores| scores.row == row_index - advance)
            .expect("cDNA rolling history")
    };
    let mut w5 = IntronCandidateWindow::default();
    let mut wc = [
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
        IntronCandidateWindow::default(),
    ];
    let mut w3 = IntronCandidateWindow::default();
    for j in 0..=target.bases.len() {
        // 5-prime UTR local EST subgraph.
        if j >= intron.min_len as usize {
            let start = j - intron.min_len as usize;
            let (score, state) = best(&[
                (current.u5m[start], State::M),
                (current.u5i[start], State::I),
                (current.u5d[start], State::D),
            ]);
            if score > 0 {
                if let Some(donor) = splice_score(
                    &target.bases,
                    start,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                ) {
                    w5.insert(IntronCandidate {
                        start,
                        score: add(score, intron.open_penalty.saturating_add(donor)),
                        state_rank: state.rank(),
                    });
                }
            }
        }
        if j > intron.max_len as usize {
            w5.expire_before(j - intron.max_len as usize);
        }
        if j >= 2 {
            if let (Some(candidate), Some(acceptor)) = (
                w5.best(),
                splice_score(
                    &target.bases,
                    j - 2,
                    SpliceType::AcceptorForward,
                    intron.force_gtag,
                ),
            ) {
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::U5M,
                    j,
                    add(candidate.score, acceptor),
                    match state_from_rank(candidate.state_rank) {
                        State::M => CdnaState::U5M,
                        State::I => CdnaState::U5I,
                        State::D => CdnaState::U5D,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::intron(
                        TraceRun {
                            transition_id: 6,
                            op: Op::Splice5,
                            query_advance: 0,
                            target_advance: 2,
                            repeats: 1,
                        },
                        TraceRun {
                            transition_id: 7,
                            op: Op::Intron,
                            query_advance: 0,
                            target_advance: (j - candidate.start - 4) as u32,
                            repeats: 1,
                        },
                        TraceRun {
                            transition_id: 8,
                            op: Op::Splice3,
                            query_advance: 0,
                            target_advance: 2,
                            repeats: 1,
                        },
                    ),
                );
            }
        }
        if row_index > 0 && j > 0 {
            let previous = prior(1);
            let (score, state) = best(&[
                (previous.u5m[j - 1], State::M),
                (previous.u5i[j - 1], State::I),
                (previous.u5d[j - 1], State::D),
            ]);
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U5M,
                j,
                add(
                    score,
                    dna_score(query.bases[row_index - 1], target.bases[j - 1], scoring),
                ),
                match state {
                    State::M => CdnaState::U5M,
                    State::I => CdnaState::U5I,
                    State::D => CdnaState::U5D,
                    State::Stop => unreachable!(),
                },
                TraceFragment::one(TraceRun {
                    transition_id: 1,
                    op: Op::Match,
                    query_advance: 1,
                    target_advance: 1,
                    repeats: 1,
                }),
            );
        }
        if row_index > 0 {
            let previous = prior(1);
            let (score, state, id) = if add(previous.u5m[j], scoring.gap_open)
                >= add(previous.u5i[j], scoring.gap_extend)
            {
                (add(previous.u5m[j], scoring.gap_open), CdnaState::U5M, 2)
            } else {
                (add(previous.u5i[j], scoring.gap_extend), CdnaState::U5I, 4)
            };
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U5I,
                j,
                score,
                state,
                TraceFragment::one(TraceRun {
                    transition_id: id,
                    op: Op::Insert,
                    query_advance: 1,
                    target_advance: 0,
                    repeats: 1,
                }),
            );
        }
        if j > 0 {
            let (score, state, id) = if add(current.u5m[j - 1], scoring.gap_open)
                >= add(current.u5d[j - 1], scoring.gap_extend)
            {
                (add(current.u5m[j - 1], scoring.gap_open), CdnaState::U5M, 3)
            } else {
                (
                    add(current.u5d[j - 1], scoring.gap_extend),
                    CdnaState::U5D,
                    5,
                )
            };
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U5D,
                j,
                score,
                state,
                TraceFragment::one(TraceRun {
                    transition_id: id,
                    op: Op::Delete,
                    query_advance: 0,
                    target_advance: 1,
                    repeats: 1,
                }),
            );
        }
        let (score, state) = best(&[
            (current.u5m[j], State::M),
            (current.u5i[j], State::I),
            (current.u5d[j], State::D),
        ]);
        cdna_row_relax(
            current,
            &mut parents,
            CdnaState::CodingM,
            j,
            score,
            match state {
                State::M => CdnaState::U5M,
                State::I => CdnaState::U5I,
                State::D => CdnaState::U5D,
                State::Stop => unreachable!(),
            },
            TraceFragment::empty(),
        );

        // CDS phase-aware target introns and codon transitions.
        if row_index >= 3 {
            let previous = prior(3);
            for (phase, window) in wc.iter_mut().enumerate() {
                let post = 3 - phase;
                if j >= post + intron.min_len as usize + phase {
                    let post_start = j - post;
                    let donor = post_start - intron.min_len as usize;
                    let source_j = donor - phase;
                    let (score, state) = best(&[
                        (previous.cm[source_j], State::M),
                        (previous.ci[source_j], State::I),
                        (previous.cd[source_j], State::D),
                    ]);
                    if score > 0 {
                        if let Some(donor_score) = splice_score(
                            &target.bases,
                            donor,
                            SpliceType::DonorForward,
                            intron.force_gtag,
                        ) {
                            window.insert(IntronCandidate {
                                start: donor,
                                score: add(score, intron.open_penalty.saturating_add(donor_score)),
                                state_rank: state.rank(),
                            });
                        }
                    }
                    if post_start > intron.max_len as usize {
                        window.expire_before(post_start - intron.max_len as usize);
                    }
                    if post_start >= 2 {
                        if let (Some(candidate), Some(acceptor)) = (
                            window.best(),
                            splice_score(
                                &target.bases,
                                post_start - 2,
                                SpliceType::AcceptorForward,
                                intron.force_gtag,
                            ),
                        ) {
                            let pre = candidate.start - phase;
                            let mut codon = Vec::with_capacity(3);
                            codon.extend_from_slice(&target.bases[pre..candidate.start]);
                            codon.extend_from_slice(&target.bases[post_start..j]);
                            if codon.len() == 3 {
                                let fragment = if phase == 0 {
                                    TraceFragment::intron_match(
                                        TraceRun {
                                            transition_id: 20,
                                            op: Op::Splice5,
                                            query_advance: 0,
                                            target_advance: 2,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 21,
                                            op: Op::Intron,
                                            query_advance: 0,
                                            target_advance: (post_start - candidate.start - 4)
                                                as u32,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 22,
                                            op: Op::Splice3,
                                            query_advance: 0,
                                            target_advance: 2,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 1,
                                            op: Op::Match,
                                            query_advance: 3,
                                            target_advance: 3,
                                            repeats: 1,
                                        },
                                    )
                                } else {
                                    TraceFragment::phase_intron(
                                        TraceRun {
                                            transition_id: 23,
                                            op: Op::SplitCodon,
                                            query_advance: phase as u32,
                                            target_advance: phase as u32,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 20,
                                            op: Op::Splice5,
                                            query_advance: 0,
                                            target_advance: 2,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 21,
                                            op: Op::Intron,
                                            query_advance: 0,
                                            target_advance: (post_start - candidate.start - 4)
                                                as u32,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 22,
                                            op: Op::Splice3,
                                            query_advance: 0,
                                            target_advance: 2,
                                            repeats: 1,
                                        },
                                        TraceRun {
                                            transition_id: 24,
                                            op: Op::SplitCodon,
                                            query_advance: post as u32,
                                            target_advance: post as u32,
                                            repeats: 1,
                                        },
                                    )
                                };
                                cdna_row_relax(
                                    current,
                                    &mut parents,
                                    CdnaState::CodingM,
                                    j,
                                    add(
                                        add(candidate.score, acceptor),
                                        protein_score(
                                            translate_dna(
                                                &query.bases[row_index - 3..row_index],
                                                0,
                                            )[0],
                                            translate_dna(&codon, 0)[0],
                                            scoring,
                                        ),
                                    ),
                                    match state_from_rank(candidate.state_rank) {
                                        State::M => CdnaState::CodingM,
                                        State::I => CdnaState::CodingI,
                                        State::D => CdnaState::CodingD,
                                        State::Stop => unreachable!(),
                                    },
                                    fragment,
                                );
                            }
                        }
                    }
                }
            }
            if j >= 3 {
                let (score, state, id) = if add(current.cm[j - 3], scoring.codon_gap_open)
                    >= add(current.cd[j - 3], scoring.codon_gap_extend)
                {
                    (
                        add(current.cm[j - 3], scoring.codon_gap_open),
                        CdnaState::CodingM,
                        3,
                    )
                } else {
                    (
                        add(current.cd[j - 3], scoring.codon_gap_extend),
                        CdnaState::CodingD,
                        5,
                    )
                };
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::CodingD,
                    j,
                    score,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id: id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 3,
                        repeats: 1,
                    }),
                );
            }
            let (score, state, id) = if add(previous.cm[j], scoring.codon_gap_open)
                >= add(previous.ci[j], scoring.codon_gap_extend)
            {
                (
                    add(previous.cm[j], scoring.codon_gap_open),
                    CdnaState::CodingM,
                    2,
                )
            } else {
                (
                    add(previous.ci[j], scoring.codon_gap_extend),
                    CdnaState::CodingI,
                    4,
                )
            };
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::CodingI,
                j,
                score,
                state,
                TraceFragment::one(TraceRun {
                    transition_id: id,
                    op: Op::Insert,
                    query_advance: 3,
                    target_advance: 0,
                    repeats: 1,
                }),
            );
            if j >= 3 {
                let (score, state) = best(&[
                    (previous.cm[j - 3], State::M),
                    (previous.ci[j - 3], State::I),
                    (previous.cd[j - 3], State::D),
                ]);
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::CodingM,
                    j,
                    add(
                        score,
                        protein_score(
                            translate_dna(&query.bases[row_index - 3..row_index], 0)[0],
                            translate_dna(&target.bases[j - 3..j], 0)[0],
                            scoring,
                        ),
                    ),
                    match state {
                        State::M => CdnaState::CodingM,
                        State::I => CdnaState::CodingI,
                        State::D => CdnaState::CodingD,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::one(TraceRun {
                        transition_id: 1,
                        op: Op::Match,
                        query_advance: 3,
                        target_advance: 3,
                        repeats: 1,
                    }),
                );
            }
        }
        for advance in [1_usize, 2, 4, 5] {
            if row_index >= advance {
                let previous = prior(advance);
                let (score, state) = best(&[
                    (previous.cm[j], State::M),
                    (previous.ci[j], State::I),
                    (previous.cd[j], State::D),
                ]);
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::CodingM,
                    j,
                    add(score, scoring.frameshift),
                    match state {
                        State::M => CdnaState::CodingM,
                        State::I => CdnaState::CodingI,
                        State::D => CdnaState::CodingD,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::one(TraceRun {
                        transition_id: 6,
                        op: Op::Frameshift,
                        query_advance: advance as u32,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j >= advance {
                let (score, state) = best(&[
                    (current.cm[j - advance], State::M),
                    (current.ci[j - advance], State::I),
                    (current.cd[j - advance], State::D),
                ]);
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::CodingM,
                    j,
                    add(score, scoring.frameshift),
                    match state {
                        State::M => CdnaState::CodingM,
                        State::I => CdnaState::CodingI,
                        State::D => CdnaState::CodingD,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::one(TraceRun {
                        transition_id: 7,
                        op: Op::Frameshift,
                        query_advance: 0,
                        target_advance: advance as u32,
                        repeats: 1,
                    }),
                );
            }
        }

        // 3-prime UTR, including the CDS epsilon exit.
        if j >= intron.min_len as usize {
            let start = j - intron.min_len as usize;
            let (score, state) = best(&[
                (current.u3m[start], State::M),
                (current.u3i[start], State::I),
                (current.u3d[start], State::D),
            ]);
            if score > 0 {
                if let Some(donor) = splice_score(
                    &target.bases,
                    start,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                ) {
                    w3.insert(IntronCandidate {
                        start,
                        score: add(score, intron.open_penalty.saturating_add(donor)),
                        state_rank: state.rank(),
                    });
                }
            }
        }
        if j > intron.max_len as usize {
            w3.expire_before(j - intron.max_len as usize);
        }
        if j >= 2 {
            if let (Some(candidate), Some(acceptor)) = (
                w3.best(),
                splice_score(
                    &target.bases,
                    j - 2,
                    SpliceType::AcceptorForward,
                    intron.force_gtag,
                ),
            ) {
                cdna_row_relax(
                    current,
                    &mut parents,
                    CdnaState::U3M,
                    j,
                    add(candidate.score, acceptor),
                    match state_from_rank(candidate.state_rank) {
                        State::M => CdnaState::U3M,
                        State::I => CdnaState::U3I,
                        State::D => CdnaState::U3D,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::intron(
                        TraceRun {
                            transition_id: 6,
                            op: Op::Splice5,
                            query_advance: 0,
                            target_advance: 2,
                            repeats: 1,
                        },
                        TraceRun {
                            transition_id: 7,
                            op: Op::Intron,
                            query_advance: 0,
                            target_advance: (j - candidate.start - 4) as u32,
                            repeats: 1,
                        },
                        TraceRun {
                            transition_id: 8,
                            op: Op::Splice3,
                            query_advance: 0,
                            target_advance: 2,
                            repeats: 1,
                        },
                    ),
                );
            }
        }
        if row_index > 0 && j > 0 {
            let previous = prior(1);
            let (score, state) = best(&[
                (previous.u3m[j - 1], State::M),
                (previous.u3i[j - 1], State::I),
                (previous.u3d[j - 1], State::D),
            ]);
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U3M,
                j,
                add(
                    score,
                    dna_score(query.bases[row_index - 1], target.bases[j - 1], scoring),
                ),
                match state {
                    State::M => CdnaState::U3M,
                    State::I => CdnaState::U3I,
                    State::D => CdnaState::U3D,
                    State::Stop => unreachable!(),
                },
                TraceFragment::one(TraceRun {
                    transition_id: 1,
                    op: Op::Match,
                    query_advance: 1,
                    target_advance: 1,
                    repeats: 1,
                }),
            );
        }
        if row_index > 0 {
            let previous = prior(1);
            let (score, state, id) = if add(previous.u3m[j], scoring.gap_open)
                >= add(previous.u3i[j], scoring.gap_extend)
            {
                (add(previous.u3m[j], scoring.gap_open), CdnaState::U3M, 2)
            } else {
                (add(previous.u3i[j], scoring.gap_extend), CdnaState::U3I, 4)
            };
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U3I,
                j,
                score,
                state,
                TraceFragment::one(TraceRun {
                    transition_id: id,
                    op: Op::Insert,
                    query_advance: 1,
                    target_advance: 0,
                    repeats: 1,
                }),
            );
        }
        if j > 0 {
            let (score, state, id) = if add(current.u3m[j - 1], scoring.gap_open)
                >= add(current.u3d[j - 1], scoring.gap_extend)
            {
                (add(current.u3m[j - 1], scoring.gap_open), CdnaState::U3M, 3)
            } else {
                (
                    add(current.u3d[j - 1], scoring.gap_extend),
                    CdnaState::U3D,
                    5,
                )
            };
            cdna_row_relax(
                current,
                &mut parents,
                CdnaState::U3D,
                j,
                score,
                state,
                TraceFragment::one(TraceRun {
                    transition_id: id,
                    op: Op::Delete,
                    query_advance: 0,
                    target_advance: 1,
                    repeats: 1,
                }),
            );
        }
        let (score, state) = best(&[
            (current.cm[j], State::M),
            (current.ci[j], State::I),
            (current.cd[j], State::D),
        ]);
        cdna_row_relax(
            current,
            &mut parents,
            CdnaState::U3M,
            j,
            score,
            match state {
                State::M => CdnaState::CodingM,
                State::I => CdnaState::CodingI,
                State::D => CdnaState::CodingD,
                State::Stop => unreachable!(),
            },
            TraceFragment::empty(),
        );
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CheckpointPlan {
    pub(crate) stride: usize,
    peak_bytes: u128,
}

/// Estimate the largest allocation made while reconstructing one section.
///
/// This counts score checkpoints, the three temporary parent matrices, and
/// the six rolling score rows used by block reconstruction.  Vec headers and
/// the alignment output are deliberately excluded: they are small relative to
/// a DP row and do not scale with the matrix area.
fn checkpoint_peak_bytes_with_parent(
    query_len: usize,
    target_len: usize,
    stride: usize,
    parent_bytes: usize,
) -> u128 {
    let cols = target_len as u128 + 1;
    let stride = stride.max(1) as u128;
    let checkpoint_rows = checkpoint_count(query_len, stride as usize) as u128;
    let score_bytes = size_of::<Score>() as u128;
    let parent_bytes = parent_bytes as u128;
    let saved_scores = checkpoint_rows * 3 * cols * score_bytes;
    let reconstruction = (stride + 1) * 3 * cols * parent_bytes;
    let rolling_scores = 6 * cols * score_bytes;
    let checkpoint_headers = checkpoint_rows * size_of::<ScoreRow>() as u128;
    saved_scores + checkpoint_headers + reconstruction + rolling_scores
}

/// Pick the lowest estimated peak that fits the requested budget.
///
/// A DP row plus one reconstruction block has a non-zero minimum footprint.
/// When no candidate fits (notably `--dpmemory 0`), the minimum-footprint
/// plan is selected instead.  The caller therefore gets an exact alignment
/// and a deterministic forced-checkpoint mode; this is not an RSS hard cap.
fn checkpoint_plan_with_parent(
    query_len: usize,
    target_len: usize,
    budget_bytes: Option<u128>,
    parent_bytes: usize,
) -> CheckpointPlan {
    let mut best = CheckpointPlan {
        stride: 1,
        peak_bytes: checkpoint_peak_bytes_with_parent(query_len, target_len, 1, parent_bytes),
    };
    for stride in 2..=query_len.max(1) {
        let peak_bytes =
            checkpoint_peak_bytes_with_parent(query_len, target_len, stride, parent_bytes);
        let candidate = CheckpointPlan { stride, peak_bytes };
        let candidate_fits = budget_bytes.is_none_or(|budget| peak_bytes <= budget);
        let best_fits = budget_bytes.is_none_or(|budget| best.peak_bytes <= budget);
        if (candidate_fits && !best_fits)
            || (candidate_fits == best_fits && peak_bytes < best.peak_bytes)
        {
            best = candidate;
        }
    }
    best
}

/// As [`checkpoint_plan_with_parent`], for models whose edges can look back
/// multiple query rows and which therefore checkpoint a score-row history.
fn checkpoint_plan_with_history(
    query_len: usize,
    target_len: usize,
    budget_bytes: Option<u128>,
    parent_bytes: usize,
    parent_rows: usize,
    score_rows: usize,
    history_rows: usize,
) -> CheckpointPlan {
    checkpoint_plan_with_history_and_checkpoint_bytes(
        query_len,
        target_len,
        budget_bytes,
        parent_bytes,
        parent_rows,
        score_rows,
        history_rows,
        0,
    )
}

/// Choose a checkpoint section size when each saved section also owns a
/// variable-size semantic snapshot (for example, monotonic intron queues).
///
/// `checkpoint_bytes` is a conservative bound for one snapshot.  It is kept
/// separate from score-row history because score history scales only with the
/// number of DP states, whereas queue contents scale with the allowed intron
/// span.  Omitting it would make a requested memory budget select an
/// unsoundly small stride.
#[allow(clippy::too_many_arguments)]
fn checkpoint_plan_with_history_and_checkpoint_bytes(
    query_len: usize,
    target_len: usize,
    budget_bytes: Option<u128>,
    parent_bytes: usize,
    parent_rows: usize,
    score_rows: usize,
    history_rows: usize,
    checkpoint_bytes: u128,
) -> CheckpointPlan {
    let peak = |stride| {
        let cols = target_len as u128 + 1;
        let checkpoints = checkpoint_count(query_len, stride) as u128;
        let score_bytes = size_of::<Score>() as u128;
        let saved_scores =
            checkpoints * history_rows as u128 * score_rows as u128 * cols * score_bytes;
        let checkpoint_headers = checkpoints * size_of::<ScoreCheckpoint>() as u128;
        let reconstruction =
            (stride as u128 + 1) * parent_rows as u128 * cols * parent_bytes as u128;
        let rolling_scores = history_rows as u128 * score_rows as u128 * cols * score_bytes;
        saved_scores
            + checkpoint_headers
            + checkpoints * checkpoint_bytes
            + reconstruction
            + rolling_scores
    };
    let mut best = CheckpointPlan {
        stride: 1,
        peak_bytes: peak(1),
    };
    for stride in 2..=query_len.max(1) {
        let candidate = CheckpointPlan {
            stride,
            peak_bytes: peak(stride),
        };
        let candidate_fits = budget_bytes.is_none_or(|budget| candidate.peak_bytes <= budget);
        let best_fits = budget_bytes.is_none_or(|budget| best.peak_bytes <= budget);
        if (candidate_fits && !best_fits)
            || (candidate_fits == best_fits && candidate.peak_bytes < best.peak_bytes)
        {
            best = candidate;
        }
    }
    best
}

/// Plan exact genome-to-genome checkpoint sections.  Every boundary retains
/// nine score-state rows and a clone of all UTR/CDS query and joint queues.
#[allow(dead_code)] // used when genome2genome dispatch is wired in
pub(crate) fn genome_checkpoint_plan(
    query_len: usize,
    target_len: usize,
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> CheckpointPlan {
    let queue_snapshot_bytes = genome_checkpoint_queue_upper_bound_bytes(target_len + 1, intron)
        + size_of::<GenomeCheckpoint>() as u128;
    checkpoint_plan_with_history_and_checkpoint_bytes(
        query_len,
        target_len,
        budget_bytes,
        size_of::<Option<GenomeParent>>(),
        9,
        9,
        genome_checkpoint_history_rows(intron),
        queue_snapshot_bytes,
    )
}

fn checkpoint_plan(
    query_len: usize,
    target_len: usize,
    budget_bytes: Option<u128>,
) -> CheckpointPlan {
    checkpoint_plan_with_parent(query_len, target_len, budget_bytes, size_of::<State>())
}

fn score_checkpoints(
    query: &Sequence,
    target: &[u8],
    model: Model,
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
    stride: usize,
) -> (Vec<ScoreRow>, EndPoint, usize) {
    let (n, m) = (query.bases.len(), target.len());
    let mut previous_m = vec![NEG_INF; m + 1];
    let mut previous_i = vec![NEG_INF; m + 1];
    let mut previous_d = vec![NEG_INF; m + 1];
    let mut current_m = vec![NEG_INF; m + 1];
    let mut current_i = vec![NEG_INF; m + 1];
    let mut current_d = vec![NEG_INF; m + 1];
    let mut checkpoints = Vec::with_capacity(checkpoint_count(n, stride));
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

fn score_est_checkpoints(
    query: &Sequence,
    bases: &[u8],
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> (Vec<ScoreRow>, EndPoint, usize) {
    let (n, m) = (query.bases.len(), bases.len());
    let stride =
        checkpoint_plan_with_parent(n, m, budget_bytes, size_of::<Option<EstParent>>()).stride;
    let mut previous_m = vec![NEG_INF; m + 1];
    let mut previous_i = vec![NEG_INF; m + 1];
    let mut previous_d = vec![NEG_INF; m + 1];
    let mut current_m = vec![NEG_INF; m + 1];
    let mut current_i = vec![NEG_INF; m + 1];
    let mut current_d = vec![NEG_INF; m + 1];
    let mut checkpoints = Vec::with_capacity(checkpoint_count(n, stride));
    let mut end = EndPoint {
        i: 0,
        j: 0,
        state: State::M,
        score: 0,
    };
    for row in 0..=n {
        fill_est_row(
            query,
            bases,
            row,
            intron,
            (row > 0).then_some((&previous_m, &previous_i, &previous_d)),
            &mut current_m,
            &mut current_i,
            &mut current_d,
            None,
        );
        for j in 0..=m {
            consider_end(&mut end, row, j, State::M, current_m[j]);
            consider_end(&mut end, row, j, State::I, current_i[j]);
            consider_end(&mut end, row, j, State::D, current_d[j]);
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

fn score_p2g_checkpoints(
    query: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    bestfit: bool,
    budget_bytes: Option<u128>,
) -> (Vec<ScoreRow>, EndPoint, usize) {
    let (n, m) = (query.bases.len(), bases.len());
    let stride =
        checkpoint_plan_with_parent(n, m, budget_bytes, size_of::<Option<P2GenomeParent>>()).stride;
    let mut previous_m = vec![0; m + 1];
    let mut previous_i = vec![NEG_INF; m + 1];
    let mut previous_d = vec![NEG_INF; m + 1];
    let mut current_m = vec![NEG_INF; m + 1];
    let mut current_i = vec![NEG_INF; m + 1];
    let mut current_d = vec![NEG_INF; m + 1];
    let mut checkpoints = Vec::with_capacity(checkpoint_count(n, stride));
    checkpoints.push(ScoreRow {
        row: 0,
        m: previous_m.clone(),
        i: previous_i.clone(),
        d: previous_d.clone(),
    });
    let mut end = EndPoint {
        i: 0,
        j: 0,
        state: State::M,
        score: if bestfit { NEG_INF } else { 0 },
    };
    for row in 1..=n {
        fill_p2g_row(
            query,
            bases,
            row,
            scoring,
            intron,
            bestfit,
            (&previous_m, &previous_i, &previous_d),
            &mut current_m,
            &mut current_i,
            &mut current_d,
            None,
        );
        if !bestfit || row == n {
            for j in 0..=m {
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

fn score_c2g_checkpoints(
    query: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> (Vec<ScoreCheckpoint>, EndPoint, usize) {
    let (n, m) = (query.bases.len(), bases.len());
    let stride = checkpoint_plan_with_history(
        n,
        m,
        budget_bytes,
        size_of::<Option<P2GenomeParent>>(),
        3,
        3,
        6,
    )
    .stride;
    let mut history = vec![ScoreRow {
        row: 0,
        m: vec![0; m + 1],
        i: vec![NEG_INF; m + 1],
        d: vec![NEG_INF; m + 1],
    }];
    let mut checkpoints = Vec::with_capacity(checkpoint_count(n, stride));
    checkpoints.push(ScoreCheckpoint {
        row: 0,
        history: history.clone(),
    });
    let mut end = EndPoint {
        i: 0,
        j: 0,
        state: State::M,
        score: 0,
    };
    for row in 1..=n {
        let mut current = ScoreRow {
            row,
            m: vec![NEG_INF; m + 1],
            i: vec![NEG_INF; m + 1],
            d: vec![NEG_INF; m + 1],
        };
        fill_c2g_row(
            query,
            bases,
            row,
            scoring,
            intron,
            &history,
            &mut current.m,
            &mut current.i,
            &mut current.d,
            None,
        );
        for j in 0..=m {
            consider_end(&mut end, row, j, State::M, current.m[j]);
            consider_end(&mut end, row, j, State::I, current.i[j]);
            consider_end(&mut end, row, j, State::D, current.d[j]);
        }
        history.push(current);
        if history.len() > 6 {
            history.remove(0);
        }
        if row % stride == 0 || row == n {
            checkpoints.push(ScoreCheckpoint {
                row,
                history: history.clone(),
            });
        }
    }
    (checkpoints, end, stride)
}

fn empty_cdna_score_row(row: usize, cols: usize) -> CdnaScoreRow {
    CdnaScoreRow {
        row,
        u5m: vec![0; cols],
        u5i: vec![NEG_INF; cols],
        u5d: vec![NEG_INF; cols],
        cm: vec![0; cols],
        ci: vec![NEG_INF; cols],
        cd: vec![NEG_INF; cols],
        u3m: vec![NEG_INF; cols],
        u3i: vec![NEG_INF; cols],
        u3d: vec![NEG_INF; cols],
    }
}

fn score_cdna_checkpoints(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> (Vec<CdnaCheckpoint>, CdnaEndPoint, usize) {
    let (n, cols) = (query.bases.len(), target.bases.len() + 1);
    let stride = checkpoint_plan_with_history(
        n,
        target.bases.len(),
        budget_bytes,
        size_of::<Option<CdnaParent>>(),
        9,
        9,
        6,
    )
    .stride;
    let mut history = Vec::with_capacity(6);
    let mut first = empty_cdna_score_row(0, cols);
    fill_cdna_row(
        query, target, 0, scoring, intron, &history, &mut first, None,
    );
    history.push(first);
    let mut checkpoints = Vec::with_capacity(checkpoint_count(n, stride));
    checkpoints.push(CdnaCheckpoint {
        row: 0,
        history: history.clone(),
    });
    let mut end = CdnaEndPoint {
        i: 0,
        j: 0,
        state: CdnaState::CodingM,
        score: 0,
    };
    for j in 0..cols {
        for (score, state) in [
            (history[0].u3m[j], CdnaState::U3M),
            (history[0].u3i[j], CdnaState::U3I),
            (history[0].u3d[j], CdnaState::U3D),
            (history[0].cm[j], CdnaState::CodingM),
            (history[0].ci[j], CdnaState::CodingI),
            (history[0].cd[j], CdnaState::CodingD),
        ] {
            if (score, 0, j) > (end.score, end.i, end.j) {
                end = CdnaEndPoint {
                    i: 0,
                    j,
                    state,
                    score,
                };
            }
        }
    }
    for row in 1..=n {
        let mut current = empty_cdna_score_row(row, cols);
        fill_cdna_row(
            query,
            target,
            row,
            scoring,
            intron,
            &history,
            &mut current,
            None,
        );
        for j in 0..cols {
            for (score, state) in [
                (current.u3m[j], CdnaState::U3M),
                (current.u3i[j], CdnaState::U3I),
                (current.u3d[j], CdnaState::U3D),
                (current.cm[j], CdnaState::CodingM),
                (current.ci[j], CdnaState::CodingI),
                (current.cd[j], CdnaState::CodingD),
            ] {
                if score > end.score {
                    end = CdnaEndPoint {
                        i: row,
                        j,
                        state,
                        score,
                    };
                }
            }
        }
        history.push(current);
        if history.len() > 6 {
            history.remove(0);
        }
        if row % stride == 0 || row == n {
            checkpoints.push(CdnaCheckpoint {
                row,
                history: history.clone(),
            });
        }
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

fn rebuild_est_block(
    query: &Sequence,
    bases: &[u8],
    intron: IntronScoring,
    checkpoint: &ScoreRow,
    end_row: usize,
) -> EstParentRows {
    let cols = bases.len() + 1;
    let rows = end_row - checkpoint.row;
    let mut pm = vec![None; (rows + 1) * cols];
    let mut pi = vec![None; (rows + 1) * cols];
    let mut pd = vec![None; (rows + 1) * cols];
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
        fill_est_row(
            query,
            bases,
            row,
            intron,
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

fn rebuild_p2g_block(
    query: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    bestfit: bool,
    checkpoint: &ScoreRow,
    end_row: usize,
) -> P2ParentRows {
    let cols = bases.len() + 1;
    let rows = end_row - checkpoint.row;
    let mut pm = vec![None; (rows + 1) * cols];
    let mut pi = vec![None; (rows + 1) * cols];
    let mut pd = vec![None; (rows + 1) * cols];
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
        fill_p2g_row(
            query,
            bases,
            row,
            scoring,
            intron,
            bestfit,
            (&previous_m, &previous_i, &previous_d),
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

fn rebuild_c2g_block(
    query: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    checkpoint: &ScoreCheckpoint,
    end_row: usize,
) -> P2ParentRows {
    let cols = bases.len() + 1;
    let rows = end_row - checkpoint.row;
    let mut pm = vec![None; (rows + 1) * cols];
    let mut pi = vec![None; (rows + 1) * cols];
    let mut pd = vec![None; (rows + 1) * cols];
    let mut history = checkpoint.history.clone();
    for row in checkpoint.row + 1..=end_row {
        let offset = (row - checkpoint.row) * cols;
        let mut current = ScoreRow {
            row,
            m: vec![NEG_INF; cols],
            i: vec![NEG_INF; cols],
            d: vec![NEG_INF; cols],
        };
        fill_c2g_row(
            query,
            bases,
            row,
            scoring,
            intron,
            &history,
            &mut current.m,
            &mut current.i,
            &mut current.d,
            Some((
                &mut pm[offset..offset + cols],
                &mut pi[offset..offset + cols],
                &mut pd[offset..offset + cols],
            )),
        );
        history.push(current);
        if history.len() > 6 {
            history.remove(0);
        }
    }
    (pm, pi, pd)
}

fn rebuild_cdna_block(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    checkpoint: &CdnaCheckpoint,
    end_row: usize,
) -> CdnaParentRows {
    let cols = target.bases.len() + 1;
    let rows = end_row - checkpoint.row;
    let mut parents = CdnaParentRows::new((rows + 1) * cols);
    let mut history = checkpoint.history.clone();
    for row in checkpoint.row + 1..=end_row {
        let mut current = empty_cdna_score_row(row, cols);
        let offset = (row - checkpoint.row) * cols;
        let mut row_parents = CdnaParentRows {
            u5m: parents.u5m[offset..offset + cols].to_vec(),
            u5i: parents.u5i[offset..offset + cols].to_vec(),
            u5d: parents.u5d[offset..offset + cols].to_vec(),
            cm: parents.cm[offset..offset + cols].to_vec(),
            ci: parents.ci[offset..offset + cols].to_vec(),
            cd: parents.cd[offset..offset + cols].to_vec(),
            u3m: parents.u3m[offset..offset + cols].to_vec(),
            u3i: parents.u3i[offset..offset + cols].to_vec(),
            u3d: parents.u3d[offset..offset + cols].to_vec(),
        };
        fill_cdna_row(
            query,
            target,
            row,
            scoring,
            intron,
            &history,
            &mut current,
            Some(&mut row_parents),
        );
        parents.u5m[offset..offset + cols].copy_from_slice(&row_parents.u5m);
        parents.u5i[offset..offset + cols].copy_from_slice(&row_parents.u5i);
        parents.u5d[offset..offset + cols].copy_from_slice(&row_parents.u5d);
        parents.cm[offset..offset + cols].copy_from_slice(&row_parents.cm);
        parents.ci[offset..offset + cols].copy_from_slice(&row_parents.ci);
        parents.cd[offset..offset + cols].copy_from_slice(&row_parents.cd);
        parents.u3m[offset..offset + cols].copy_from_slice(&row_parents.u3m);
        parents.u3i[offset..offset + cols].copy_from_slice(&row_parents.u3i);
        parents.u3d[offset..offset + cols].copy_from_slice(&row_parents.u3d);
        history.push(current);
        if history.len() > 6 {
            history.remove(0);
        }
    }
    parents
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
    memory_budget: Option<u128>,
) -> Alignment {
    if model == Model::Ungapped {
        return align_with_scorer(query, target, model, scoring, strand, scorer);
    }
    let reverse_target = (strand == Strand::Reverse).then(|| reverse_complement(&target.bases));
    let oriented_target = reverse_target.as_deref().unwrap_or(&target.bases);
    let plan = checkpoint_plan(query.bases.len(), oriented_target.len(), memory_budget);
    let (checkpoints, end, stride) =
        score_checkpoints(query, oriented_target, model, scoring, scorer, plan.stride);
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
            oriented_target,
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
        oriented_target,
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
    align_checkpointed_with_scorer(query, target, model, scoring, strand, dna_score, None)
}

/// Exact checkpointed EST-to-genome alignment with bounded target introns.
///
/// Score rows and traceback parents are reconstructed one query section at a
/// time, avoiding the three full score matrices and three full parent grids
/// used by [`align_est2genome_stranded`].
pub fn align_est2genome_low_memory(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    both_strands: bool,
) -> Alignment {
    align_est2genome_low_memory_with_budget(query, target, intron, both_strands, None)
}

fn align_est2genome_low_memory_with_budget(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    both_strands: bool,
    budget_bytes: Option<u128>,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let forward = align_est2genome_checkpointed_one(
        query,
        target,
        &target.bases,
        intron,
        Strand::Forward,
        budget_bytes,
    );
    if !both_strands {
        return forward;
    }
    let reverse_bases = reverse_complement(&target.bases);
    let reverse = align_est2genome_checkpointed_one(
        query,
        target,
        &reverse_bases,
        intron,
        Strand::Reverse,
        budget_bytes,
    );
    if reverse.score > forward.score {
        reverse
    } else {
        forward
    }
}

fn align_est2genome_checkpointed_one(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    intron: IntronScoring,
    strand: Strand,
    budget_bytes: Option<u128>,
) -> Alignment {
    let (checkpoints, end, stride) = score_est_checkpoints(query, bases, intron, budget_bytes);
    let cols = bases.len() + 1;
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let final_state = state;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop && i > 0 {
        let block_start = ((i - 1) / stride) * stride;
        let checkpoint = checkpoints
            .iter()
            .find(|checkpoint| checkpoint.row == block_start)
            .expect("EST checkpoint row");
        let (pm, pi, pd) = rebuild_est_block(query, bases, intron, checkpoint, i);
        loop {
            if state == State::Stop || i == block_start {
                break;
            }
            let offset = (i - block_start) * cols + j;
            let parent = match state {
                State::M => pm[offset],
                State::I => pi[offset],
                State::D => pd[offset],
                State::Stop => None,
            };
            let Some(parent) = parent else {
                state = State::Stop;
                break;
            };
            let (query_advance, target_advance) = parent.fragment.advances();
            i -= query_advance;
            j -= target_advance;
            fragments.push((parent.fragment, parent.raw_fragment));
            state = parent.state;
        }
    }
    fragments.reverse();
    let mut trace = Vec::new();
    let mut raw_trace = vec![RawStep {
        transition_id: 0,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for (fragment, raw_fragment) in fragments {
        fragment.append_to(&mut trace);
        raw_fragment.append_to(&mut raw_trace);
    }
    if final_state == State::I {
        raw_trace.push(RawStep {
            transition_id: 9,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    if final_state == State::D {
        raw_trace.push(RawStep {
            transition_id: 10,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    raw_trace.push(RawStep {
        transition_id: 11,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    let (target_start, target_end) = if strand == Strand::Forward {
        (j as u64, oriented_target_end)
    } else {
        (
            target.bases.len() as u64 - j as u64,
            target.bases.len() as u64 - oriented_target_end,
        )
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end,
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

/// Exact checkpointed local or best-fit protein-to-genome alignment.
///
/// This recomputes parent pointers in bounded query-row sections, retaining
/// the phase-aware intron semantics and raw trace of the full matrix backend.
pub fn align_protein_to_genome_low_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    bestfit: bool,
) -> Vec<Alignment> {
    align_protein_to_genome_low_memory_with_budget(
        query,
        target,
        scoring,
        intron,
        both_strands,
        bestfit,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn align_protein_to_genome_low_memory_with_budget(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    bestfit: bool,
    budget_bytes: Option<u128>,
) -> Vec<Alignment> {
    let intron = canonical_intron_scoring(intron);
    let mut out = vec![align_p2g_checkpointed_one(
        query,
        target,
        &target.bases,
        scoring,
        intron,
        Strand::Forward,
        bestfit,
        budget_bytes,
    )];
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.push(align_p2g_checkpointed_one(
            query,
            target,
            &reverse,
            scoring,
            intron,
            Strand::Reverse,
            bestfit,
            budget_bytes,
        ));
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn align_p2g_checkpointed_one(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    bestfit: bool,
    budget_bytes: Option<u128>,
) -> Alignment {
    let (checkpoints, end, stride) =
        score_p2g_checkpoints(query, bases, scoring, intron, bestfit, budget_bytes);
    let cols = bases.len() + 1;
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let final_state = state;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop && i > 0 {
        let block_start = ((i - 1) / stride) * stride;
        let checkpoint = checkpoints
            .iter()
            .find(|checkpoint| checkpoint.row == block_start)
            .expect("protein2genome checkpoint row");
        let (pm, pi, pd) = rebuild_p2g_block(query, bases, scoring, intron, bestfit, checkpoint, i);
        loop {
            if state == State::Stop || i == block_start {
                break;
            }
            let offset = (i - block_start) * cols + j;
            let parent = match state {
                State::M => pm[offset],
                State::I => pi[offset],
                State::D => pd[offset],
                State::Stop => None,
            };
            let Some(parent) = parent else {
                state = State::Stop;
                break;
            };
            let raw_fragment =
                protein2genome_raw_fragment(parent, state, i, j, query, bases, scoring, intron);
            let (query_advance, target_advance) = parent.fragment.advances();
            debug_assert!(i >= query_advance && j >= target_advance);
            i -= query_advance;
            j -= target_advance;
            fragments.push((parent.fragment, raw_fragment));
            state = parent.state;
        }
    }
    fragments.reverse();
    let mut trace = Vec::new();
    let mut raw_trace = vec![RawStep {
        transition_id: 0,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for (fragment, raw_fragment) in fragments {
        fragment.append_to(&mut trace);
        raw_trace.extend(raw_fragment);
    }
    if final_state == State::I {
        raw_trace.push(RawStep {
            transition_id: 6,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    if final_state == State::D {
        raw_trace.push(RawStep {
            transition_id: 7,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    raw_trace.push(RawStep {
        transition_id: 12,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    let (target_start, target_end) = match strand {
        Strand::Forward => (j as u64, oriented_target_end),
        Strand::Reverse => (
            (target.bases.len() as u64).saturating_sub(j as u64),
            (target.bases.len() as u64).saturating_sub(oriented_target_end),
        ),
        Strand::Unknown => unreachable!("protein2genome uses a concrete target strand"),
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end,
        query_strand: Strand::Unknown,
        target_start,
        target_end,
        target_len: target.bases.len() as u64,
        target_strand: strand,
        score: end.score,
        raw_trace,
        trace,
    }
}

pub fn align_coding_to_genome_low_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    align_coding_to_genome_low_memory_with_budget(
        query,
        target,
        scoring,
        intron,
        both_strands,
        None,
    )
}

fn align_coding_to_genome_low_memory_with_budget(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    budget_bytes: Option<u128>,
) -> Vec<Alignment> {
    let intron = canonical_intron_scoring(intron);
    let mut out = vec![align_c2g_checkpointed_one(
        query,
        target,
        &target.bases,
        scoring,
        intron,
        Strand::Forward,
        budget_bytes,
    )];
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.push(align_c2g_checkpointed_one(
            query,
            target,
            &reverse,
            scoring,
            intron,
            Strand::Reverse,
            budget_bytes,
        ));

        let reverse_query = Sequence {
            id: query.id.clone(),
            bases: reverse_complement(&query.bases),
        };
        for (bases, strand) in [
            (&target.bases, Strand::Forward),
            (&reverse, Strand::Reverse),
        ] {
            let mut alignment = align_c2g_checkpointed_one(
                &reverse_query,
                target,
                bases,
                scoring,
                intron,
                strand,
                budget_bytes,
            );
            alignment.query_start = query.bases.len() as u64 - alignment.query_start;
            alignment.query_end = query.bases.len() as u64 - alignment.query_end;
            alignment.query_strand = Strand::Reverse;
            out.push(alignment);
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn align_c2g_checkpointed_one(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    budget_bytes: Option<u128>,
) -> Alignment {
    let (checkpoints, end, stride) =
        score_c2g_checkpoints(query, bases, scoring, intron, budget_bytes);
    let cols = bases.len() + 1;
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let final_state = state;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop && i > 0 {
        let block_start = ((i - 1) / stride) * stride;
        let checkpoint = checkpoints
            .iter()
            .find(|checkpoint| checkpoint.row == block_start)
            .expect("coding2genome checkpoint row");
        let (pm, pi, pd) = rebuild_c2g_block(query, bases, scoring, intron, checkpoint, i);
        loop {
            if state == State::Stop || i <= block_start {
                break;
            }
            let offset = (i - block_start) * cols + j;
            let parent = match state {
                State::M => pm[offset],
                State::I => pi[offset],
                State::D => pd[offset],
                State::Stop => None,
            };
            let Some(parent) = parent else {
                state = State::Stop;
                break;
            };
            let raw_fragment =
                coding2genome_raw_fragment(parent, state, i, j, query, bases, scoring, intron);
            let (query_advance, target_advance) = parent.fragment.advances();
            debug_assert!(i >= query_advance && j >= target_advance);
            i -= query_advance;
            j -= target_advance;
            fragments.push((parent.fragment, raw_fragment));
            state = parent.state;
        }
    }
    fragments.reverse();
    let mut trace = Vec::new();
    let mut raw_trace = vec![RawStep {
        transition_id: 0,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for (fragment, raw_fragment) in fragments {
        fragment.append_to(&mut trace);
        raw_trace.extend(raw_fragment);
    }
    if final_state == State::I {
        raw_trace.push(RawStep {
            transition_id: 4,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    if final_state == State::D {
        raw_trace.push(RawStep {
            transition_id: 7,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    raw_trace.push(RawStep {
        transition_id: 16,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    let (target_start, target_end) = match strand {
        Strand::Forward => (j as u64, oriented_target_end),
        Strand::Reverse => (
            (target.bases.len() as u64).saturating_sub(j as u64),
            (target.bases.len() as u64).saturating_sub(oriented_target_end),
        ),
        Strand::Unknown => unreachable!("coding2genome uses a concrete target strand"),
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end,
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

pub fn align_cdna_to_genome_low_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Alignment {
    align_cdna_to_genome_low_memory_with_budget(query, target, scoring, intron, None)
}

fn align_cdna_to_genome_low_memory_with_budget(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let (checkpoints, end, stride) =
        score_cdna_checkpoints(query, target, scoring, intron, budget_bytes);
    let cols = target.bases.len() + 1;
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let final_state = state;
    let (query_end, target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    'traceback: while i > 0 {
        let block_start = ((i - 1) / stride) * stride;
        let checkpoint = checkpoints
            .iter()
            .find(|checkpoint| checkpoint.row == block_start)
            .expect("cDNA checkpoint row");
        let parents = rebuild_cdna_block(query, target, scoring, intron, checkpoint, i);
        loop {
            if i <= block_start {
                break;
            }
            let offset = (i - block_start) * cols + j;
            let Some(parent) = parents.get(state, offset) else {
                break 'traceback;
            };
            let raw_fragment =
                cdna_raw_fragment(parent, state, i, j, query, target, scoring, intron);
            let (query_advance, target_advance) = parent.fragment.advances();
            debug_assert!(i >= query_advance && j >= target_advance);
            i -= query_advance;
            j -= target_advance;
            fragments.push((parent.fragment, raw_fragment));
            state = parent.state;
        }
        if i == 0 {
            break;
        }
    }
    // The composed graph has zero-query-advance transitions at row zero
    // (notably CodingM -> U3M).  Rebuild that row once to retain them.
    if i == 0 {
        let mut row_zero = empty_cdna_score_row(0, cols);
        let mut zero_parents = CdnaParentRows::new(cols);
        fill_cdna_row(
            query,
            target,
            0,
            scoring,
            intron,
            &[],
            &mut row_zero,
            Some(&mut zero_parents),
        );
        while let Some(parent) = zero_parents.get(state, j) {
            let raw_fragment =
                cdna_raw_fragment(parent, state, i, j, query, target, scoring, intron);
            let (query_advance, target_advance) = parent.fragment.advances();
            if query_advance > i || target_advance > j {
                break;
            }
            i -= query_advance;
            j -= target_advance;
            fragments.push((parent.fragment, raw_fragment));
            state = parent.state;
        }
    }
    fragments.reverse();
    let mut trace = Vec::new();
    let mut raw_trace = vec![RawStep {
        transition_id: 300,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for (fragment, raw_fragment) in fragments {
        for atom in fragment.atoms.into_iter().flatten() {
            if atom.op == Op::Intron && atom.query_advance > 0 && atom.target_advance > 0 {
                TraceFragment::one(TraceRun {
                    query_advance: 0,
                    ..atom
                })
                .append_to(&mut trace);
                TraceFragment::one(TraceRun {
                    target_advance: 0,
                    ..atom
                })
                .append_to(&mut trace);
            } else {
                TraceFragment::one(atom).append_to(&mut trace);
            }
        }
        raw_trace.extend(raw_fragment);
    }
    let terminal = match final_state {
        CdnaState::U5I => Some(304),
        CdnaState::U5D => Some(307),
        CdnaState::CodingI => Some(4),
        CdnaState::CodingD => Some(7),
        CdnaState::U3I => Some(404),
        CdnaState::U3D => Some(407),
        _ => None,
    };
    if let Some(transition_id) = terminal {
        raw_trace.push(RawStep {
            transition_id,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    let end_id = match final_state {
        CdnaState::U5M | CdnaState::U5I | CdnaState::U5D => 308,
        CdnaState::CodingM | CdnaState::CodingI | CdnaState::CodingD => 16,
        CdnaState::U3M | CdnaState::U3I | CdnaState::U3D => 408,
    };
    raw_trace.push(RawStep {
        transition_id: end_id,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end,
        query_strand: Strand::Forward,
        target_start: j as u64,
        target_end,
        target_len: target.bases.len() as u64,
        target_strand: Strand::Forward,
        score: end.score,
        raw_trace,
        trace,
    }
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
        None,
    );
    alignment.query_strand = Strand::Unknown;
    alignment.target_strand = Strand::Unknown;
    alignment
}

fn full_affine_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * (3 * size_of::<Score>() as u128 + 3 * size_of::<State>() as u128)
}

fn full_est2genome_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * (3 * size_of::<Score>() as u128 + 3 * size_of::<Option<EstParent>>() as u128)
}

fn full_p2g_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * (3 * size_of::<Score>() as u128 + 4 * size_of::<Option<P2GenomeParent>>() as u128)
}

fn full_cdna_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * 9 * (size_of::<Score>() as u128 + size_of::<Option<CdnaParent>>() as u128)
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
        align_checkpointed_with_scorer(
            query,
            target,
            model,
            scoring,
            strand,
            dna_score,
            Some((memory_mb as u128) << 20),
        )
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
        let mut alignment = align_checkpointed_with_scorer(
            query,
            target,
            model,
            scoring,
            Strand::Forward,
            protein_score,
            Some((memory_mb as u128) << 20),
        );
        alignment.query_strand = Strand::Unknown;
        alignment.target_strand = Strand::Unknown;
        alignment
    } else {
        align_protein(query, target, model, scoring)
    }
}

pub fn align_est2genome_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    both_strands: bool,
    memory_mb: usize,
) -> Alignment {
    if full_est2genome_dp_bytes(query.bases.len(), target.bases.len()) > ((memory_mb as u128) << 20)
    {
        align_est2genome_low_memory_with_budget(
            query,
            target,
            intron,
            both_strands,
            Some((memory_mb as u128) << 20),
        )
    } else {
        align_est2genome_stranded(query, target, intron, both_strands)
    }
}

pub fn align_cdna_to_genome_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    memory_mb: usize,
) -> Alignment {
    if full_cdna_dp_bytes(query.bases.len(), target.bases.len()) > ((memory_mb as u128) << 20) {
        align_cdna_to_genome_low_memory_with_budget(
            query,
            target,
            scoring,
            intron,
            Some((memory_mb as u128) << 20),
        )
    } else {
        align_cdna_to_genome(query, target, scoring, intron)
    }
}

pub fn align_coding_to_genome_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    if full_p2g_dp_bytes(query.bases.len(), target.bases.len()) > ((memory_mb as u128) << 20) {
        align_coding_to_genome_low_memory_with_budget(
            query,
            target,
            scoring,
            intron,
            both_strands,
            Some((memory_mb as u128) << 20),
        )
    } else {
        align_coding_to_genome(query, target, scoring, intron, both_strands)
    }
}

pub fn align_protein_to_genome_with_dp_memory(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    bestfit: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    if full_p2g_dp_bytes(query.bases.len(), target.bases.len()) > ((memory_mb as u128) << 20) {
        align_protein_to_genome_low_memory_with_budget(
            query,
            target,
            scoring,
            intron,
            both_strands,
            bestfit,
            Some((memory_mb as u128) << 20),
        )
    } else if bestfit {
        align_protein_to_genome_bestfit(query, target, scoring, intron, both_strands)
    } else {
        align_protein_to_genome(query, target, scoring, intron, both_strands)
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

pub fn align_est2genome_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    intron: IntronScoring,
    both_strands: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets.iter().map(move |target| {
                align_est2genome_with_dp_memory(query, target, intron, both_strands, memory_mb)
            })
        })
        .collect()
}

pub fn align_cdna_to_genome_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    memory_mb: usize,
) -> Vec<Alignment> {
    align_cdna_to_genome_database_with_dp_memory_stranded(
        queries, targets, scoring, intron, false, memory_mb,
    )
}

/// cDNA-to-genome database alignment over the upstream four orientation
/// combinations when `both_strands` is enabled.
pub fn align_cdna_to_genome_database_with_dp_memory_stranded(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_cdna_to_genome_with_dp_memory(
                query, target, scoring, intron, memory_mb,
            ));
            if !both_strands {
                continue;
            }

            let reverse_target = Sequence {
                id: target.id.clone(),
                bases: reverse_complement(&target.bases),
            };
            let reverse_query = Sequence {
                id: query.id.clone(),
                bases: reverse_complement(&query.bases),
            };
            for (oriented_query, query_reverse) in [(query, false), (&reverse_query, true)] {
                for (oriented_target, target_reverse) in [(target, false), (&reverse_target, true)]
                {
                    if !query_reverse && !target_reverse {
                        continue;
                    }
                    let mut alignment = align_cdna_to_genome_with_dp_memory(
                        oriented_query,
                        oriented_target,
                        scoring,
                        intron,
                        memory_mb,
                    );
                    alignment.target_id = target.id.clone();
                    if target_reverse {
                        alignment.target_start = target.bases.len() as u64 - alignment.target_start;
                        alignment.target_end = target.bases.len() as u64 - alignment.target_end;
                        alignment.target_strand = Strand::Reverse;
                    }
                    if query_reverse {
                        alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                        alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                        alignment.query_strand = Strand::Reverse;
                    }
                    out.push(alignment);
                }
            }
        }
    }
    out
}

pub fn align_protein_to_genome_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    bestfit: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets.iter().flat_map(move |target| {
                align_protein_to_genome_with_dp_memory(
                    query,
                    target,
                    scoring,
                    intron,
                    both_strands,
                    bestfit,
                    memory_mb,
                )
            })
        })
        .collect()
}

pub fn align_coding_to_genome_database_with_dp_memory(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    memory_mb: usize,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets.iter().flat_map(move |target| {
                align_coding_to_genome_with_dp_memory(
                    query,
                    target,
                    scoring,
                    intron,
                    both_strands,
                    memory_mb,
                )
            })
        })
        .collect()
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
    fn checkpoint_plan_minimizes_estimated_peak_and_honors_a_feasible_budget() {
        let (query_len, target_len) = (1_024, 10_000);
        let unlimited = checkpoint_plan(query_len, target_len, None);
        for stride in 1..=query_len {
            assert!(
                unlimited.peak_bytes
                    <= checkpoint_peak_bytes_with_parent(
                        query_len,
                        target_len,
                        stride,
                        size_of::<State>(),
                    )
            );
        }

        let budget = unlimited.peak_bytes;
        let limited = checkpoint_plan(query_len, target_len, Some(budget));
        assert!(limited.peak_bytes <= budget);
        assert_eq!(limited, unlimited);

        // A physical zero-byte cap cannot contain even one DP row.  It has
        // defined forced-checkpoint semantics: select the minimum footprint.
        assert_eq!(checkpoint_plan(query_len, target_len, Some(0)), unlimited,);

        let est = checkpoint_plan_with_parent(
            query_len,
            target_len,
            None,
            size_of::<Option<EstParent>>(),
        );
        for stride in 1..=query_len {
            assert!(
                est.peak_bytes
                    <= checkpoint_peak_bytes_with_parent(
                        query_len,
                        target_len,
                        stride,
                        size_of::<Option<EstParent>>(),
                    )
            );
        }
        assert_eq!(
            checkpoint_plan_with_history(
                query_len,
                target_len,
                Some(0),
                size_of::<Option<P2GenomeParent>>(),
                3,
                3,
                6,
            ),
            checkpoint_plan_with_history(
                query_len,
                target_len,
                None,
                size_of::<Option<P2GenomeParent>>(),
                3,
                3,
                6,
            ),
        );

        // With a one-row query there is only one possible stride, so the
        // snapshot contribution can be checked exactly.
        let base = checkpoint_plan_with_history_and_checkpoint_bytes(
            1,
            9,
            None,
            size_of::<Option<P2GenomeParent>>(),
            6,
            6,
            6,
            0,
        );
        let queue_snapshot_bytes = 123_u128;
        let queued = checkpoint_plan_with_history_and_checkpoint_bytes(
            1,
            9,
            None,
            size_of::<Option<P2GenomeParent>>(),
            6,
            6,
            6,
            queue_snapshot_bytes,
        );
        assert_eq!(base.stride, 1);
        assert_eq!(queued.stride, 1);
        assert_eq!(
            queued.peak_bytes - base.peak_bytes,
            checkpoint_count(1, 1) as u128 * queue_snapshot_bytes
        );
    }

    #[test]
    fn genome_checkpoint_plan_budgets_all_state_rows_and_queue_snapshots() {
        let intron = IntronScoring {
            min_len: 30,
            max_len: 80,
            ..IntronScoring::default()
        };
        let (query_len, target_len) = (120, 90);
        let plan = genome_checkpoint_plan(query_len, target_len, intron, None);
        assert!((1..=query_len).contains(&plan.stride));
        let forced = genome_checkpoint_plan(query_len, target_len, intron, Some(0));
        assert_eq!(forced, plan);
        assert!(
            plan.peak_bytes >= genome_checkpoint_queue_upper_bound_bytes(target_len + 1, intron,)
        );
    }

    #[test]
    fn genome_checkpoint_rows_preserve_local_entry_and_all_parent_lanes() {
        let row = empty_genome_score_row(7, 3);
        assert_eq!(row.row, 7);
        assert_eq!(row.u5m, vec![0; 3]);
        assert!(
            row.u5i
                .iter()
                .chain(&row.cm)
                .chain(&row.u3m)
                .all(|score| *score == NEG_INF)
        );

        let mut parents = GenomeParentRows::new(2);
        parents.u3d[1] = Some(GenomeParent {
            state: GenomeState::CdsD,
            fragment: CompactGenomeFragment::from_trace_fragment(TraceFragment::empty()),
        });
        assert!(parents.get(GenomeState::U3D, 1).is_some());
        assert!(parents.get(GenomeState::UtrM, 1).is_none());
    }

    #[test]
    fn checkpoint_storage_is_preallocated_to_the_planned_row_count() {
        let query = sequence("q", &"ACGT".repeat(20));
        let target = sequence("t", &"TGCA".repeat(10));
        let stride = 7;
        let (checkpoints, _, actual_stride) = score_checkpoints(
            &query,
            &target.bases,
            Model::Local,
            Scoring::default(),
            dna_score,
            stride,
        );
        assert_eq!(actual_stride, stride);
        assert_eq!(
            checkpoints.len(),
            checkpoint_count(query.bases.len(), stride)
        );
        assert_eq!(checkpoints.capacity(), checkpoints.len());
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

    #[test]
    fn checkpointed_est2genome_traceback_is_identical() {
        let query = sequence("q", "ACGTACGTACGTACGT");
        let target = sequence("t", "ACGTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGACGTACGTACGT");
        let intron = IntronScoring {
            min_len: 30,
            max_len: 200,
            open_penalty: -1,
            force_gtag: true,
        };
        for both_strands in [false, true] {
            let full = align_est2genome_stranded(&query, &target, intron, both_strands);
            let low = align_est2genome_low_memory(&query, &target, intron, both_strands);
            assert_same(&low, &full);
        }
    }

    #[test]
    fn zero_memory_limit_forces_checkpointed_est2genome() {
        let query = sequence("q", "ACGTACGTACGTACGT");
        let target = sequence("t", "ACGTGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGACGTACGTACGT");
        let intron = IntronScoring {
            min_len: 30,
            max_len: 200,
            open_penalty: -1,
            force_gtag: true,
        };
        assert_same(
            &align_est2genome_with_dp_memory(&query, &target, intron, false, 0),
            &align_est2genome_stranded(&query, &target, intron, false),
        );
    }

    #[test]
    fn checkpointed_protein2genome_traceback_is_identical() {
        let query = sequence("protein", "MADQLTEQIAEFKEAFSLFDKDGDGTITT");
        let target = sequence(
            "genome",
            "ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGGAGGCCTTCTCCCTCTTTGACAAGGATGGAGATGGCACTATTACCACC",
        );
        for bestfit in [false, true] {
            for both_strands in [false, true] {
                let full = if bestfit {
                    align_protein_to_genome_bestfit(
                        &query,
                        &target,
                        Scoring::default(),
                        IntronScoring::default(),
                        both_strands,
                    )
                } else {
                    align_protein_to_genome(
                        &query,
                        &target,
                        Scoring::default(),
                        IntronScoring::default(),
                        both_strands,
                    )
                };
                let low = align_protein_to_genome_low_memory(
                    &query,
                    &target,
                    Scoring::default(),
                    IntronScoring::default(),
                    both_strands,
                    bestfit,
                );
                assert_eq!(low.len(), full.len());
                for (low, full) in low.iter().zip(&full) {
                    assert_same(low, full);
                }
            }
        }
    }

    #[test]
    fn zero_memory_limit_forces_checkpointed_protein2genome() {
        let query = sequence("protein", "MADQLTEQIAEFKEAFSLFDKDGDGTITT");
        let target = sequence(
            "genome",
            "ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGGAGGCCTTCTCCCTCTTTGACAAGGATGGAGATGGCACTATTACCACC",
        );
        let full = align_protein_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
        );
        let low = align_protein_to_genome_with_dp_memory(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
            false,
            0,
        );
        assert_eq!(low.len(), full.len());
        for (low, full) in low.iter().zip(&full) {
            assert_same(low, full);
        }
    }

    #[test]
    fn checkpointed_coding2genome_traceback_is_identical() {
        let query = sequence(
            "coding",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let target = sequence(
            "genome",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        for both_strands in [false, true] {
            let full = align_coding_to_genome(
                &query,
                &target,
                Scoring::default(),
                IntronScoring {
                    min_len: 30,
                    max_len: 200,
                    open_penalty: -1,
                    force_gtag: false,
                },
                both_strands,
            );
            let low = align_coding_to_genome_low_memory(
                &query,
                &target,
                Scoring::default(),
                IntronScoring {
                    min_len: 30,
                    max_len: 200,
                    open_penalty: -1,
                    force_gtag: false,
                },
                both_strands,
            );
            assert_eq!(low.len(), full.len());
            for (low, full) in low.iter().zip(&full) {
                assert_same(low, full);
            }
        }
    }

    #[test]
    fn checkpointed_coding2genome_matches_full_on_frameshift_boundaries() {
        let queries = ["A", "AT", "ATG", "ATGC", "ATGCA", "ATGCATG"];
        let targets = ["A", "AT", "ATG", "ATGC", "ATGCA", "CATGCAT"];
        let intron = IntronScoring {
            min_len: 4,
            max_len: 16,
            open_penalty: -1,
            force_gtag: false,
        };
        for query_bases in queries {
            for target_bases in targets {
                let query = sequence("q", query_bases);
                let target = sequence("t", target_bases);
                for both_strands in [false, true] {
                    let full = align_coding_to_genome(
                        &query,
                        &target,
                        Scoring::default(),
                        intron,
                        both_strands,
                    );
                    let low = align_coding_to_genome_low_memory(
                        &query,
                        &target,
                        Scoring::default(),
                        intron,
                        both_strands,
                    );
                    assert_eq!(low.len(), full.len());
                    for (low, full) in low.iter().zip(&full) {
                        assert_same(low, full);
                    }
                }
            }
        }
    }

    #[test]
    fn checkpointed_cdna2genome_traceback_is_identical() {
        let query = sequence(
            "coding",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let target = sequence(
            "genome",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let intron = IntronScoring {
            min_len: 30,
            max_len: 200,
            open_penalty: -1,
            force_gtag: false,
        };
        assert_same(
            &align_cdna_to_genome_low_memory(&query, &target, Scoring::default(), intron),
            &align_cdna_to_genome(&query, &target, Scoring::default(), intron),
        );
    }

    #[test]
    fn zero_memory_limit_forces_checkpointed_cdna2genome() {
        let query = sequence(
            "coding",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let target = sequence(
            "genome",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let intron = IntronScoring {
            min_len: 30,
            max_len: 200,
            open_penalty: -1,
            force_gtag: false,
        };
        assert_same(
            &align_cdna_to_genome_with_dp_memory(&query, &target, Scoring::default(), intron, 0),
            &align_cdna_to_genome(&query, &target, Scoring::default(), intron),
        );
    }
}
