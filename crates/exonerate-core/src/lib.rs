//! The stable, scalar reference backend for the first Exonerate-RS milestone.
//! It deliberately exposes a trace independent of any textual report format.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt;
use std::mem::size_of;

mod low_memory;
use low_memory::{GenomeCheckpoint, GenomeScoreRow, genome_checkpoint_plan};
pub use low_memory::{
    align_cdna_to_genome_database_with_dp_memory,
    align_cdna_to_genome_database_with_dp_memory_stranded, align_cdna_to_genome_low_memory,
    align_cdna_to_genome_with_dp_memory, align_coding_to_genome_database_with_dp_memory,
    align_coding_to_genome_low_memory, align_coding_to_genome_with_dp_memory,
    align_database_with_dp_memory, align_est2genome_database_with_dp_memory,
    align_est2genome_low_memory, align_est2genome_with_dp_memory, align_low_memory,
    align_protein_database_with_dp_memory, align_protein_low_memory,
    align_protein_to_genome_database_with_dp_memory, align_protein_to_genome_low_memory,
    align_protein_to_genome_with_dp_memory, align_protein_with_dp_memory, align_with_dp_memory,
};

/// Generic model graph, independent from a particular DP memory layout.
pub mod model {
    use super::Score;
    use std::collections::VecDeque;
    pub type StateId = u16;
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub enum Scope {
        Corner,
        Query,
        Anywhere,
        Edge,
    }
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub enum TiePolicy {
        Earliest,
        Latest,
    }
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub enum Label {
        None,
        Match,
        Gap,
        Ner,
        Splice5,
        Splice3,
        Intron,
        SplitCodon,
        Frameshift,
    }
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub enum ScoreKernel {
        Constant(Score),
        GapOpen,
        GapExtend,
        DnaSubstitution,
        ProteinSubstitution,
        CodonSubstitution,
        SpliceSite,
        IntronOpen,
        IntronClose {
            min_len: u32,
            max_len: u32,
        },
        Frameshift,
        Phase,
        NerSpan {
            min_len: u32,
            max_len: u32,
            open: Score,
        },
    }
    #[derive(Clone, Debug)]
    pub struct Transition {
        pub from: StateId,
        pub to: StateId,
        pub query_advance: u32,
        pub target_advance: u32,
        pub kernel: ScoreKernel,
        pub label: Label,
    }
    #[derive(Clone, Debug)]
    pub struct ModelIr {
        pub scope: Scope,
        pub tie_policy: TiePolicy,
        pub state_count: StateId,
        pub start: StateId,
        pub end: StateId,
        pub transitions: Box<[Transition]>,
    }
    #[derive(Clone, Debug, Eq, PartialEq)]
    pub enum ModelError {
        Empty,
        UnknownState,
        UnreachableEnd,
        EpsilonCycle,
        InvalidIntronBounds,
    }
    impl ModelIr {
        pub fn validate(&self) -> Result<(), ModelError> {
            if self.state_count == 0 {
                return Err(ModelError::Empty);
            }
            if self.start >= self.state_count || self.end >= self.state_count {
                return Err(ModelError::UnknownState);
            }
            for e in self.transitions.iter() {
                if e.from >= self.state_count || e.to >= self.state_count {
                    return Err(ModelError::UnknownState);
                }
                match e.kernel {
                    ScoreKernel::IntronClose { min_len, max_len }
                    | ScoreKernel::NerSpan {
                        min_len, max_len, ..
                    } if min_len > max_len => {
                        return Err(ModelError::InvalidIntronBounds);
                    }
                    _ => {}
                }
            }
            let mut seen = vec![false; self.state_count as usize];
            let mut queue = VecDeque::from([self.start]);
            while let Some(s) = queue.pop_front() {
                if seen[s as usize] {
                    continue;
                }
                seen[s as usize] = true;
                for e in self.transitions.iter().filter(|e| e.from == s) {
                    queue.push_back(e.to);
                }
            }
            if !seen[self.end as usize] {
                return Err(ModelError::UnreachableEnd);
            }
            let mut color = vec![0_u8; self.state_count as usize];
            fn visit(s: StateId, model: &ModelIr, color: &mut [u8]) -> bool {
                color[s as usize] = 1;
                for e in model
                    .transitions
                    .iter()
                    .filter(|e| e.from == s && e.query_advance == 0 && e.target_advance == 0)
                {
                    if color[e.to as usize] == 1
                        || (color[e.to as usize] == 0 && visit(e.to, model, color))
                    {
                        return true;
                    }
                }
                color[s as usize] = 2;
                false
            }
            for s in 0..self.state_count {
                if color[s as usize] == 0 && visit(s, self, &mut color) {
                    return Err(ModelError::EpsilonCycle);
                }
            }
            Ok(())
        }
    }
    pub fn affine(scope: Scope) -> ModelIr {
        let x = |from, to, query_advance, target_advance, kernel, label| Transition {
            from,
            to,
            query_advance,
            target_advance,
            kernel,
            label,
        };
        ModelIr {
            scope,
            tie_policy: TiePolicy::Latest,
            state_count: 5,
            start: 0,
            end: 4,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 1, 1, ScoreKernel::DnaSubstitution, Label::Match),
                x(1, 2, 1, 0, ScoreKernel::GapOpen, Label::Gap),
                x(1, 3, 0, 1, ScoreKernel::GapOpen, Label::Gap),
                x(2, 2, 1, 0, ScoreKernel::GapExtend, Label::Gap),
                x(3, 3, 0, 1, ScoreKernel::GapExtend, Label::Gap),
                x(2, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(3, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 4, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }
    /// Local affine model with one bounded non-equivalenced-region span.
    pub fn ner(min_len: u32, max_len: u32, open: Score) -> ModelIr {
        let mut model = affine(Scope::Anywhere);
        let mut transitions = model.transitions.into_vec();
        transitions.push(Transition {
            from: 1,
            to: 1,
            query_advance: 1,
            target_advance: 1,
            kernel: ScoreKernel::NerSpan {
                min_len,
                max_len,
                open,
            },
            label: Label::Ner,
        });
        model.transitions = transitions.into_boxed_slice();
        model
    }

    /// Structural IR for the two-strand local `est2genome` model. Each match
    /// orientation has its own target-intron path (5'SS/3'SS or reverse).
    pub fn est2genome(min_intron: u32, max_intron: u32) -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: 7,
            start: 0,
            end: 6,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(0, 2, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 1, 1, ScoreKernel::DnaSubstitution, Label::Match),
                x(2, 2, 1, 1, ScoreKernel::DnaSubstitution, Label::Match),
                x(1, 3, 0, 2, ScoreKernel::IntronOpen, Label::Splice5),
                x(3, 3, 0, 1, ScoreKernel::Constant(0), Label::Intron),
                x(
                    3,
                    1,
                    0,
                    2,
                    ScoreKernel::IntronClose {
                        min_len: min_intron,
                        max_len: max_intron,
                    },
                    Label::Splice3,
                ),
                x(2, 4, 0, 2, ScoreKernel::IntronOpen, Label::Splice3),
                x(4, 4, 0, 1, ScoreKernel::Constant(0), Label::Intron),
                x(
                    4,
                    2,
                    0,
                    2,
                    ScoreKernel::IntronClose {
                        min_len: min_intron,
                        max_len: max_intron,
                    },
                    Label::Splice5,
                ),
                x(1, 6, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(2, 6, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }

    /// Minimal local query-intron graph used by generic C4 composition and tests.
    pub fn query_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: 4,
            start: 0,
            end: 3,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 1, 1, ScoreKernel::DnaSubstitution, Label::Match),
                x(1, 2, 2, 0, ScoreKernel::IntronOpen, Label::Splice5),
                x(2, 2, 1, 0, ScoreKernel::Constant(0), Label::Intron),
                x(
                    2,
                    1,
                    2,
                    0,
                    ScoreKernel::IntronClose {
                        min_len: min_intron,
                        max_len: max_intron,
                    },
                    Label::Splice3,
                ),
                x(1, 3, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }

    /// Minimal local joint-intron graph. A compiled joint long state allows
    /// query and target intron lengths to vary independently within bounds.
    pub fn joint_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: 4,
            start: 0,
            end: 3,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 1, 1, ScoreKernel::DnaSubstitution, Label::Match),
                x(1, 2, 2, 2, ScoreKernel::IntronOpen, Label::Splice5),
                x(2, 2, 1, 1, ScoreKernel::Constant(0), Label::Intron),
                x(
                    2,
                    1,
                    2,
                    2,
                    ScoreKernel::IntronClose {
                        min_len: min_intron,
                        max_len: max_intron,
                    },
                    Label::Splice3,
                ),
                x(1, 3, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }

    /// Codon-to-codon phase model with an intron on the query, or jointly on
    /// query and target. `Phase` post edges recover the pre-intron bases from
    /// the long-state parent shadow before scoring the split codon.
    fn codon_phase_intron(on_target: bool, min_intron: u32, max_intron: u32) -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        let mut transitions = vec![
            x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
            x(1, 1, 3, 3, ScoreKernel::CodonSubstitution, Label::Match),
            x(1, 2, 0, 0, ScoreKernel::Constant(0), Label::None),
        ];
        let mut next = 3;
        for phase in 0..3 {
            let pre = if phase == 0 {
                1
            } else {
                let state = next;
                next += 1;
                transitions.push(x(
                    1,
                    state,
                    phase,
                    phase,
                    ScoreKernel::Constant(0),
                    Label::SplitCodon,
                ));
                state
            };
            let donor = next;
            let loop_state = next + 1;
            let acceptor = next + 2;
            next += 3;
            let target_advance = u32::from(on_target) * 2;
            transitions.push(x(
                pre,
                donor,
                2,
                target_advance,
                ScoreKernel::IntronOpen,
                Label::Splice5,
            ));
            transitions.push(x(
                donor,
                loop_state,
                0,
                0,
                ScoreKernel::Constant(0),
                Label::None,
            ));
            transitions.push(x(
                loop_state,
                loop_state,
                1,
                u32::from(on_target),
                ScoreKernel::Constant(0),
                Label::Intron,
            ));
            transitions.push(x(
                loop_state,
                acceptor,
                2,
                target_advance,
                ScoreKernel::IntronClose {
                    min_len: min_intron,
                    max_len: max_intron,
                },
                Label::Splice3,
            ));
            if phase == 0 {
                transitions.push(x(acceptor, 1, 0, 0, ScoreKernel::Constant(0), Label::None));
            } else {
                transitions.push(x(
                    acceptor,
                    1,
                    3 - phase,
                    3 - phase,
                    ScoreKernel::Phase,
                    Label::SplitCodon,
                ));
            }
        }
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: next,
            start: 0,
            end: 2,
            transitions: transitions.into_boxed_slice(),
        }
    }

    pub fn query_codon_phase_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        codon_phase_intron(false, min_intron, max_intron)
    }

    pub fn joint_codon_phase_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        codon_phase_intron(true, min_intron, max_intron)
    }

    /// Target-phased intron subgraph used by protein2genome. It contains the
    /// upstream 0:0, 1:2, and 2:1 phase paths.
    pub fn target_phase_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        let intron = |pre, donor, loop_state, acceptor, post, target_advance| {
            [
                x(pre, donor, 0, 2, ScoreKernel::IntronOpen, Label::Splice5),
                x(
                    donor,
                    loop_state,
                    0,
                    0,
                    ScoreKernel::Constant(0),
                    Label::None,
                ),
                x(
                    loop_state,
                    loop_state,
                    0,
                    1,
                    ScoreKernel::Constant(0),
                    Label::Intron,
                ),
                x(
                    loop_state,
                    acceptor,
                    0,
                    2,
                    ScoreKernel::IntronClose {
                        min_len: min_intron,
                        max_len: max_intron,
                    },
                    Label::Splice3,
                ),
                x(
                    acceptor,
                    post,
                    1,
                    target_advance,
                    ScoreKernel::CodonSubstitution,
                    Label::SplitCodon,
                ),
            ]
        };
        let mut transitions = Vec::new();
        // phase 0: no pre-intron bases, then a full codon after the intron
        transitions.push(x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None));
        transitions.extend(intron(1, 2, 3, 4, 10, 3));
        transitions.last_mut().expect("phase-0 post edge").label = Label::Match;
        // phase 1: one base before the intron, two after it
        transitions.push(x(0, 5, 0, 1, ScoreKernel::Constant(0), Label::SplitCodon));
        transitions.extend(intron(5, 6, 7, 8, 10, 2));
        // phase 2: two bases before the intron, one after it
        transitions.push(x(0, 9, 0, 2, ScoreKernel::Constant(0), Label::SplitCodon));
        transitions.extend(intron(9, 11, 12, 13, 10, 1));
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: 14,
            start: 0,
            end: 10,
            transitions: transitions.into_boxed_slice(),
        }
    }

    /// Protein-to-target-phase graph with a positive coding prefix before an intron.
    pub fn protein_phase_intron(min_intron: u32, max_intron: u32) -> ModelIr {
        let mut model = target_phase_intron(min_intron, max_intron);
        let mut transitions = model.transitions.into_vec();
        transitions.push(Transition {
            from: 1,
            to: 1,
            query_advance: 1,
            target_advance: 3,
            kernel: ScoreKernel::CodonSubstitution,
            label: Label::Match,
        });
        transitions.push(Transition {
            from: 1,
            to: 5,
            query_advance: 0,
            target_advance: 1,
            kernel: ScoreKernel::Constant(0),
            label: Label::SplitCodon,
        });
        transitions.push(Transition {
            from: 1,
            to: 9,
            query_advance: 0,
            target_advance: 2,
            kernel: ScoreKernel::Constant(0),
            label: Label::SplitCodon,
        });
        model.transitions = transitions.into_boxed_slice();
        model
    }

    /// Local codon-to-codon affine graph with frameshifts on both DNA chains.
    pub fn coding_to_coding() -> ModelIr {
        let x = |from, to, q, t, kernel, label| Transition {
            from,
            to,
            query_advance: q,
            target_advance: t,
            kernel,
            label,
        };
        ModelIr {
            scope: Scope::Anywhere,
            tie_policy: TiePolicy::Latest,
            state_count: 8,
            start: 0,
            end: 7,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 3, 3, ScoreKernel::CodonSubstitution, Label::Match),
                x(1, 2, 3, 0, ScoreKernel::GapOpen, Label::Gap),
                x(2, 2, 3, 0, ScoreKernel::GapExtend, Label::Gap),
                x(2, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 3, 0, 3, ScoreKernel::GapOpen, Label::Gap),
                x(3, 3, 0, 3, ScoreKernel::GapExtend, Label::Gap),
                x(3, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 4, 1, 0, ScoreKernel::Frameshift, Label::Frameshift),
                x(1, 4, 2, 0, ScoreKernel::Frameshift, Label::Frameshift),
                x(4, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(4, 1, 3, 0, ScoreKernel::Constant(0), Label::Frameshift),
                x(1, 5, 0, 1, ScoreKernel::Frameshift, Label::Frameshift),
                x(1, 5, 0, 2, ScoreKernel::Frameshift, Label::Frameshift),
                x(5, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(5, 1, 0, 3, ScoreKernel::Constant(0), Label::Frameshift),
                x(1, 7, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }

    /// Upstream `protein2dna` transition graph. The scalar backend does not
    /// execute this graph yet; it is the authoritative IR for the forthcoming
    /// codon/frameshift Viterbi kernel.
    pub fn protein_to_dna(scope: Scope) -> ModelIr {
        let x = |from, to, query_advance, target_advance, kernel, label| Transition {
            from,
            to,
            query_advance,
            target_advance,
            kernel,
            label,
        };
        ModelIr {
            scope,
            tie_policy: TiePolicy::Latest,
            state_count: 6,
            start: 0,
            end: 5,
            transitions: [
                x(0, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 1, 1, 3, ScoreKernel::CodonSubstitution, Label::Match),
                x(1, 2, 1, 0, ScoreKernel::GapOpen, Label::Gap),
                x(1, 3, 0, 3, ScoreKernel::GapOpen, Label::Gap),
                x(2, 2, 1, 0, ScoreKernel::GapExtend, Label::Gap),
                x(3, 3, 0, 3, ScoreKernel::GapExtend, Label::Gap),
                x(2, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(3, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(1, 4, 0, 1, ScoreKernel::Frameshift, Label::Frameshift),
                x(1, 4, 0, 2, ScoreKernel::Frameshift, Label::Frameshift),
                x(4, 1, 0, 0, ScoreKernel::Constant(0), Label::None),
                x(4, 1, 0, 3, ScoreKernel::Constant(0), Label::Frameshift),
                x(1, 5, 0, 0, ScoreKernel::Constant(0), Label::None),
            ]
            .into(),
        }
    }

    pub fn ungapped_translated(scope: Scope) -> ModelIr {
        ModelIr {
            scope,
            tie_policy: TiePolicy::Earliest,
            state_count: 3,
            start: 0,
            end: 2,
            transitions: [
                Transition {
                    from: 0,
                    to: 1,
                    query_advance: 0,
                    target_advance: 0,
                    kernel: ScoreKernel::Constant(0),
                    label: Label::None,
                },
                Transition {
                    from: 1,
                    to: 1,
                    query_advance: 3,
                    target_advance: 3,
                    kernel: ScoreKernel::CodonSubstitution,
                    label: Label::Match,
                },
                Transition {
                    from: 1,
                    to: 2,
                    query_advance: 0,
                    target_advance: 0,
                    kernel: ScoreKernel::Constant(0),
                    label: Label::None,
                },
            ]
            .into(),
        }
    }

    pub fn ungapped(scope: Scope) -> ModelIr {
        ModelIr {
            scope,
            tie_policy: TiePolicy::Earliest,
            state_count: 3,
            start: 0,
            end: 2,
            transitions: [
                Transition {
                    from: 0,
                    to: 1,
                    query_advance: 0,
                    target_advance: 0,
                    kernel: ScoreKernel::Constant(0),
                    label: Label::None,
                },
                Transition {
                    from: 1,
                    to: 1,
                    query_advance: 1,
                    target_advance: 1,
                    kernel: ScoreKernel::DnaSubstitution,
                    label: Label::Match,
                },
                Transition {
                    from: 1,
                    to: 2,
                    query_advance: 0,
                    target_advance: 0,
                    kernel: ScoreKernel::Constant(0),
                    label: Label::None,
                },
            ]
            .into(),
        }
    }
}
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

mod nucleic_matrix {
    include!("nucleic_matrix.rs");
}
mod blosum62_matrix {
    include!("blosum62_matrix.rs");
}

pub type Score = i32;
const NEG_INF: Score = i32::MIN / 4;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sequence {
    pub id: String,
    pub bases: Vec<u8>,
}

#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    Fasta(String),
    InvalidModel(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Fasta(e) | Self::InvalidModel(e) => f.write_str(e),
        }
    }
}
impl std::error::Error for Error {}
impl From<io::Error> for Error {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

/// Parse FASTA records without treating sequence bytes as UTF-8.
pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<Sequence>, Error> {
    let reader = BufReader::new(File::open(path)?);
    let mut records = Vec::new();
    let mut id: Option<String> = None;
    let mut bases = Vec::new();
    for line in reader.split(b'\n') {
        let mut line = line?;
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            if let Some(id) = id.take() {
                if bases.is_empty() {
                    return Err(Error::Fasta(format!("empty FASTA record {id}")));
                }
                records.push(Sequence {
                    id,
                    bases: std::mem::take(&mut bases),
                });
            }
            let header = &line[1..];
            let token = header
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if token.is_empty() {
                return Err(Error::Fasta("FASTA header has no identifier".into()));
            }
            id = Some(String::from_utf8_lossy(token).into_owned());
        } else if !line.is_empty() {
            if id.is_none() {
                return Err(Error::Fasta("sequence before first FASTA header".into()));
            }
            bases.extend(line.into_iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    if let Some(id) = id {
        if bases.is_empty() {
            return Err(Error::Fasta(format!("empty FASTA record {id}")));
        }
        records.push(Sequence { id, bases });
    }
    if records.is_empty() {
        return Err(Error::Fasta("no FASTA records found".into()));
    }
    Ok(records)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Model {
    Ungapped,
    Global,
    BestFit,
    Local,
    Overlap,
}
impl std::str::FromStr for Model {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ungapped" | "u" => Ok(Self::Ungapped),
            "affine:global" | "a:g" => Ok(Self::Global),
            "affine:bestfit" | "a:b" => Ok(Self::BestFit),
            "affine:local" | "a:l" => Ok(Self::Local),
            "affine:overlap" | "a:o" => Ok(Self::Overlap),
            _ => Err(Error::InvalidModel(format!(
                "unsupported model {s:?}; expected ungapped or affine:{{global,bestfit,local,overlap}}"
            ))),
        }
    }
}

impl Model {
    pub fn ir(self) -> model::ModelIr {
        match self {
            Self::Ungapped => model::ungapped(model::Scope::Anywhere),
            Self::Global => model::affine(model::Scope::Corner),
            Self::BestFit => model::affine(model::Scope::Query),
            Self::Local => model::affine(model::Scope::Anywhere),
            Self::Overlap => model::affine(model::Scope::Edge),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Scoring {
    pub gap_open: Score,
    pub gap_extend: Score,
    pub codon_gap_open: Score,
    pub codon_gap_extend: Score,
    pub frameshift: Score,
}
impl Default for Scoring {
    fn default() -> Self {
        Self {
            gap_open: -12,
            gap_extend: -4,
            codon_gap_open: -18,
            codon_gap_extend: -8,
            frameshift: -28,
        }
    }
}

/// Monotonic maximum queue for one target-intron DP row. Callers insert a
/// candidate when it becomes eligible; `expire_before` then enforces the
/// maximum intron length without scanning every possible start.
#[derive(Debug, Default)]
pub struct IntronWindow {
    entries: VecDeque<(usize, Score)>,
}
impl IntronWindow {
    pub fn insert(&mut self, start: usize, score: Score) {
        while self.entries.back().is_some_and(|entry| entry.1 <= score) {
            self.entries.pop_back();
        }
        self.entries.push_back((start, score));
    }
    pub fn expire_before(&mut self, start: usize) {
        while self.entries.front().is_some_and(|entry| entry.0 < start) {
            self.entries.pop_front();
        }
    }
    pub fn best(&self) -> Option<(usize, Score)> {
        self.entries.front().copied()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IntronScoring {
    pub min_len: u32,
    pub max_len: u32,
    pub open_penalty: Score,
    pub force_gtag: bool,
}
impl Default for IntronScoring {
    fn default() -> Self {
        Self {
            min_len: 30,
            max_len: 200_000,
            open_penalty: -30,
            force_gtag: false,
        }
    }
}

/// A spliced intron contains two-base donor and acceptor signals, so its
/// canonical traceback cannot represent a span shorter than four bases.
/// Upstream accepts smaller `--minintron` values but simply has no such
/// traceback candidate; clamp the internal eligibility span to avoid unsigned
/// underflow while preserving the caller's maximum length.
fn canonical_intron_scoring(mut intron: IntronScoring) -> IntronScoring {
    intron.min_len = intron.min_len.max(4);
    intron
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    Unknown,
    Forward,
    Reverse,
}
impl Strand {
    pub fn symbol(self) -> char {
        match self {
            Self::Unknown => '.',
            Self::Forward => '+',
            Self::Reverse => '-',
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Op {
    Match,
    Insert,
    Delete,
    Frameshift,
    Splice5,
    Splice3,
    Intron,
    SplitCodon,
    Ner,
}
impl Op {
    fn advances(self) -> (u32, u32) {
        match self {
            Self::Match => (1, 1),
            Self::Insert => (1, 0),
            Self::Delete => (0, 1),
            // Frameshift runs carry their explicit nucleotide advance in TraceRun.
            Self::Frameshift
            | Self::Splice5
            | Self::Splice3
            | Self::Intron
            | Self::SplitCodon
            | Self::Ner => (0, 0),
        }
    }
    fn cigar(self, query_advance: u32, target_advance: u32) -> char {
        if query_advance == 0 {
            'D'
        } else if target_advance == 0 {
            'I'
        } else {
            'M'
        }
    }
    fn vulgar(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insert | Self::Delete => 'G',
            Self::Frameshift => 'F',
            Self::Splice5 => '5',
            Self::Splice3 => '3',
            Self::Intron => 'I',
            Self::SplitCodon => 'S',
            Self::Ner => 'N',
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TraceRun {
    pub transition_id: u16,
    pub op: Op,
    pub query_advance: u32,
    pub target_advance: u32,
    pub repeats: u64,
}

/// One uncompressed model transition in traceback order. Unlike `TraceRun`,
/// this preserves zero-length epsilon transitions and the score contributed by
/// each individual transition, as required by Exonerate's `%P` RYO fields.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RawStep {
    pub transition_id: u16,
    pub query_advance: u32,
    pub target_advance: u32,
    pub score: Score,
}

/// A fixed-size group of atomic transitions stored by a DP parent. Intron
/// edges use all three slots (5' splice, loop, 3' splice); ordinary edges
/// use one slot. This keeps traceback storage compact and allocation-free.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TraceFragment {
    atoms: [Option<TraceRun>; 5],
}
impl TraceFragment {
    pub const fn empty() -> Self {
        Self {
            atoms: [None, None, None, None, None],
        }
    }
    pub const fn one(atom: TraceRun) -> Self {
        Self {
            atoms: [Some(atom), None, None, None, None],
        }
    }
    pub const fn intron(splice5: TraceRun, loop_run: TraceRun, splice3: TraceRun) -> Self {
        Self {
            atoms: [Some(splice5), Some(loop_run), Some(splice3), None, None],
        }
    }
    pub const fn phase_intron(
        pre: TraceRun,
        splice5: TraceRun,
        loop_run: TraceRun,
        splice3: TraceRun,
        post: TraceRun,
    ) -> Self {
        Self {
            atoms: [
                Some(pre),
                Some(splice5),
                Some(loop_run),
                Some(splice3),
                Some(post),
            ],
        }
    }
    pub const fn intron_match(
        splice5: TraceRun,
        loop_run: TraceRun,
        splice3: TraceRun,
        matched: TraceRun,
    ) -> Self {
        Self {
            atoms: [
                Some(splice5),
                Some(loop_run),
                Some(splice3),
                Some(matched),
                None,
            ],
        }
    }
    pub fn advances(self) -> (usize, usize) {
        self.atoms
            .into_iter()
            .flatten()
            .fold((0usize, 0usize), |(query, target), run| {
                (
                    query + run.query_advance as usize * run.repeats as usize,
                    target + run.target_advance as usize * run.repeats as usize,
                )
            })
    }
    pub fn append_to(self, trace: &mut Vec<TraceRun>) {
        for atom in self.atoms.into_iter().flatten() {
            if let Some(last) = trace.last_mut()
                && last.op == atom.op
                && last.transition_id == atom.transition_id
                && last.query_advance == atom.query_advance
                && last.target_advance == atom.target_advance
            {
                last.repeats += atom.repeats;
            } else {
                trace.push(atom);
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct Alignment {
    pub query_id: String,
    pub target_id: String,
    pub query_start: u64,
    pub query_end: u64,
    pub query_strand: Strand,
    pub target_start: u64,
    pub target_end: u64,
    pub target_len: u64,
    pub target_strand: Strand,
    pub score: Score,
    pub raw_trace: Vec<RawStep>,
    pub trace: Vec<TraceRun>,
}

impl Alignment {
    pub fn sugar(&self) -> String {
        format!(
            "sugar: {} {} {} {} {} {} {} {} {}",
            self.query_id,
            self.query_start,
            self.query_end,
            self.query_strand.symbol(),
            self.target_id,
            self.target_start,
            self.target_end,
            self.target_strand.symbol(),
            self.score
        )
    }
    pub fn cigar(&self) -> String {
        let mut ops: Vec<(char, u64)> = Vec::new();
        for run in &self.trace {
            let op = run.op.cigar(run.query_advance, run.target_advance);
            let length = u64::from(run.query_advance.max(run.target_advance)) * run.repeats;
            if let Some((last, count)) = ops.last_mut()
                && *last == op
            {
                *count += length;
            } else {
                ops.push((op, length));
            }
        }
        let body = ops
            .into_iter()
            .map(|(op, count)| format!("{} {}", op, count))
            .collect::<Vec<_>>()
            .join(" ");
        if body.is_empty() {
            self.sugar().replacen("sugar:", "cigar:", 1)
        } else {
            // Upstream's implicit zero-length start operation leaves two
            // spaces before a leading match, but only one before a leading
            // gap in a global traceback.
            let separator = if matches!(self.trace.first().map(|run| run.op), Some(Op::Match)) {
                "  "
            } else {
                " "
            };
            format!(
                "{}{}{}",
                self.sugar().replacen("sugar:", "cigar:", 1),
                separator,
                body
            )
        }
    }

    /// CIGAR rendering for the composed cDNA-to-genome model.
    ///
    /// Upstream keeps a separate CIGAR `M` segment when a nucleotide match
    /// meets a 3:3 coding match, even though both render as `M`.  This leaves
    /// the generic CIGAR contract unchanged for the other models.
    pub fn cdna_cigar(&self) -> String {
        let mut ops: Vec<(char, u64, u32, u32)> = Vec::new();
        for run in &self.trace {
            let op = run.op.cigar(run.query_advance, run.target_advance);
            let length = u64::from(run.query_advance.max(run.target_advance)) * run.repeats;
            if let Some((last, count, query_advance, target_advance)) = ops.last_mut()
                && *last == op
                && (op != 'M'
                    || (*query_advance == run.query_advance
                        && *target_advance == run.target_advance))
            {
                *count += length;
            } else {
                ops.push((op, length, run.query_advance, run.target_advance));
            }
        }
        let body = ops
            .into_iter()
            .map(|(op, count, _, _)| format!("{} {}", op, count))
            .collect::<Vec<_>>()
            .join(" ");
        if body.is_empty() {
            self.sugar().replacen("sugar:", "cigar:", 1)
        } else {
            let separator = if matches!(self.trace.first().map(|run| run.op), Some(Op::Match)) {
                "  "
            } else {
                " "
            };
            format!(
                "{}{}{}",
                self.sugar().replacen("sugar:", "cigar:", 1),
                separator,
                body
            )
        }
    }
    /// Render a compact GFF3 match record plus contiguous aligned parts.
    pub fn gff3(&self) -> String {
        let strand = self.target_strand.symbol();
        let (start, end) = if self.target_start <= self.target_end {
            (self.target_start + 1, self.target_end)
        } else {
            (self.target_end + 1, self.target_start)
        };
        let id = format!("match_{}_{}", self.query_id, self.target_id);
        let mut lines = vec![format!(
            "{}\texonerate-rs\tmatch\t{}\t{}\t{}\t{}\t.\tID={};Target={} {} {}",
            self.target_id,
            start,
            end,
            self.score,
            strand,
            id,
            self.query_id,
            self.query_start + 1,
            self.query_end
        )];
        let mut q = self.query_start;
        let mut t = if self.target_strand == Strand::Reverse {
            self.target_start.max(self.target_end)
        } else {
            self.target_start.min(self.target_end)
        };
        let mut part = 0_u32;
        for run in &self.trace {
            let qa = u64::from(run.query_advance) * run.repeats;
            let ta = u64::from(run.target_advance) * run.repeats;
            if qa > 0 && ta > 0 {
                part += 1;
                let (a, b) = if self.target_strand == Strand::Reverse {
                    (t.saturating_sub(ta) + 1, t)
                } else {
                    (t + 1, t + ta)
                };
                lines.push(format!("{}\texonerate-rs\tmatch_part\t{}\t{}\t.\t{}\t.\tID={}.part{};Parent={};Target={} {} {}", self.target_id, a.min(b), a.max(b), strand, id, part, id, self.query_id, q + 1, q + qa));
            }
            q += qa;
            if self.target_strand == Strand::Reverse {
                t = t.saturating_sub(ta);
            } else {
                t += ta;
            }
        }
        lines.join("\n")
    }
    /** Render the same alignment as GFF3 features projected onto the query. */
    pub fn query_gff3(&self) -> String {
        let strand = self.query_strand.symbol();
        let (start, end) = if self.query_start <= self.query_end {
            (self.query_start + 1, self.query_end)
        } else {
            (self.query_end + 1, self.query_start)
        };
        let id = format!("match_{}_{}", self.target_id, self.query_id);
        let mut lines = vec![format!(
            "{}\texonerate-rs\tmatch\t{}\t{}\t{}\t{}\t.\tID={};Target={} {} {} {}",
            self.query_id,
            start,
            end,
            self.score,
            strand,
            id,
            self.target_id,
            self.target_start.min(self.target_end) + 1,
            self.target_start.max(self.target_end),
            self.target_strand.symbol()
        )];
        let mut q = if self.query_strand == Strand::Reverse {
            self.query_start.max(self.query_end)
        } else {
            self.query_start.min(self.query_end)
        };
        let mut t = if self.target_strand == Strand::Reverse {
            self.target_start.max(self.target_end)
        } else {
            self.target_start.min(self.target_end)
        };
        let mut part = 0_u32;
        for run in &self.trace {
            let qa = u64::from(run.query_advance) * run.repeats;
            let ta = u64::from(run.target_advance) * run.repeats;
            if qa > 0 && ta > 0 {
                part += 1;
                let q_next = if self.query_strand == Strand::Reverse {
                    q.saturating_sub(qa)
                } else {
                    q + qa
                };
                let t_next = if self.target_strand == Strand::Reverse {
                    t.saturating_sub(ta)
                } else {
                    t + ta
                };
                lines.push(format!(
                    "{}\texonerate-rs\tmatch_part\t{}\t{}\t.\t{}\t.\tID={}.part{};Parent={};Target={} {} {} {}",
                    self.query_id,
                    q.min(q_next) + 1,
                    q.max(q_next),
                    strand,
                    id,
                    part,
                    id,
                    self.target_id,
                    t.min(t_next) + 1,
                    t.max(t_next),
                    self.target_strand.symbol()
                ));
            }
            q = if self.query_strand == Strand::Reverse {
                q.saturating_sub(qa)
            } else {
                q + qa
            };
            t = if self.target_strand == Strand::Reverse {
                t.saturating_sub(ta)
            } else {
                t + ta
            };
        }
        lines.join("\n")
    }

    /// Render a deterministic human-readable alignment block. Complex biological
    /// transitions are expanded according to their coordinate advances.
    pub fn pretty(&self, query: &Sequence, target: &Sequence) -> String {
        let query_bases = if self.query_strand == Strand::Reverse {
            reverse_complement(&query.bases)
        } else {
            query.bases.clone()
        };
        let target_bases = if self.target_strand == Strand::Reverse {
            reverse_complement(&target.bases)
        } else {
            target.bases.clone()
        };
        let mut query_position = if self.query_strand == Strand::Reverse {
            query
                .bases
                .len()
                .saturating_sub(self.query_start.max(self.query_end) as usize)
        } else {
            self.query_start.min(self.query_end) as usize
        };
        let mut target_position = if self.target_strand == Strand::Reverse {
            target
                .bases
                .len()
                .saturating_sub(self.target_start.max(self.target_end) as usize)
        } else {
            self.target_start.min(self.target_end) as usize
        };
        let mut query_row = Vec::new();
        let mut target_row = Vec::new();
        for run in &self.trace {
            for _ in 0..run.repeats {
                let query_advance = run.query_advance as usize;
                let target_advance = run.target_advance as usize;
                let width = query_advance.max(target_advance);
                for offset in 0..width {
                    query_row.push(if offset < query_advance {
                        query_bases
                            .get(query_position + offset)
                            .copied()
                            .unwrap_or(b'?')
                    } else {
                        b'-'
                    });
                    target_row.push(if offset < target_advance {
                        target_bases
                            .get(target_position + offset)
                            .copied()
                            .unwrap_or(b'?')
                    } else {
                        b'-'
                    });
                }
                query_position += query_advance;
                target_position += target_advance;
            }
        }
        let middle = query_row
            .iter()
            .zip(&target_row)
            .map(|(query_base, target_base)| {
                if query_base.eq_ignore_ascii_case(target_base) && *query_base != b'-' {
                    b'|'
                } else {
                    b' '
                }
            })
            .collect::<Vec<_>>();
        let mut output = format!(
            "C4 Alignment:\n------------\n         Query: {}\n        Target: {}\n     Raw score: {}\n   Query range: {} -> {}\n  Target range: {} -> {}\n",
            self.query_id,
            self.target_id,
            self.score,
            self.query_start,
            self.query_end,
            self.target_start,
            self.target_end
        );
        for offset in (0..query_row.len()).step_by(60) {
            let end = (offset + 60).min(query_row.len());
            output.push_str(&format!(
                "\n{}\n{}\n{}\n",
                String::from_utf8_lossy(&query_row[offset..end]),
                String::from_utf8_lossy(&middle[offset..end]),
                String::from_utf8_lossy(&target_row[offset..end])
            ));
        }
        output
    }

    pub fn vulgar(&self) -> String {
        let mut s = self.sugar().replacen("sugar:", "vulgar:", 1);
        let mut current: Option<(char, u64, u64)> = None;
        for run in &self.trace {
            let (query, target) = (run.query_advance as u64, run.target_advance as u64);
            let label = if run.op == Op::Match && run.query_advance == 3 && run.target_advance == 3
            {
                'C'
            } else {
                run.op.vulgar()
            };
            let item = (label, query * run.repeats, target * run.repeats);
            match current {
                Some((last, q, t))
                    if last == item.0 && (q > 0 || item.1 == 0) && (t > 0 || item.2 == 0) =>
                {
                    current = Some((last, q + item.1, t + item.2))
                }
                Some((last, q, t)) => {
                    s.push_str(&format!(" {} {} {}", last, q, t));
                    current = Some(item);
                }
                None => current = Some(item),
            }
        }
        if let Some((label, query, target)) = current {
            s.push_str(&format!(" {} {} {}", label, query, target));
        }
        s
    }
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        x => x,
    }
}
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&b| complement(b)).collect()
}

/// Direction-specific splice-site predictor matching Exonerate's bundled
/// primate 5'/3' frequency matrices. Scores use the same log transform and
/// C-style rounding as `SplicePredictor_predict_array_int`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SpliceType {
    DonorForward,
    AcceptorForward,
    DonorReverse,
    AcceptorReverse,
}
const SPLICE_5: [[u8; 4]; 9] = [
    [28, 40, 17, 14],
    [59, 14, 13, 14],
    [8, 5, 81, 6],
    [0, 0, 100, 0],
    [0, 0, 0, 100],
    [54, 2, 42, 2],
    [74, 8, 11, 8],
    [5, 6, 85, 4],
    [16, 18, 21, 45],
];
const SPLICE_3: [[u8; 4]; 15] = [
    [10, 31, 14, 44],
    [8, 36, 14, 43],
    [6, 34, 12, 48],
    [6, 34, 8, 52],
    [9, 37, 9, 45],
    [9, 38, 10, 44],
    [8, 44, 9, 40],
    [9, 41, 8, 41],
    [6, 44, 6, 45],
    [6, 40, 6, 48],
    [23, 28, 26, 23],
    [2, 79, 1, 18],
    [100, 0, 0, 0],
    [0, 0, 100, 0],
    [28, 14, 47, 11],
];
fn splice_base_index(base: u8, reverse: bool) -> Option<usize> {
    match (base.to_ascii_uppercase(), reverse) {
        (b'A', false) | (b'T', true) => Some(0),
        (b'C', false) | (b'G', true) => Some(1),
        (b'G', false) | (b'C', true) => Some(2),
        (b'T' | b'U', false) | (b'A', true) => Some(3),
        _ => None,
    }
}
/// Score a splice boundary at `position`; `None` denotes an invalid coordinate
/// or a forced non-canonical dinucleotide.
pub fn splice_score(
    sequence: &[u8],
    position: usize,
    kind: SpliceType,
    force_gtag: bool,
) -> Option<Score> {
    let (matrix_5, reverse, splice_after, expected) = match kind {
        SpliceType::DonorForward => (true, false, 3usize, *b"GT"),
        SpliceType::AcceptorForward => (false, false, 12usize, *b"AG"),
        SpliceType::DonorReverse => (true, true, 4usize, *b"AC"),
        SpliceType::AcceptorReverse => (false, true, 1usize, *b"CT"),
    };
    if position + 2 > sequence.len() {
        return None;
    }
    if force_gtag
        && [
            sequence[position].to_ascii_uppercase(),
            sequence[position + 1].to_ascii_uppercase(),
        ] != expected
    {
        return None;
    }
    let matrix: &[[u8; 4]] = if matrix_5 { &SPLICE_5 } else { &SPLICE_3 };
    let start = position.saturating_sub(splice_after);
    let offset = splice_after.saturating_sub(position);
    let mut value = 0.0_f64;
    for row in offset..matrix.len() {
        let sequence_pos = start + row - offset;
        if sequence_pos >= sequence.len() {
            break;
        }
        let source_row = if reverse { matrix.len() - 1 - row } else { row };
        if let Some(base) = splice_base_index(sequence[sequence_pos], reverse) {
            value += (((f64::from(matrix[source_row][base]) + 1.0) / 26.0).ln()) * 1.5;
        }
    }
    Some(if value < 0.0 {
        (value - 0.5) as Score
    } else {
        (value + 0.5) as Score
    })
}

pub fn translate_dna(bases: &[u8], frame: u8) -> Vec<u8> {
    const CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    fn base_index(base: u8) -> Option<usize> {
        match base.to_ascii_uppercase() {
            b'T' | b'U' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    let mut protein = Vec::new();
    let mut pos = frame as usize;
    while pos + 3 <= bases.len() {
        let amino = match (
            base_index(bases[pos]),
            base_index(bases[pos + 1]),
            base_index(bases[pos + 2]),
        ) {
            (Some(a), Some(b), Some(c)) => CODE[a * 16 + b * 4 + c],
            _ => b'X',
        };
        protein.push(amino);
        pos += 3;
    }
    protein
}

/// Score a codon split by an intron. `left` and `right` must jointly contain
/// three nucleotides; the returned score uses the same translated amino acid
/// and BLOSUM62 lookup as protein2dna match transitions.
pub fn split_codon_score(query_amino: u8, left: &[u8], right: &[u8]) -> Option<Score> {
    if left.len() + right.len() != 3 {
        return None;
    }
    let mut codon = [0_u8; 3];
    codon[..left.len()].copy_from_slice(left);
    codon[left.len()..].copy_from_slice(right);
    Some(protein_score(
        query_amino,
        translate_dna(&codon, 0)[0],
        Scoring::default(),
    ))
}

fn nucleic_index(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'R' => 1,
        b'N' => 2,
        b'D' => 3,
        b'C' | b'U' => 4,
        b'Q' => 5,
        b'E' => 6,
        b'G' => 7,
        b'H' => 8,
        b'I' => 9,
        b'L' => 10,
        b'K' => 11,
        b'M' => 12,
        b'F' => 13,
        b'P' => 14,
        b'S' => 15,
        b'T' => 16,
        b'W' => 17,
        b'Y' => 18,
        b'V' => 19,
        b'B' => 20,
        b'Z' => 21,
        b'*' => 23,
        _ => 22,
    }
}
fn dna_score(a: u8, b: u8, _scoring: Scoring) -> Score {
    nucleic_matrix::NUCLEIC[nucleic_index(a)][nucleic_index(b)]
}
/// Return the upstream nucleic substitution score used by RYO similarity
/// reporting. Gap penalties are deliberately not part of this lookup.
pub fn dna_substitution_score(a: u8, b: u8) -> Score {
    dna_score(a, b, Scoring::default())
}
fn protein_index(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'R' => 1,
        b'N' => 2,
        b'D' => 3,
        b'C' => 4,
        b'Q' => 5,
        b'E' => 6,
        b'G' => 7,
        b'H' => 8,
        b'I' => 9,
        b'L' => 10,
        b'K' => 11,
        b'M' => 12,
        b'F' => 13,
        b'P' => 14,
        b'S' => 15,
        b'T' => 16,
        b'W' => 17,
        b'Y' => 18,
        b'V' => 19,
        b'B' => 20,
        b'Z' => 21,
        b'*' => 23,
        _ => 22,
    }
}
fn protein_score(a: u8, b: u8, _scoring: Scoring) -> Score {
    blosum62_matrix::BLOSUM62[protein_index(a)][protein_index(b)]
}
/// Return the upstream BLOSUM62 substitution score used by RYO similarity
/// reporting. Gap penalties are deliberately not part of this lookup.
pub fn protein_substitution_score(a: u8, b: u8) -> Score {
    protein_score(a, b, Scoring::default())
}
pub fn dna_self_score(sequence: &Sequence) -> Score {
    sequence
        .bases
        .iter()
        .map(|&base| dna_score(base, base, Scoring::default()))
        .sum()
}
pub fn protein_self_score(sequence: &Sequence) -> Score {
    sequence
        .bases
        .iter()
        .map(|&base| protein_score(base, base, Scoring::default()))
        .sum()
}

/** Self-comparison score used by upstream advance-3 translated matches. */
pub fn translated_self_score(sequence: &Sequence) -> Score {
    sequence
        .bases
        .chunks_exact(3)
        .map(|codon| {
            let amino_acid = translate_dna(codon, 0)[0];
            protein_score(amino_acid, amino_acid, Scoring::default())
        })
        .sum()
}

fn kernel_score(kernel: model::ScoreKernel, a: u8, b: u8, scoring: Scoring) -> Score {
    match kernel {
        model::ScoreKernel::Constant(value) => value,
        model::ScoreKernel::GapOpen => scoring.gap_open,
        model::ScoreKernel::GapExtend => scoring.gap_extend,
        model::ScoreKernel::DnaSubstitution => dna_score(a, b, scoring),
        _ => NEG_INF,
    }
}
fn idx(i: usize, j: usize, cols: usize) -> usize {
    i * cols + j
}
fn add(a: Score, b: Score) -> Score {
    if a <= NEG_INF / 2 {
        NEG_INF
    } else {
        a.saturating_add(b)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum State {
    M,
    I,
    D,
    Stop,
}
impl State {
    fn rank(self) -> u8 {
        match self {
            Self::M => 3,
            Self::I => 2,
            Self::D => 1,
            Self::Stop => 0,
        }
    }
}
fn best(candidates: &[(Score, State)]) -> (Score, State) {
    *candidates
        .iter()
        .max_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.rank().cmp(&b.1.rank())))
        .expect("nonempty")
}

#[derive(Clone, Copy)]
struct P2Parent {
    state: State,
    op: Op,
    transition_id: u16,
    query_advance: u32,
    target_advance: u32,
}

#[derive(Clone, Copy)]
struct P2GenomeParent {
    state: State,
    fragment: TraceFragment,
}

fn p2_update(
    current: &mut (Score, State, Option<P2Parent>),
    score: Score,
    state: State,
    parent: P2Parent,
) {
    if score > current.0 || (score == current.0 && state.rank() > current.1.rank()) {
        *current = (score, state, Some(parent));
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IntronCandidate {
    pub start: usize,
    pub score: Score,
    pub state_rank: u8,
}

/// Bounded monotonic queue retaining the full predecessor identity needed by
/// an intron traceback, rather than only its score.
#[derive(Clone, Default)]
pub struct IntronCandidateWindow {
    entries: VecDeque<IntronCandidate>,
}
impl IntronCandidateWindow {
    pub fn insert(&mut self, candidate: IntronCandidate) {
        while self.entries.back().is_some_and(|entry| {
            entry.score < candidate.score
                || (entry.score == candidate.score && entry.state_rank <= candidate.state_rank)
        }) {
            self.entries.pop_back();
        }
        self.entries.push_back(candidate);
    }
    pub fn expire_before(&mut self, start: usize) {
        while self
            .entries
            .front()
            .is_some_and(|entry| entry.start < start)
        {
            self.entries.pop_front();
        }
    }
    pub fn best(&self) -> Option<IntronCandidate> {
        self.entries.front().copied()
    }

    /// Bytes needed to clone the live candidates into a checkpoint.
    ///
    /// `VecDeque::clone` only needs storage for live elements; spare ring
    /// capacity is not semantic state and is deliberately excluded.
    #[allow(dead_code)] // consumed when genome2genome checkpoint snapshots are wired in
    fn checkpoint_bytes(&self) -> usize {
        size_of::<Self>() + self.entries.len() * size_of::<IntronCandidate>()
    }
}

#[derive(Clone, Copy)]
struct EstParent {
    state: State,
    fragment: TraceFragment,
    raw_fragment: RawFragment,
}

fn state_from_rank(rank: u8) -> State {
    match rank {
        3 => State::M,
        2 => State::I,
        1 => State::D,
        _ => State::Stop,
    }
}

/// Failure while compiling or executing a generic model graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ExecutionError {
    InvalidModel(model::ModelError),
    UnsupportedKernel(model::ScoreKernel),
    UnsupportedCodonAdvance { query: u32, target: u32 },
    TooManyTransitions(usize),
}
impl fmt::Display for ExecutionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModel(error) => write!(f, "invalid model: {error:?}"),
            Self::UnsupportedKernel(kernel) => write!(f, "unsupported generic kernel: {kernel:?}"),
            Self::UnsupportedCodonAdvance { query, target } => {
                write!(f, "unsupported codon advance {query}:{target}")
            }
            Self::TooManyTransitions(count) => {
                write!(
                    f,
                    "model has {count} transitions; traceback IDs are limited to 65535"
                )
            }
        }
    }
}
impl std::error::Error for ExecutionError {}

#[derive(Clone, Copy)]
struct RawFragment {
    atoms: [Option<RawStep>; 5],
}
impl RawFragment {
    const fn one(step: RawStep) -> Self {
        Self {
            atoms: [Some(step), None, None, None, None],
        }
    }
    const fn three(first: RawStep, middle: RawStep, last: RawStep) -> Self {
        Self {
            atoms: [Some(first), Some(middle), Some(last), None, None],
        }
    }
    fn with_prefix(prefix: Option<u16>, steps: &[RawStep]) -> Self {
        let mut atoms = [None; 5];
        let mut index = 0;
        if let Some(transition_id) = prefix {
            atoms[0] = Some(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
            index = 1;
        }
        for step in steps.iter().copied() {
            atoms[index] = Some(step);
            index += 1;
        }
        Self { atoms }
    }
    fn append_to(self, trace: &mut Vec<RawStep>) {
        trace.extend(self.atoms.into_iter().flatten());
    }
}

#[derive(Clone, Copy)]
struct GenericParent {
    state: model::StateId,
    query: usize,
    target: usize,
    fragment: TraceFragment,
    raw_fragment: RawFragment,
}

#[derive(Clone, Copy)]
struct TargetLongState {
    source: model::StateId,
    destination: model::StateId,
    open_transition: u16,
    loop_transition: u16,
    close_transition: u16,
    min_len: usize,
    max_len: usize,
    donor: SpliceType,
    acceptor: SpliceType,
}

#[derive(Clone, Copy)]
struct NerLongState {
    source: model::StateId,
    destination: model::StateId,
    transition: u16,
    min_len: usize,
    max_len: usize,
    open: Score,
}

fn ner_long_states(ir: &model::ModelIr) -> Vec<NerLongState> {
    ir.transitions
        .iter()
        .enumerate()
        .filter_map(|(transition, edge)| {
            let model::ScoreKernel::NerSpan {
                min_len,
                max_len,
                open,
            } = edge.kernel
            else {
                return None;
            };
            Some(NerLongState {
                source: edge.from,
                destination: edge.to,
                transition: transition as u16,
                min_len: (min_len as usize).max(4),
                max_len: max_len as usize,
                open,
            })
        })
        .collect()
}

#[derive(Clone, Copy)]
struct JointLongState {
    source: model::StateId,
    destination: model::StateId,
    open_transition: u16,
    loop_transition: u16,
    close_transition: u16,
    min_len: usize,
    max_len: usize,
    donor: SpliceType,
    acceptor: SpliceType,
}

#[derive(Clone, Copy)]
struct QueryLongState {
    source: model::StateId,
    destination: model::StateId,
    open_transition: u16,
    loop_transition: u16,
    close_transition: u16,
    min_len: usize,
    max_len: usize,
    donor: SpliceType,
    acceptor: SpliceType,
}

fn epsilon_reaches(
    ir: &model::ModelIr,
    start: model::StateId,
    destination: model::StateId,
) -> bool {
    let mut seen = vec![false; ir.state_count as usize];
    let mut queue = VecDeque::from([start]);
    while let Some(state) = queue.pop_front() {
        if state == destination {
            return true;
        }
        if std::mem::replace(&mut seen[state as usize], true) {
            continue;
        }
        for edge in ir.transitions.iter().filter(|edge| {
            edge.from == state && edge.query_advance == 0 && edge.target_advance == 0
        }) {
            queue.push_back(edge.to);
        }
    }
    false
}

fn target_long_states(ir: &model::ModelIr) -> Vec<TargetLongState> {
    let mut result = Vec::new();
    for (close_id, close) in ir.transitions.iter().enumerate() {
        let model::ScoreKernel::IntronClose { min_len, max_len } = close.kernel else {
            continue;
        };
        if (close.query_advance, close.target_advance) != (0, 2) {
            continue;
        }
        let Some((loop_id, _)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            edge.from == close.from
                && edge.to == close.from
                && (edge.query_advance, edge.target_advance) == (0, 1)
                && edge.label == model::Label::Intron
        }) else {
            continue;
        };
        let Some((open_id, open)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            matches!(edge.kernel, model::ScoreKernel::IntronOpen)
                && (edge.query_advance, edge.target_advance) == (0, 2)
                && epsilon_reaches(ir, edge.to, close.from)
        }) else {
            continue;
        };
        let kinds = match (open.label, close.label) {
            (model::Label::Splice5, model::Label::Splice3) => {
                Some((SpliceType::DonorForward, SpliceType::AcceptorForward))
            }
            (model::Label::Splice3, model::Label::Splice5) => {
                Some((SpliceType::AcceptorReverse, SpliceType::DonorReverse))
            }
            _ => None,
        };
        if let Some((donor, acceptor)) = kinds {
            result.push(TargetLongState {
                source: open.from,
                destination: close.to,
                open_transition: open_id as u16,
                loop_transition: loop_id as u16,
                close_transition: close_id as u16,
                min_len: (min_len as usize).max(4),
                max_len: max_len as usize,
                donor,
                acceptor,
            });
        }
    }
    result
}

fn joint_long_states(ir: &model::ModelIr) -> Vec<JointLongState> {
    let mut result = Vec::new();
    for (close_id, close) in ir.transitions.iter().enumerate() {
        let model::ScoreKernel::IntronClose { min_len, max_len } = close.kernel else {
            continue;
        };
        if (close.query_advance, close.target_advance) != (2, 2) {
            continue;
        }
        let Some((loop_id, _)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            edge.from == close.from
                && edge.to == close.from
                && (edge.query_advance, edge.target_advance) == (1, 1)
                && edge.label == model::Label::Intron
        }) else {
            continue;
        };
        let Some((open_id, open)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            matches!(edge.kernel, model::ScoreKernel::IntronOpen)
                && (edge.query_advance, edge.target_advance) == (2, 2)
                && epsilon_reaches(ir, edge.to, close.from)
        }) else {
            continue;
        };
        let kinds = match (open.label, close.label) {
            (model::Label::Splice5, model::Label::Splice3) => {
                Some((SpliceType::DonorForward, SpliceType::AcceptorForward))
            }
            (model::Label::Splice3, model::Label::Splice5) => {
                Some((SpliceType::AcceptorReverse, SpliceType::DonorReverse))
            }
            _ => None,
        };
        if let Some((donor, acceptor)) = kinds {
            result.push(JointLongState {
                source: open.from,
                destination: close.to,
                open_transition: open_id as u16,
                loop_transition: loop_id as u16,
                close_transition: close_id as u16,
                min_len: (min_len as usize).max(4),
                max_len: max_len as usize,
                donor,
                acceptor,
            });
        }
    }
    result
}

fn query_long_states(ir: &model::ModelIr) -> Vec<QueryLongState> {
    let mut result = Vec::new();
    for (close_id, close) in ir.transitions.iter().enumerate() {
        let model::ScoreKernel::IntronClose { min_len, max_len } = close.kernel else {
            continue;
        };
        if (close.query_advance, close.target_advance) != (2, 0) {
            continue;
        }
        let Some((loop_id, _)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            edge.from == close.from
                && edge.to == close.from
                && (edge.query_advance, edge.target_advance) == (1, 0)
                && edge.label == model::Label::Intron
        }) else {
            continue;
        };
        let Some((open_id, open)) = ir.transitions.iter().enumerate().find(|(_, edge)| {
            matches!(edge.kernel, model::ScoreKernel::IntronOpen)
                && (edge.query_advance, edge.target_advance) == (2, 0)
                && epsilon_reaches(ir, edge.to, close.from)
        }) else {
            continue;
        };
        let kinds = match (open.label, close.label) {
            (model::Label::Splice5, model::Label::Splice3) => {
                Some((SpliceType::DonorForward, SpliceType::AcceptorForward))
            }
            (model::Label::Splice3, model::Label::Splice5) => {
                Some((SpliceType::AcceptorReverse, SpliceType::DonorReverse))
            }
            _ => None,
        };
        if let Some((donor, acceptor)) = kinds {
            result.push(QueryLongState {
                source: open.from,
                destination: close.to,
                open_transition: open_id as u16,
                loop_transition: loop_id as u16,
                close_transition: close_id as u16,
                min_len: min_len as usize,
                max_len: max_len as usize,
                donor,
                acceptor,
            });
        }
    }
    result
}

fn epsilon_order(ir: &model::ModelIr) -> Vec<model::StateId> {
    let mut indegree = vec![0usize; ir.state_count as usize];
    for edge in ir
        .transitions
        .iter()
        .filter(|edge| edge.query_advance == 0 && edge.target_advance == 0)
    {
        indegree[edge.to as usize] += 1;
    }
    let mut queue: VecDeque<_> = (0..ir.state_count)
        .filter(|state| indegree[*state as usize] == 0)
        .collect();
    let mut order = Vec::with_capacity(ir.state_count as usize);
    while let Some(state) = queue.pop_front() {
        order.push(state);
        for edge in ir.transitions.iter().filter(|edge| {
            edge.from == state && edge.query_advance == 0 && edge.target_advance == 0
        }) {
            indegree[edge.to as usize] -= 1;
            if indegree[edge.to as usize] == 0 {
                queue.push_back(edge.to);
            }
        }
    }
    order
}

fn generic_kernel_score(
    edge: &model::Transition,
    query: &[u8],
    target: &[u8],
    i: usize,
    j: usize,
    scoring: Scoring,
) -> Result<Score, ExecutionError> {
    let qa = edge.query_advance as usize;
    let ta = edge.target_advance as usize;
    Ok(match edge.kernel {
        model::ScoreKernel::Constant(value) => value,
        model::ScoreKernel::GapOpen => {
            if qa >= 3 || ta >= 3 {
                scoring.codon_gap_open
            } else {
                scoring.gap_open
            }
        }
        model::ScoreKernel::GapExtend => {
            if qa >= 3 || ta >= 3 {
                scoring.codon_gap_extend
            } else {
                scoring.gap_extend
            }
        }
        model::ScoreKernel::DnaSubstitution => dna_score(query[i - 1], target[j - 1], scoring),
        model::ScoreKernel::ProteinSubstitution => {
            protein_score(query[i - 1], target[j - 1], scoring)
        }
        model::ScoreKernel::CodonSubstitution => match (qa, ta) {
            (1, 3) => protein_score(
                query[i - 1],
                translate_dna(&target[j - 3..j], 0)[0],
                scoring,
            ),
            (3, 1) => protein_score(
                translate_dna(&query[i - 3..i], 0)[0],
                target[j - 1],
                scoring,
            ),
            (3, 3) => protein_score(
                translate_dna(&query[i - 3..i], 0)[0],
                translate_dna(&target[j - 3..j], 0)[0],
                scoring,
            ),
            _ => {
                return Err(ExecutionError::UnsupportedCodonAdvance {
                    query: edge.query_advance,
                    target: edge.target_advance,
                });
            }
        },
        model::ScoreKernel::Frameshift => scoring.frameshift,
        model::ScoreKernel::NerSpan { .. }
        | model::ScoreKernel::SpliceSite
        | model::ScoreKernel::IntronOpen
        | model::ScoreKernel::IntronClose { .. }
        | model::ScoreKernel::Phase => return Err(ExecutionError::UnsupportedKernel(edge.kernel)),
    })
}

fn generic_state_priority(ir: &model::ModelIr, state: model::StateId) -> u8 {
    let mut priority = 0;
    for edge in ir
        .transitions
        .iter()
        .filter(|edge| edge.from == state && edge.to == state)
    {
        priority = priority.max(
            match (edge.label, edge.query_advance, edge.target_advance) {
                (model::Label::Match, q, t) if q > 0 && t > 0 => 3,
                (model::Label::Gap, q, 0) if q > 0 => 2,
                (model::Label::Gap, 0, t) if t > 0 => 1,
                _ => 0,
            },
        );
    }
    priority
}

fn generic_op(edge: &model::Transition) -> Option<Op> {
    match edge.label {
        model::Label::None => None,
        model::Label::Match => Some(Op::Match),
        model::Label::Gap => {
            if edge.target_advance == 0 {
                Some(Op::Insert)
            } else {
                Some(Op::Delete)
            }
        }
        model::Label::Ner => Some(Op::Ner),
        model::Label::Splice5 => Some(Op::Splice5),
        model::Label::Splice3 => Some(Op::Splice3),
        model::Label::Intron => Some(Op::Intron),
        model::Label::SplitCodon => Some(Op::SplitCodon),
        model::Label::Frameshift => Some(Op::Frameshift),
    }
}

/// Execute a validated finite C4-style graph with deterministic Viterbi
/// traceback. Epsilon edges may form a DAG. Bounded introns are compiled into monotonic long states; phase post edges use
/// parent-coordinate shadows to score non-contiguous split codons.
pub fn align_model_ir(
    query: &Sequence,
    target: &Sequence,
    ir: &model::ModelIr,
    scoring: Scoring,
    strand: Strand,
) -> Result<Alignment, ExecutionError> {
    align_model_ir_with_intron(query, target, ir, scoring, IntronScoring::default(), strand)
}

pub fn align_model_ir_with_intron(
    query: &Sequence,
    target: &Sequence,
    ir: &model::ModelIr,
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
) -> Result<Alignment, ExecutionError> {
    align_model_ir_with_intron_forbidden(query, target, ir, scoring, intron, strand, None)
}

fn align_model_ir_with_intron_forbidden(
    query: &Sequence,
    target: &Sequence,
    ir: &model::ModelIr,
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Result<Alignment, ExecutionError> {
    ir.validate().map_err(ExecutionError::InvalidModel)?;
    if ir.transitions.len() > u16::MAX as usize {
        return Err(ExecutionError::TooManyTransitions(ir.transitions.len()));
    }
    let target_bases = if strand == Strand::Reverse {
        reverse_complement(&target.bases)
    } else {
        target.bases.clone()
    };
    let protein_pair = ir
        .transitions
        .iter()
        .any(|edge| matches!(edge.kernel, model::ScoreKernel::ProteinSubstitution));
    let protein_query = protein_pair
        || ir.transitions.iter().any(|edge| {
            matches!(edge.kernel, model::ScoreKernel::CodonSubstitution)
                && edge.query_advance == 1
                && edge.target_advance == 3
        });
    let protein_target = protein_pair
        || ir.transitions.iter().any(|edge| {
            matches!(edge.kernel, model::ScoreKernel::CodonSubstitution)
                && edge.query_advance == 3
                && edge.target_advance == 1
        });
    let (n, m, states) = (
        query.bases.len(),
        target_bases.len(),
        ir.state_count as usize,
    );
    let cells = (n + 1) * (m + 1);
    let mut scores = vec![NEG_INF; cells * states];
    let mut parents: Vec<Option<GenericParent>> = vec![None; cells * states];
    let slot =
        |i: usize, j: usize, state: model::StateId| idx(i, j, m + 1) * states + state as usize;
    let legal_start = |i: usize, j: usize| match ir.scope {
        model::Scope::Corner => i == 0 && j == 0,
        model::Scope::Query => i == 0,
        model::Scope::Anywhere => true,
        model::Scope::Edge => i == 0 || j == 0,
    };
    for i in 0..=n {
        for j in 0..=m {
            if legal_start(i, j) {
                scores[slot(i, j, ir.start)] = 0;
            }
        }
    }
    let order = epsilon_order(ir);
    let long_states = target_long_states(ir);
    let query_long_states = query_long_states(ir);
    let joint_long_states = joint_long_states(ir);
    let ner_long_states = ner_long_states(ir);
    let mut compiled = vec![false; ir.transitions.len()];
    for long in &long_states {
        compiled[long.open_transition as usize] = true;
        compiled[long.loop_transition as usize] = true;
        compiled[long.close_transition as usize] = true;
    }
    for long in &query_long_states {
        compiled[long.open_transition as usize] = true;
        compiled[long.loop_transition as usize] = true;
        compiled[long.close_transition as usize] = true;
    }
    for long in &joint_long_states {
        compiled[long.open_transition as usize] = true;
        compiled[long.loop_transition as usize] = true;
        compiled[long.close_transition as usize] = true;
    }
    for long in &ner_long_states {
        compiled[long.transition as usize] = true;
    }
    let mut outgoing = vec![Vec::<(u16, &model::Transition)>::new(); states];
    for (transition, edge) in ir.transitions.iter().enumerate() {
        outgoing[edge.from as usize].push((transition as u16, edge));
    }
    let mut ner_long_columns: Vec<Vec<JointIntronWindow>> = ner_long_states
        .iter()
        .map(|_| (0..=m).map(|_| JointIntronWindow::default()).collect())
        .collect();
    let mut joint_long_columns: Vec<Vec<JointIntronWindow>> = joint_long_states
        .iter()
        .map(|_| (0..=m).map(|_| JointIntronWindow::default()).collect())
        .collect();
    let mut query_long_windows: Vec<Vec<IntronCandidateWindow>> = query_long_states
        .iter()
        .map(|_| (0..=m).map(|_| IntronCandidateWindow::default()).collect())
        .collect();
    for i in 0..=n {
        for (long_index, long) in ner_long_states.iter().enumerate() {
            if i >= long.min_len {
                let query_start = i - long.min_len;
                for target_start in 0..=m {
                    let source_score = scores[slot(query_start, target_start, long.source)];
                    if source_score > NEG_INF / 2 {
                        ner_long_columns[long_index][target_start].insert(JointIntronCandidate {
                            query_start,
                            target_start,
                            score: add(source_score, long.open),
                            state_rank: 0,
                        });
                    }
                }
            }
            if i > long.max_len {
                for column in &mut ner_long_columns[long_index] {
                    column.expire_before(i - long.max_len);
                }
            }
        }
        for (long_index, long) in joint_long_states.iter().enumerate() {
            if i >= long.min_len {
                let query_start = i - long.min_len;
                if let Some(query_donor) =
                    splice_score(&query.bases, query_start, long.donor, intron.force_gtag)
                {
                    for target_start in 0..=m {
                        let source_score = scores[slot(query_start, target_start, long.source)];
                        if source_score > 0 {
                            if let Some(target_donor) = splice_score(
                                &target_bases,
                                target_start,
                                long.donor,
                                intron.force_gtag,
                            ) {
                                joint_long_columns[long_index][target_start].insert(
                                    JointIntronCandidate {
                                        query_start,
                                        target_start,
                                        score: add(
                                            source_score,
                                            intron
                                                .open_penalty
                                                .saturating_add(query_donor)
                                                .saturating_add(target_donor),
                                        ),
                                        state_rank: 0,
                                    },
                                );
                            }
                        }
                    }
                }
            }
            if i > long.max_len {
                for column in &mut joint_long_columns[long_index] {
                    column.expire_before(i - long.max_len);
                }
            }
        }
        let mut ner_long_targets: Vec<JointIntronWindow> = ner_long_states
            .iter()
            .map(|_| JointIntronWindow::default())
            .collect();
        let mut joint_long_targets: Vec<JointIntronWindow> = joint_long_states
            .iter()
            .map(|_| JointIntronWindow::default())
            .collect();
        let mut long_windows: Vec<IntronCandidateWindow> = long_states
            .iter()
            .map(|_| IntronCandidateWindow::default())
            .collect();
        for j in 0..=m {
            for (long_index, long) in ner_long_states.iter().enumerate() {
                if j >= long.min_len {
                    let target_start = j - long.min_len;
                    if let Some(candidate) = ner_long_columns[long_index][target_start].best() {
                        ner_long_targets[long_index].insert(candidate);
                    }
                }
                if j > long.max_len {
                    while ner_long_targets[long_index]
                        .entries
                        .front()
                        .is_some_and(|entry| entry.target_start < j - long.max_len)
                    {
                        ner_long_targets[long_index].entries.pop_front();
                    }
                }
                if let Some(candidate) = ner_long_targets[long_index].best() {
                    let destination = slot(i, j, long.destination);
                    if candidate.score > scores[destination] {
                        scores[destination] = candidate.score;
                        parents[destination] = Some(GenericParent {
                            state: long.source,
                            query: candidate.query_start,
                            target: candidate.target_start,
                            fragment: TraceFragment::one(TraceRun {
                                transition_id: long.transition,
                                op: Op::Ner,
                                query_advance: (i - candidate.query_start) as u32,
                                target_advance: (j - candidate.target_start) as u32,
                                repeats: 1,
                            }),
                            raw_fragment: RawFragment::one(RawStep {
                                transition_id: long.transition,
                                query_advance: (i - candidate.query_start) as u32,
                                target_advance: (j - candidate.target_start) as u32,
                                score: long.open,
                            }),
                        });
                    }
                }
            }
            for (long_index, long) in joint_long_states.iter().enumerate() {
                if j >= long.min_len {
                    let target_start = j - long.min_len;
                    if let Some(candidate) = joint_long_columns[long_index][target_start].best() {
                        joint_long_targets[long_index].insert(candidate);
                    }
                }
                if j > long.max_len {
                    while joint_long_targets[long_index]
                        .entries
                        .front()
                        .is_some_and(|entry| entry.target_start < j - long.max_len)
                    {
                        joint_long_targets[long_index].entries.pop_front();
                    }
                }
                if i >= 2 && j >= 2 {
                    if let (Some(candidate), Some(query_acceptor), Some(target_acceptor)) = (
                        joint_long_targets[long_index].best(),
                        splice_score(&query.bases, i - 2, long.acceptor, intron.force_gtag),
                        splice_score(&target_bases, j - 2, long.acceptor, intron.force_gtag),
                    ) {
                        let destination = slot(i, j, long.destination);
                        let score = add(add(candidate.score, query_acceptor), target_acceptor);
                        if score > scores[destination] {
                            scores[destination] = score;
                            parents[destination] = Some(GenericParent {
                                state: long.source,
                                query: candidate.query_start,
                                target: candidate.target_start,
                                fragment: TraceFragment::intron(
                                    TraceRun {
                                        transition_id: long.open_transition,
                                        op: if long.donor == SpliceType::DonorForward {
                                            Op::Splice5
                                        } else {
                                            Op::Splice3
                                        },
                                        query_advance: 2,
                                        target_advance: 2,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.loop_transition,
                                        op: Op::Intron,
                                        query_advance: (i - candidate.query_start - 4) as u32,
                                        target_advance: (j - candidate.target_start - 4) as u32,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.close_transition,
                                        op: if long.acceptor == SpliceType::AcceptorForward {
                                            Op::Splice3
                                        } else {
                                            Op::Splice5
                                        },
                                        query_advance: 2,
                                        target_advance: 2,
                                        repeats: 1,
                                    },
                                ),
                                raw_fragment: RawFragment::three(
                                    RawStep {
                                        transition_id: long.open_transition,
                                        query_advance: 2,
                                        target_advance: 2,
                                        score: intron
                                            .open_penalty
                                            .saturating_add(
                                                splice_score(
                                                    &query.bases,
                                                    candidate.query_start,
                                                    long.donor,
                                                    intron.force_gtag,
                                                )
                                                .unwrap_or(0),
                                            )
                                            .saturating_add(
                                                splice_score(
                                                    &target_bases,
                                                    candidate.target_start,
                                                    long.donor,
                                                    intron.force_gtag,
                                                )
                                                .unwrap_or(0),
                                            ),
                                    },
                                    RawStep {
                                        transition_id: long.loop_transition,
                                        query_advance: (i - candidate.query_start - 4) as u32,
                                        target_advance: (j - candidate.target_start - 4) as u32,
                                        score: 0,
                                    },
                                    RawStep {
                                        transition_id: long.close_transition,
                                        query_advance: 2,
                                        target_advance: 2,
                                        score: query_acceptor.saturating_add(target_acceptor),
                                    },
                                ),
                            });
                        }
                    }
                }
            }
            for (long_index, long) in query_long_states.iter().enumerate() {
                if i >= long.min_len {
                    let source_i = i - long.min_len;
                    let source_score = scores[slot(source_i, j, long.source)];
                    if source_score > 0 {
                        if let Some(donor) =
                            splice_score(&query.bases, source_i, long.donor, intron.force_gtag)
                        {
                            query_long_windows[long_index][j].insert(IntronCandidate {
                                start: source_i,
                                score: add(source_score, intron.open_penalty.saturating_add(donor)),
                                state_rank: 0,
                            });
                        }
                    }
                }
                if i > long.max_len {
                    query_long_windows[long_index][j].expire_before(i - long.max_len);
                }
                if i >= 2 {
                    if let (Some(candidate), Some(acceptor)) = (
                        query_long_windows[long_index][j].best(),
                        splice_score(&query.bases, i - 2, long.acceptor, intron.force_gtag),
                    ) {
                        let destination = slot(i, j, long.destination);
                        let score = add(candidate.score, acceptor);
                        if score > scores[destination] {
                            scores[destination] = score;
                            parents[destination] = Some(GenericParent {
                                state: long.source,
                                query: candidate.start,
                                target: j,
                                fragment: TraceFragment::intron(
                                    TraceRun {
                                        transition_id: long.open_transition,
                                        op: if long.donor == SpliceType::DonorForward {
                                            Op::Splice5
                                        } else {
                                            Op::Splice3
                                        },
                                        query_advance: 2,
                                        target_advance: 0,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.loop_transition,
                                        op: Op::Intron,
                                        query_advance: (i - candidate.start - 4) as u32,
                                        target_advance: 0,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.close_transition,
                                        op: if long.acceptor == SpliceType::AcceptorForward {
                                            Op::Splice3
                                        } else {
                                            Op::Splice5
                                        },
                                        query_advance: 2,
                                        target_advance: 0,
                                        repeats: 1,
                                    },
                                ),
                                raw_fragment: RawFragment::three(
                                    RawStep {
                                        transition_id: long.open_transition,
                                        query_advance: 2,
                                        target_advance: 0,
                                        score: intron.open_penalty.saturating_add(
                                            splice_score(
                                                &query.bases,
                                                candidate.start,
                                                long.donor,
                                                intron.force_gtag,
                                            )
                                            .unwrap_or(0),
                                        ),
                                    },
                                    RawStep {
                                        transition_id: long.loop_transition,
                                        query_advance: (i - candidate.start - 4) as u32,
                                        target_advance: 0,
                                        score: 0,
                                    },
                                    RawStep {
                                        transition_id: long.close_transition,
                                        query_advance: 2,
                                        target_advance: 0,
                                        score: acceptor,
                                    },
                                ),
                            });
                        }
                    }
                }
            }
            for (long_index, long) in long_states.iter().enumerate() {
                if j >= long.min_len {
                    let source_j = j - long.min_len;
                    let source_score = scores[slot(i, source_j, long.source)];
                    if source_score > 0 {
                        if let Some(donor) =
                            splice_score(&target_bases, source_j, long.donor, intron.force_gtag)
                        {
                            long_windows[long_index].insert(IntronCandidate {
                                start: source_j,
                                score: add(source_score, intron.open_penalty.saturating_add(donor)),
                                state_rank: 0,
                            });
                        }
                    }
                }
                if j > long.max_len {
                    long_windows[long_index].expire_before(j - long.max_len);
                }
                if j >= 2 {
                    if let (Some(candidate), Some(acceptor)) = (
                        long_windows[long_index].best(),
                        splice_score(&target_bases, j - 2, long.acceptor, intron.force_gtag),
                    ) {
                        let destination = slot(i, j, long.destination);
                        let score = add(candidate.score, acceptor);
                        if score > scores[destination] {
                            scores[destination] = score;
                            parents[destination] = Some(GenericParent {
                                state: long.source,
                                query: i,
                                target: candidate.start,
                                fragment: TraceFragment::intron(
                                    TraceRun {
                                        transition_id: long.open_transition,
                                        op: if long.donor == SpliceType::DonorForward {
                                            Op::Splice5
                                        } else {
                                            Op::Splice3
                                        },
                                        query_advance: 0,
                                        target_advance: 2,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.loop_transition,
                                        op: Op::Intron,
                                        query_advance: 0,
                                        target_advance: (j - candidate.start - 4) as u32,
                                        repeats: 1,
                                    },
                                    TraceRun {
                                        transition_id: long.close_transition,
                                        op: if long.acceptor == SpliceType::AcceptorForward {
                                            Op::Splice3
                                        } else {
                                            Op::Splice5
                                        },
                                        query_advance: 0,
                                        target_advance: 2,
                                        repeats: 1,
                                    },
                                ),
                                raw_fragment: RawFragment::three(
                                    RawStep {
                                        transition_id: long.open_transition,
                                        query_advance: 0,
                                        target_advance: 2,
                                        score: intron.open_penalty.saturating_add(
                                            splice_score(
                                                &target_bases,
                                                candidate.start,
                                                long.donor,
                                                intron.force_gtag,
                                            )
                                            .unwrap_or(0),
                                        ),
                                    },
                                    RawStep {
                                        transition_id: long.loop_transition,
                                        query_advance: 0,
                                        target_advance: (j - candidate.start - 4) as u32,
                                        score: 0,
                                    },
                                    RawStep {
                                        transition_id: long.close_transition,
                                        query_advance: 0,
                                        target_advance: 2,
                                        score: acceptor,
                                    },
                                ),
                            });
                        }
                    }
                }
            }
            for &state in &order {
                let source_score = scores[slot(i, j, state)];
                if source_score <= NEG_INF / 2 {
                    continue;
                }
                'edges: for &(transition, edge) in &outgoing[state as usize] {
                    if compiled[transition as usize] {
                        continue;
                    }
                    let ni = i + edge.query_advance as usize;
                    let nj = j + edge.target_advance as usize;
                    if ni > n || nj > m {
                        continue;
                    }
                    let scores_pair = edge.label == model::Label::Match
                        || (edge.label == model::Label::SplitCodon
                            && matches!(
                                edge.kernel,
                                model::ScoreKernel::Phase | model::ScoreKernel::CodonSubstitution
                            ));
                    if scores_pair
                        && forbidden_match_transition(
                            forbidden,
                            i,
                            j,
                            edge.query_advance as usize,
                            edge.target_advance as usize,
                        )
                    {
                        continue;
                    }
                    let is_partial_codon = matches!(edge.kernel, model::ScoreKernel::Phase)
                        || (matches!(edge.kernel, model::ScoreKernel::CodonSubstitution)
                            && !matches!(
                                (edge.query_advance, edge.target_advance),
                                (1, 3) | (3, 1) | (3, 3)
                            ));
                    let edge_score = if is_partial_codon {
                        let (mut shadow_i, mut shadow_j, mut shadow_state) = (i, j, state);
                        let (donor_i, donor_j) = loop {
                            let Some(parent) = parents[slot(shadow_i, shadow_j, shadow_state)]
                            else {
                                continue 'edges;
                            };
                            if parent.query != shadow_i || parent.target != shadow_j {
                                break (parent.query, parent.target);
                            }
                            (shadow_i, shadow_j, shadow_state) =
                                (parent.query, parent.target, parent.state);
                        };
                        let query_is_protein = edge.query_advance == 1 && protein_query;
                        let target_is_protein = edge.target_advance == 1 && protein_target;
                        let query_amino = if query_is_protein {
                            query.bases[ni - 1]
                        } else {
                            let pre = 3 - edge.query_advance as usize;
                            if donor_i < pre {
                                continue 'edges;
                            }
                            let mut codon = Vec::with_capacity(3);
                            codon.extend_from_slice(&query.bases[donor_i - pre..donor_i]);
                            codon.extend_from_slice(&query.bases[i..ni]);
                            translate_dna(&codon, 0)[0]
                        };
                        let target_amino = if target_is_protein {
                            target_bases[nj - 1]
                        } else {
                            let pre = 3 - edge.target_advance as usize;
                            if donor_j < pre {
                                continue 'edges;
                            }
                            let mut codon = Vec::with_capacity(3);
                            codon.extend_from_slice(&target_bases[donor_j - pre..donor_j]);
                            codon.extend_from_slice(&target_bases[j..nj]);
                            translate_dna(&codon, 0)[0]
                        };
                        protein_score(query_amino, target_amino, scoring)
                    } else {
                        generic_kernel_score(edge, &query.bases, &target_bases, ni, nj, scoring)?
                    };
                    let candidate = add(source_score, edge_score);
                    let destination = slot(ni, nj, edge.to);
                    let current_priority = parents[destination]
                        .map_or(0, |parent| generic_state_priority(ir, parent.state));
                    let candidate_priority = generic_state_priority(ir, state);
                    if candidate > scores[destination]
                        || (candidate == scores[destination]
                            && candidate_priority > current_priority)
                    {
                        scores[destination] = candidate;
                        let fragment = generic_op(edge).map_or_else(TraceFragment::empty, |op| {
                            TraceFragment::one(TraceRun {
                                transition_id: transition,
                                op,
                                query_advance: edge.query_advance,
                                target_advance: edge.target_advance,
                                repeats: 1,
                            })
                        });
                        parents[destination] = Some(GenericParent {
                            state,
                            query: i,
                            target: j,
                            fragment,
                            raw_fragment: RawFragment::one(RawStep {
                                transition_id: transition,
                                query_advance: edge.query_advance,
                                target_advance: edge.target_advance,
                                score: edge_score,
                            }),
                        });
                    }
                }
            }
        }
    }
    let legal_end = |i: usize, j: usize| match ir.scope {
        model::Scope::Corner => i == n && j == m,
        model::Scope::Query => i == n,
        model::Scope::Anywhere => true,
        model::Scope::Edge => i == n || j == m,
    };
    let mut end = (0usize, 0usize, NEG_INF);
    for i in 0..=n {
        for j in 0..=m {
            if legal_end(i, j) {
                let score = scores[slot(i, j, ir.end)];
                let replace = score > end.2
                    || (score == end.2
                        && ir.tie_policy == model::TiePolicy::Latest
                        && (i, j) > (end.0, end.1));
                if replace {
                    end = (i, j, score);
                }
            }
        }
    }
    let (mut i, mut j, score) = end;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut state = ir.end;
    let mut reversed = Vec::new();
    while let Some(parent) = parents[slot(i, j, state)] {
        reversed.push((parent.fragment, parent.raw_fragment));
        (i, j, state) = (parent.query, parent.target, parent.state);
    }
    reversed.reverse();
    let mut trace = Vec::new();
    let mut raw_trace = Vec::new();
    for (fragment, raw_fragment) in reversed {
        fragment.append_to(&mut trace);
        raw_fragment.append_to(&mut raw_trace);
    }
    let (target_start, target_end) = if strand == Strand::Reverse {
        (
            target.bases.len() as u64 - j as u64,
            target.bases.len() as u64 - oriented_target_end,
        )
    } else {
        (j as u64, oriented_target_end)
    };
    Ok(Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end,
        query_strand: if protein_query {
            Strand::Unknown
        } else {
            Strand::Forward
        },
        target_start,
        target_end,
        target_len: target.bases.len() as u64,
        target_strand: if protein_target {
            Strand::Unknown
        } else {
            strand
        },
        score,
        raw_trace,
        trace,
    })
}

/// Exhaustive DNA alignment.  The returned trace is a canonical run-length
/// representation, suitable for every output formatter.
pub fn align(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
) -> Alignment {
    align_with_scorer(query, target, model, scoring, strand, dna_score)
}

pub fn align_protein(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
) -> Alignment {
    let mut alignment = align_with_scorer(
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

fn align_with_scorer(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
    scorer: fn(u8, u8, Scoring) -> Score,
) -> Alignment {
    align_with_scorer_forbidden(query, target, model, scoring, strand, scorer, None)
}

fn align_with_scorer_forbidden(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
    scorer: fn(u8, u8, Scoring) -> Score,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let ir = model.ir();
    debug_assert!(ir.validate().is_ok(), "built-in model must be valid");
    let t = if strand == Strand::Forward {
        target.bases.clone()
    } else {
        reverse_complement(&target.bases)
    };
    let (n, m, cols) = (query.bases.len(), t.len(), t.len() + 1);
    if model == Model::Ungapped {
        return ungapped(query, target, &t, scoring, strand, scorer);
    }
    let gap_open = kernel_score(ir.transitions[2].kernel, b'A', b'A', scoring);
    let gap_extend = kernel_score(ir.transitions[4].kernel, b'A', b'A', scoring);
    let size = (n + 1) * (m + 1);
    let mut mm = vec![NEG_INF; size];
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut pm = vec![State::Stop; size];
    let mut pi = vec![State::Stop; size];
    let mut pd = vec![State::Stop; size];
    let local = model == Model::Local;
    let bestfit = model == Model::BestFit;
    let overlap = model == Model::Overlap;
    mm[0] = 0;
    for i in 0..=n {
        for j in 0..=m {
            let k = idx(i, j, cols);
            if i == 0 && j == 0 {
                continue;
            }
            if (local && (i == 0 || j == 0))
                || (bestfit && i == 0)
                || (overlap && (i == 0 || j == 0))
            {
                mm[k] = 0;
                pm[k] = State::Stop;
            }
            if model == Model::Global {
                if j == 0 && i > 0 {
                    let p = idx(i - 1, j, cols);
                    let (v, st) = best(&[
                        (add(mm[p], gap_open), State::M),
                        (add(ii[p], gap_extend), State::I),
                    ]);
                    ii[k] = v;
                    pi[k] = st;
                }
                if i == 0 && j > 0 {
                    let p = idx(i, j - 1, cols);
                    let (v, st) = best(&[
                        (add(mm[p], gap_open), State::M),
                        (add(dd[p], gap_extend), State::D),
                    ]);
                    dd[k] = v;
                    pd[k] = st;
                }
            }
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let p = idx(i - 1, j - 1, cols);
                let sub = scorer(query.bases[i - 1], t[j - 1], scoring);
                let (mut v, mut st) = best(&[
                    (add(mm[p], sub), State::M),
                    (add(ii[p], sub), State::I),
                    (add(dd[p], sub), State::D),
                ]);
                if local && v < 0 {
                    v = 0;
                    st = State::Stop;
                }
                mm[k] = v;
                pm[k] = st;
            }
            if i > 0 && !(model == Model::Global && j == 0) {
                let p = idx(i - 1, j, cols);
                let (mut v, mut st) = best(&[
                    (add(mm[p], gap_open), State::M),
                    (add(ii[p], gap_extend), State::I),
                ]);
                if local && v < 0 {
                    v = NEG_INF;
                    st = State::Stop;
                }
                ii[k] = v;
                pi[k] = st;
            }
            if j > 0 && !(model == Model::Global && i == 0) {
                let p = idx(i, j - 1, cols);
                let (mut v, mut st) = best(&[
                    (add(mm[p], gap_open), State::M),
                    (add(dd[p], gap_extend), State::D),
                ]);
                if local && v < 0 {
                    v = NEG_INF;
                    st = State::Stop;
                }
                dd[k] = v;
                pd[k] = st;
            }
        }
    }
    let mut end = (0usize, 0usize, State::M, NEG_INF);
    for i in 0..=n {
        for j in 0..=m {
            let valid = match model {
                Model::Global => i == n && j == m,
                Model::BestFit => i == n,
                Model::Local => true,
                Model::Overlap => i == n || j == m,
                Model::Ungapped => false,
            };
            if !valid {
                continue;
            }
            let k = idx(i, j, cols);
            for (v, st) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if v > end.3 || (v == end.3 && (i, j, st.rank()) > (end.0, end.1, end.2.rank())) {
                    end = (i, j, st, v);
                }
            }
        }
    }
    traceback(
        query, target, &t, scoring, scorer, strand, end, cols, &pm, &pi, &pd,
    )
}

#[allow(clippy::too_many_arguments)]
fn traceback(
    query: &Sequence,
    target: &Sequence,
    oriented_target: &[u8],
    scoring: Scoring,
    scorer: fn(u8, u8, Scoring) -> Score,
    strand: Strand,
    end: (usize, usize, State, Score),
    cols: usize,
    pm: &[State],
    pi: &[State],
    pd: &[State],
) -> Alignment {
    let (mut i, mut j, mut st, score) = end;
    let (qe, te) = (i as u64, j as u64);
    let mut ops = Vec::new();
    while st != State::Stop {
        let k = idx(i, j, cols);
        let prev = match st {
            State::M => {
                if i == 0 || j == 0 || pm[k] == State::Stop {
                    break;
                }
                i -= 1;
                j -= 1;
                pm[k]
            }
            State::I => {
                if i == 0 {
                    break;
                }
                i -= 1;
                pi[k]
            }
            State::D => {
                if j == 0 {
                    break;
                }
                j -= 1;
                pd[k]
            }
            State::Stop => break,
        };
        let op = match st {
            State::M => Op::Match,
            State::I => Op::Insert,
            State::D => Op::Delete,
            State::Stop => unreachable!(),
        };
        let transition_id = match (st, prev) {
            (State::M, _) => 1,
            (State::I, State::M) => 2,
            (State::I, _) => 4,
            (State::D, State::M) => 3,
            (State::D, _) => 5,
            (State::Stop, _) => unreachable!(),
        };
        ops.push((op, transition_id));
        st = prev;
    }
    ops.reverse();
    let (query_start, oriented_target_start) = (i, j);
    let mut raw_trace = Vec::new();
    if !ops.is_empty() {
        raw_trace.push(RawStep {
            transition_id: 0,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
        let (mut query_position, mut target_position) = (query_start, oriented_target_start);
        let mut previous_op = None;
        for &(op, transition_id) in &ops {
            let epsilon_id = match previous_op {
                Some(Op::Insert) if op != Op::Insert => Some(6),
                Some(Op::Delete) if op != Op::Delete => Some(7),
                _ => None,
            };
            if let Some(transition_id) = epsilon_id {
                raw_trace.push(RawStep {
                    transition_id,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                });
            }
            let (query_advance, target_advance) = op.advances();
            let step_score = match op {
                Op::Match => scorer(
                    query.bases[query_position],
                    oriented_target[target_position],
                    scoring,
                ),
                Op::Insert | Op::Delete if matches!(transition_id, 2 | 3) => scoring.gap_open,
                Op::Insert | Op::Delete => scoring.gap_extend,
                _ => unreachable!("affine traceback contains only match and gap operations"),
            };
            raw_trace.push(RawStep {
                transition_id,
                query_advance,
                target_advance,
                score: step_score,
            });
            query_position += query_advance as usize;
            target_position += target_advance as usize;
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
    let (ts, te) = if strand == Strand::Forward {
        (j as u64, te)
    } else {
        (
            target.bases.len() as u64 - j as u64,
            target.bases.len() as u64 - te,
        )
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: i as u64,
        query_end: qe,
        query_strand: Strand::Forward,
        target_start: ts,
        target_end: te,
        target_len: target.bases.len() as u64,
        target_strand: strand,
        score,
        raw_trace,
        trace,
    }
}

fn ungapped(
    query: &Sequence,
    target: &Sequence,
    t: &[u8],
    scoring: Scoring,
    strand: Strand,
    scorer: fn(u8, u8, Scoring) -> Score,
) -> Alignment {
    let mut best = (0usize, 0usize, 0usize, 0usize, 0);
    for qi in 0..query.bases.len() {
        for tj in 0..t.len() {
            let (mut i, mut j, mut score, mut qstart, mut tstart) = (qi, tj, 0, qi, tj);
            while i < query.bases.len() && j < t.len() {
                score += scorer(query.bases[i], t[j], scoring);
                if score < 0 {
                    score = 0;
                    qstart = i + 1;
                    tstart = j + 1;
                }
                if score > best.4 {
                    best = (qstart, tstart, i + 1, j + 1, score);
                }
                i += 1;
                j += 1;
            }
        }
    }
    let (qs, oriented_ts, qe, oriented_te, score) = best;
    let mut raw_trace = Vec::new();
    if score > 0 {
        raw_trace.push(RawStep {
            transition_id: 0,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
        for offset in 0..(qe - qs) {
            raw_trace.push(RawStep {
                transition_id: 1,
                query_advance: 1,
                target_advance: 1,
                score: scorer(query.bases[qs + offset], t[oriented_ts + offset], scoring),
            });
        }
        raw_trace.push(RawStep {
            transition_id: 2,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    let (ts, te) = if strand == Strand::Forward {
        (oriented_ts as u64, oriented_te as u64)
    } else {
        (
            target.bases.len() as u64 - oriented_ts as u64,
            target.bases.len() as u64 - oriented_te as u64,
        )
    };
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: qs as u64,
        query_end: qe as u64,
        query_strand: Strand::Forward,
        target_start: ts,
        target_end: te,
        target_len: target.bases.len() as u64,
        target_strand: strand,
        score,
        raw_trace,
        trace: if score > 0 {
            vec![TraceRun {
                transition_id: 1,
                op: Op::Match,
                query_advance: 1,
                target_advance: 1,
                repeats: (qe - qs) as u64,
            }]
        } else {
            vec![]
        },
    }
}

/// Local EST-to-genome alignment with affine gaps, bounded target introns,
/// splice-site scoring, and a full canonical traceback.
pub fn align_est2genome(query: &Sequence, target: &Sequence, intron: IntronScoring) -> Alignment {
    align_est2genome_stranded(query, target, intron, true)
}

/// EST-to-genome alignment, optionally restricted to the target forward strand.
pub fn align_est2genome_stranded(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    both_strands: bool,
) -> Alignment {
    let forward = align_est2genome_one(query, target, &target.bases, intron, Strand::Forward, None);
    if !both_strands {
        return forward;
    }
    let reverse_bases = reverse_complement(&target.bases);
    let reverse =
        align_est2genome_one(query, target, &reverse_bases, intron, Strand::Reverse, None);
    if reverse.score > forward.score {
        reverse
    } else {
        forward
    }
}

fn align_est2genome_one(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    intron: IntronScoring,
    strand: Strand,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let (n, m, cols) = (query.bases.len(), bases.len(), bases.len() + 1);
    let size = (n + 1) * (m + 1);
    let mut mm = vec![0; size];
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut pm: Vec<Option<EstParent>> = vec![None; size];
    let mut pi: Vec<Option<EstParent>> = vec![None; size];
    let mut pd: Vec<Option<EstParent>> = vec![None; size];
    let mut end = (0usize, 0usize, State::M, 0);
    for i in 0..=n {
        let mut window = IntronCandidateWindow::default();
        for j in 0..=m {
            let k = idx(i, j, cols);
            if j >= intron.min_len as usize {
                let start = j - intron.min_len as usize;
                let source_k = idx(i, start, cols);
                let (source, source_state) = best(&[
                    (mm[source_k], State::M),
                    (ii[source_k], State::I),
                    (dd[source_k], State::D),
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
                    if score > mm[k] {
                        let length = (j - candidate.start) as u32;
                        pm[k] = Some(EstParent {
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
                        mm[k] = score;
                    }
                }
            }
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let p = idx(i - 1, j - 1, cols);
                let (previous, previous_state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let score = add(
                    previous,
                    dna_score(query.bases[i - 1], bases[j - 1], Scoring::default()),
                );
                if score > mm[k] {
                    mm[k] = score;
                    pm[k] = Some(EstParent {
                        state: previous_state,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 1,
                            op: Op::Match,
                            query_advance: 1,
                            target_advance: 1,
                            repeats: 1,
                        }),
                        raw_fragment: RawFragment::with_prefix(
                            match previous_state {
                                State::I => Some(9),
                                State::D => Some(10),
                                _ => None,
                            },
                            &[RawStep {
                                transition_id: 1,
                                query_advance: 1,
                                target_advance: 1,
                                score: dna_score(
                                    query.bases[i - 1],
                                    bases[j - 1],
                                    Scoring::default(),
                                ),
                            }],
                        ),
                    });
                }
            }
            if i > 0 {
                let p = idx(i - 1, j, cols);
                let (score, previous_state, transition_id) = {
                    let from_m = add(mm[p], -12);
                    let from_i = add(ii[p], -4);
                    if from_m >= from_i {
                        (from_m, State::M, 2)
                    } else {
                        (from_i, State::I, 4)
                    }
                };
                if score >= 0 {
                    ii[k] = score;
                    pi[k] = Some(EstParent {
                        state: previous_state,
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
            if j > 0 {
                let p = idx(i, j - 1, cols);
                let (score, previous_state, transition_id) = {
                    let from_m = add(mm[p], -12);
                    let from_d = add(dd[p], -4);
                    if from_m >= from_d {
                        (from_m, State::M, 3)
                    } else {
                        (from_d, State::D, 5)
                    }
                };
                if score >= 0 {
                    dd[k] = score;
                    pd[k] = Some(EstParent {
                        state: previous_state,
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
            for (score, state) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if (score, i, j, state.rank()) > (end.3, end.0, end.1, end.2.rank()) {
                    end = (i, j, state, score);
                }
            }
        }
    }
    let (mut i, mut j, mut state, score) = end;
    let final_state = state;
    let (query_end, target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop {
        let parent = match state {
            State::M => pm[idx(i, j, cols)],
            State::I => pi[idx(i, j, cols)],
            State::D => pd[idx(i, j, cols)],
            State::Stop => None,
        };
        let Some(parent) = parent else { break };
        let (query_advance, target_advance) = parent.fragment.advances();
        i -= query_advance;
        j -= target_advance;
        fragments.push((parent.fragment, parent.raw_fragment));
        state = parent.state;
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
        (j as u64, target_end)
    } else {
        (
            target.bases.len() as u64 - j as u64,
            target.bases.len() as u64 - target_end,
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
        score,
        raw_trace,
        trace,
    }
}

/// Pair-disjoint suboptimal EST-to-genome alignments. Each target orientation
/// owns an independent forbidden-cell set, matching upstream strand semantics.
pub fn align_est2genome_suboptimal(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_est2genome_one(
            query,
            target,
            &target.bases,
            intron,
            Strand::Forward,
            Some(forbidden),
        )
    });
    if both_strands {
        let reverse_bases = reverse_complement(&target.bases);
        out.extend(enumerate_suboptimal_pair(
            query,
            target,
            threshold,
            |forbidden| {
                align_est2genome_one(
                    query,
                    target,
                    &reverse_bases,
                    intron,
                    Strand::Reverse,
                    Some(forbidden),
                )
            },
        ));
    }
    out
}

pub fn align_est2genome_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_est2genome_suboptimal(
                query,
                target,
                intron,
                threshold,
                both_strands,
            ));
        }
    }
    out
}

/// Score-only compatibility wrapper for callers that do not need a traceback.
pub fn est2genome_score(query: &Sequence, target: &Sequence, intron: IntronScoring) -> Score {
    align_est2genome(query, target, intron).score
}

fn align_model_ir_suboptimal_pair(
    query: &Sequence,
    target: &Sequence,
    ir: &model::ModelIr,
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    threshold: Score,
) -> Vec<Alignment> {
    enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_model_ir_with_intron_forbidden(
            query,
            target,
            ir,
            scoring,
            intron,
            strand,
            Some(forbidden),
        )
        .expect("validated built-in model must execute")
    })
}

pub fn align_ungapped_translated_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let ir = model::ungapped_translated(model::Scope::Anywhere);
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_model_ir_suboptimal_pair(
                query,
                target,
                &ir,
                scoring,
                IntronScoring::default(),
                Strand::Forward,
                threshold,
            ));
            if both_strands {
                let reverse_query = Sequence {
                    id: query.id.clone(),
                    bases: reverse_complement(&query.bases),
                };
                let mut reverse = align_model_ir_suboptimal_pair(
                    &reverse_query,
                    target,
                    &ir,
                    scoring,
                    IntronScoring::default(),
                    Strand::Forward,
                    threshold,
                );
                for alignment in &mut reverse {
                    alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                    alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                    alignment.query_strand = Strand::Reverse;
                }
                out.extend(reverse);
            }
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub fn align_ner_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    min_len: u32,
    max_len: u32,
    open: Score,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let ir = model::ner(min_len, max_len, open);
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_model_ir_suboptimal_pair(
                query,
                target,
                &ir,
                scoring,
                IntronScoring::default(),
                Strand::Forward,
                threshold,
            ));
            if both_strands {
                out.extend(align_model_ir_suboptimal_pair(
                    query,
                    target,
                    &ir,
                    scoring,
                    IntronScoring::default(),
                    Strand::Reverse,
                    threshold,
                ));
            }
        }
    }
    out
}

pub fn align_ungapped_translated_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let ir = model::ungapped_translated(model::Scope::Anywhere);
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(
                align_model_ir(query, target, &ir, scoring, Strand::Forward)
                    .expect("translated ungapped model must execute"),
            );
            if both_strands {
                let reverse_query = Sequence {
                    id: query.id.clone(),
                    bases: reverse_complement(&query.bases),
                };
                let mut reverse =
                    align_model_ir(&reverse_query, target, &ir, scoring, Strand::Forward)
                        .expect("reverse translated ungapped model must execute");
                let start = query.bases.len() as u64 - reverse.query_start;
                let end = query.bases.len() as u64 - reverse.query_end;
                reverse.query_start = start;
                reverse.query_end = end;
                reverse.query_strand = Strand::Reverse;
                out.push(reverse);
            }
        }
    }
    out
}

pub fn align_ner_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    min_len: u32,
    max_len: u32,
    open: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let ir = model::ner(min_len, max_len, open);
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(
                align_model_ir(query, target, &ir, scoring, Strand::Forward)
                    .expect("built-in NER model must execute"),
            );
            if both_strands {
                out.push(
                    align_model_ir(query, target, &ir, scoring, Strand::Reverse)
                        .expect("built-in reverse NER model must execute"),
                );
            }
        }
    }
    out
}

pub fn align_est2genome_database(
    queries: &[Sequence],
    targets: &[Sequence],
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets
                .iter()
                .map(move |target| align_est2genome_stranded(query, target, intron, both_strands))
        })
        .collect()
}

pub fn align_coding2coding_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets
                .iter()
                .map(move |target| align_coding2coding(query, target, scoring))
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
pub struct HeuristicConfig {
    pub word_len: usize,
    pub padding: usize,
    pub max_word_occurrences: usize,
}
impl Default for HeuristicConfig {
    fn default() -> Self {
        Self {
            word_len: 12,
            padding: 256,
            max_word_occurrences: 128,
        }
    }
}

#[derive(Clone, Copy)]
struct SeedCluster {
    count: usize,
    query_min: usize,
    query_max: usize,
    target_min: usize,
    target_max: usize,
}

fn encoded_kmer(bases: &[u8]) -> Option<u64> {
    let mut value = 0_u64;
    for &base in bases {
        value = (value << 2)
            | match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' | b'U' => 3,
                _ => return None,
            };
    }
    Some(value)
}

fn heuristic_region(
    query: &[u8],
    target: &[u8],
    config: HeuristicConfig,
) -> Option<(usize, usize, usize, usize)> {
    let k = config.word_len;
    if k == 0 || k > 31 || query.len() < k || target.len() < k {
        return None;
    }
    let mut index: HashMap<u64, Vec<usize>> = HashMap::new();
    for position in 0..=target.len() - k {
        if let Some(word) = encoded_kmer(&target[position..position + k]) {
            let positions = index.entry(word).or_default();
            if positions.len() < config.max_word_occurrences {
                positions.push(position);
            }
        }
    }
    let mut clusters: HashMap<i64, SeedCluster> = HashMap::new();
    for query_position in 0..=query.len() - k {
        let Some(word) = encoded_kmer(&query[query_position..query_position + k]) else {
            continue;
        };
        let Some(target_positions) = index.get(&word) else {
            continue;
        };
        for &target_position in target_positions {
            let diagonal = (query_position as i64 - target_position as i64).div_euclid(8);
            clusters
                .entry(diagonal)
                .and_modify(|cluster| {
                    cluster.count += 1;
                    cluster.query_min = cluster.query_min.min(query_position);
                    cluster.query_max = cluster.query_max.max(query_position + k);
                    cluster.target_min = cluster.target_min.min(target_position);
                    cluster.target_max = cluster.target_max.max(target_position + k);
                })
                .or_insert(SeedCluster {
                    count: 1,
                    query_min: query_position,
                    query_max: query_position + k,
                    target_min: target_position,
                    target_max: target_position + k,
                });
        }
    }
    let cluster = clusters.into_values().max_by_key(|cluster| {
        (
            cluster.count,
            cluster.query_max - cluster.query_min,
            cluster.target_max - cluster.target_min,
        )
    })?;
    Some((
        cluster.query_min.saturating_sub(config.padding),
        (cluster.query_max + config.padding).min(query.len()),
        cluster.target_min.saturating_sub(config.padding),
        (cluster.target_max + config.padding).min(target.len()),
    ))
}

fn protein_to_dna_target_region(
    query: &[u8],
    target: &[u8],
    config: HeuristicConfig,
    max_intron: usize,
) -> Option<(usize, usize)> {
    let k = config.word_len;
    if k == 0 || query.len() < k || target.len() < 3 * k {
        return None;
    }
    let mut index: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    for frame in 0..3 {
        let translated = translate_dna(target, frame);
        if translated.len() < k {
            continue;
        }
        for position in 0..=translated.len() - k {
            let word = translated[position..position + k].to_vec();
            let positions = index.entry(word).or_default();
            if positions.len() < config.max_word_occurrences {
                positions.push(frame as usize + 3 * position);
            }
        }
    }
    let mut clusters: HashMap<i64, SeedCluster> = HashMap::new();
    for query_position in 0..=query.len() - k {
        let word = query[query_position..query_position + k]
            .iter()
            .map(u8::to_ascii_uppercase)
            .collect::<Vec<_>>();
        let Some(target_positions) = index.get(&word) else {
            continue;
        };
        for &target_position in target_positions {
            let diagonal = (3 * query_position as i64 - target_position as i64).div_euclid(24);
            clusters
                .entry(diagonal)
                .and_modify(|cluster| {
                    cluster.count += 1;
                    cluster.query_min = cluster.query_min.min(query_position);
                    cluster.query_max = cluster.query_max.max(query_position + k);
                    cluster.target_min = cluster.target_min.min(target_position);
                    cluster.target_max = cluster.target_max.max(target_position + 3 * k);
                })
                .or_insert(SeedCluster {
                    count: 1,
                    query_min: query_position,
                    query_max: query_position + k,
                    target_min: target_position,
                    target_max: target_position + 3 * k,
                });
        }
    }
    let cluster = clusters.into_values().max_by_key(|cluster| {
        (
            cluster.count,
            cluster.query_max - cluster.query_min,
            cluster.target_max - cluster.target_min,
        )
    })?;
    let left_context = 3 * cluster.query_min + config.padding + max_intron;
    let right_context =
        3 * query.len().saturating_sub(cluster.query_max) + config.padding + max_intron;
    Some((
        cluster.target_min.saturating_sub(left_context),
        cluster
            .target_max
            .saturating_add(right_context)
            .min(target.len()),
    ))
}

fn spliced_target_region(
    query: &[u8],
    target: &[u8],
    config: HeuristicConfig,
    max_intron: usize,
) -> Option<(usize, usize)> {
    let (query_start, query_end, target_start, target_end) =
        heuristic_region(query, target, config)?;
    let left_context = query_start.saturating_add(max_intron);
    let right_context = query
        .len()
        .saturating_sub(query_end)
        .saturating_add(max_intron);
    Some((
        target_start.saturating_sub(left_context),
        target_end.saturating_add(right_context).min(target.len()),
    ))
}

fn align_est2genome_heuristic_pair(
    query: &Sequence,
    target: &Sequence,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Alignment {
    let Some((target_start, target_end)) =
        spliced_target_region(&query.bases, &target.bases, config, intron.max_len as usize)
    else {
        return align_est2genome_stranded(query, target, intron, both_strands);
    };
    let target_slice = Sequence {
        id: target.id.clone(),
        bases: target.bases[target_start..target_end].to_vec(),
    };
    let mut alignment = align_est2genome_stranded(query, &target_slice, intron, both_strands);
    remap_target_slice(&mut alignment, target_start, target.bases.len());
    alignment
}

fn remap_target_slice(alignment: &mut Alignment, target_start: usize, target_len: usize) {
    alignment.target_start += target_start as u64;
    alignment.target_end += target_start as u64;
    alignment.target_len = target_len as u64;
}

fn align_protein_to_genome_heuristic_oriented(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    config: HeuristicConfig,
    strand: Strand,
    bestfit: bool,
) -> Alignment {
    let oriented_target = if strand == Strand::Reverse {
        reverse_complement(&target.bases)
    } else {
        target.bases.clone()
    };
    let Some((oriented_start, oriented_end)) = protein_to_dna_target_region(
        &query.bases,
        &oriented_target,
        config,
        intron.max_len as usize,
    ) else {
        return protein2genome_alignment_oriented(
            query,
            target,
            &oriented_target,
            scoring,
            intron,
            strand,
            bestfit,
            None,
        );
    };
    let (target_start, target_end) = if strand == Strand::Reverse {
        (
            target.bases.len() - oriented_end,
            target.bases.len() - oriented_start,
        )
    } else {
        (oriented_start, oriented_end)
    };
    let target_slice = Sequence {
        id: target.id.clone(),
        bases: target.bases[target_start..target_end].to_vec(),
    };
    let oriented_slice = if strand == Strand::Reverse {
        reverse_complement(&target_slice.bases)
    } else {
        target_slice.bases.clone()
    };
    let mut alignment = protein2genome_alignment_oriented(
        query,
        &target_slice,
        &oriented_slice,
        scoring,
        intron,
        strand,
        bestfit,
        None,
    );
    remap_target_slice(&mut alignment, target_start, target.bases.len());
    alignment
}

fn align_protein_to_genome_database_heuristic_impl(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
    bestfit: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_protein_to_genome_heuristic_oriented(
                query,
                target,
                scoring,
                intron,
                config,
                Strand::Forward,
                bestfit,
            ));
            if both_strands {
                out.push(align_protein_to_genome_heuristic_oriented(
                    query,
                    target,
                    scoring,
                    intron,
                    config,
                    Strand::Reverse,
                    bestfit,
                ));
            }
        }
    }
    out
}

pub fn align_protein_to_genome_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    align_protein_to_genome_database_heuristic_impl(
        queries,
        targets,
        scoring,
        intron,
        both_strands,
        config,
        false,
    )
}

pub fn align_protein_to_genome_bestfit_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    align_protein_to_genome_database_heuristic_impl(
        queries,
        targets,
        scoring,
        intron,
        both_strands,
        config,
        true,
    )
}

pub fn align_est2genome_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets.iter().map(move |target| {
                align_est2genome_heuristic_pair(query, target, intron, both_strands, config)
            })
        })
        .collect()
}

pub fn align_cdna_to_genome_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            let Some((target_start, target_end)) =
                spliced_target_region(&query.bases, &target.bases, config, intron.max_len as usize)
            else {
                out.push(align_cdna_to_genome(query, target, scoring, intron));
                continue;
            };
            let target_slice = Sequence {
                id: target.id.clone(),
                bases: target.bases[target_start..target_end].to_vec(),
            };
            let mut alignment = align_cdna_to_genome(query, &target_slice, scoring, intron);
            remap_target_slice(&mut alignment, target_start, target.bases.len());
            out.push(alignment);
        }
    }
    out
}

pub fn align_genome_to_genome_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    align_genome_to_genome_database_heuristic_with_rolling_scores(
        queries, targets, scoring, intron, config, false,
    )
}

fn align_genome_to_genome_database_heuristic_with_rolling_scores(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    config: HeuristicConfig,
    rolling_scores: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            let Some((target_start, target_end)) =
                spliced_target_region(&query.bases, &target.bases, config, intron.max_len as usize)
            else {
                out.push(align_genome_to_genome_memory_aware(
                    query,
                    target,
                    scoring,
                    intron,
                    rolling_scores,
                    None,
                ));
                continue;
            };
            let target_slice = Sequence {
                id: target.id.clone(),
                bases: target.bases[target_start..target_end].to_vec(),
            };
            let mut alignment = align_genome_to_genome_memory_aware(
                query,
                &target_slice,
                scoring,
                intron,
                rolling_scores,
                None,
            );
            remap_target_slice(&mut alignment, target_start, target.bases.len());
            out.push(alignment);
        }
    }
    out
}

/// Heuristic genome-to-genome search over the same four orientation pairs as
/// exhaustive search.  Each oriented pair keeps its own seed refinement.
pub fn align_genome_to_genome_database_heuristic_stranded(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    align_genome_to_genome_database_heuristic_stranded_with_rolling_scores(
        queries,
        targets,
        scoring,
        intron,
        both_strands,
        config,
        false,
    )
}

/// Genome-to-genome heuristic search with rolling score rows.  This is used
/// by the CLI's forced low-memory mode while checkpoint parent replay is being
/// completed.
pub fn align_genome_to_genome_database_heuristic_stranded_with_rolling_scores(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
    rolling_scores: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(
                align_genome_to_genome_database_heuristic_with_rolling_scores(
                    std::slice::from_ref(query),
                    std::slice::from_ref(target),
                    scoring,
                    intron,
                    config,
                    rolling_scores,
                ),
            );
            if !both_strands {
                continue;
            }
            let reverse_query = Sequence {
                id: query.id.clone(),
                bases: reverse_complement(&query.bases),
            };
            let reverse_target = Sequence {
                id: target.id.clone(),
                bases: reverse_complement(&target.bases),
            };
            for (oriented_query, query_reverse) in [(query, false), (&reverse_query, true)] {
                for (oriented_target, target_reverse) in [(target, false), (&reverse_target, true)]
                {
                    if !query_reverse && !target_reverse {
                        continue;
                    }
                    let mut alignment =
                        align_genome_to_genome_database_heuristic_with_rolling_scores(
                            std::slice::from_ref(oriented_query),
                            std::slice::from_ref(oriented_target),
                            scoring,
                            intron,
                            config,
                            rolling_scores,
                        )
                        .pop()
                        .expect("one oriented genome pair produces one alignment");
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

pub fn align_coding_to_genome_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            let Some((target_start, target_end)) =
                spliced_target_region(&query.bases, &target.bases, config, intron.max_len as usize)
            else {
                out.extend(align_coding_to_genome(
                    query,
                    target,
                    scoring,
                    intron,
                    both_strands,
                ));
                continue;
            };
            let target_slice = Sequence {
                id: target.id.clone(),
                bases: target.bases[target_start..target_end].to_vec(),
            };
            let mut alignments =
                align_coding_to_genome(query, &target_slice, scoring, intron, both_strands);
            for alignment in &mut alignments {
                remap_target_slice(alignment, target_start, target.bases.len());
            }
            out.extend(alignments);
        }
    }
    out
}

fn align_heuristic_pair(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    config: HeuristicConfig,
) -> Alignment {
    let Some((query_start, query_end, target_start, target_end)) =
        heuristic_region(&query.bases, &target.bases, config)
    else {
        return align(query, target, model, scoring, Strand::Forward);
    };
    let query_slice = Sequence {
        id: query.id.clone(),
        bases: query.bases[query_start..query_end].to_vec(),
    };
    let target_slice = Sequence {
        id: target.id.clone(),
        bases: target.bases[target_start..target_end].to_vec(),
    };
    let mut alignment = align(&query_slice, &target_slice, model, scoring, Strand::Forward);
    alignment.query_start += query_start as u64;
    alignment.query_end += query_start as u64;
    alignment.target_start += target_start as u64;
    alignment.target_end += target_start as u64;
    alignment.target_len = target.bases.len() as u64;
    alignment
}

pub fn align_database_heuristic(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    both_strands: bool,
    config: HeuristicConfig,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_heuristic_pair(query, target, model, scoring, config));
            if both_strands {
                let reverse_query = Sequence {
                    id: query.id.clone(),
                    bases: reverse_complement(&query.bases),
                };
                let mut alignment =
                    align_heuristic_pair(&reverse_query, target, model, scoring, config);
                alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                alignment.query_strand = Strand::Reverse;
                out.push(alignment);
            }
        }
    }
    out
}

pub fn align_database(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for q in queries {
        for t in targets {
            out.push(align(q, t, model, scoring, Strand::Forward));
            if both_strands {
                let reverse_query = Sequence {
                    id: q.id.clone(),
                    bases: reverse_complement(&q.bases),
                };
                let mut reverse = align(&reverse_query, t, model, scoring, Strand::Forward);
                reverse.query_start = q.bases.len() as u64 - reverse.query_start;
                reverse.query_end = q.bases.len() as u64 - reverse.query_end;
                reverse.query_strand = Strand::Reverse;
                out.push(reverse);
            }
        }
    }
    out
}

fn forbid_alignment_pairs(
    alignment: &Alignment,
    query_len: usize,
    target_len: usize,
    forbidden: &mut HashSet<(usize, usize)>,
) -> usize {
    let mut query_position = if alignment.query_strand == Strand::Reverse {
        query_len.saturating_sub(alignment.query_start as usize)
    } else {
        alignment.query_start as usize
    };
    let mut target_position = if alignment.target_strand == Strand::Reverse {
        target_len.saturating_sub(alignment.target_start as usize)
    } else {
        alignment.target_start as usize
    };
    let before = forbidden.len();
    for run in &alignment.trace {
        for _ in 0..run.repeats {
            if run.op == Op::Match
                || (run.op == Op::SplitCodon && run.query_advance > 0 && run.target_advance > 0)
            {
                let mut divisor = run.query_advance.min(run.target_advance) as usize;
                let mut remainder = run.query_advance.max(run.target_advance) as usize;
                while divisor > 0 {
                    let next = remainder % divisor;
                    remainder = divisor;
                    divisor = next;
                }
                let steps = remainder.max(1);
                let query_step = run.query_advance as usize / steps;
                let target_step = run.target_advance as usize / steps;
                for step in 0..steps {
                    forbidden.insert((
                        query_position + step * query_step,
                        target_position + step * target_step,
                    ));
                }
            }
            query_position += run.query_advance as usize;
            target_position += run.target_advance as usize;
        }
    }
    forbidden.len() - before
}

fn forbidden_match_transition(
    forbidden: Option<&HashSet<(usize, usize)>>,
    query_start: usize,
    target_start: usize,
    query_advance: usize,
    target_advance: usize,
) -> bool {
    let Some(forbidden) = forbidden else {
        return false;
    };
    if query_advance == 0 || target_advance == 0 {
        return false;
    }
    let mut divisor = query_advance.min(target_advance);
    let mut remainder = query_advance.max(target_advance);
    while divisor > 0 {
        let next = remainder % divisor;
        remainder = divisor;
        divisor = next;
    }
    let steps = remainder.max(1);
    let query_step = query_advance / steps;
    let target_step = target_advance / steps;
    (0..steps).any(|step| {
        forbidden.contains(&(
            query_start + step * query_step,
            target_start + step * target_step,
        ))
    })
}

fn forbidden_split_codon_transition(
    forbidden: Option<&HashSet<(usize, usize)>>,
    query_start: usize,
    pre_target_start: usize,
    pre_len: usize,
    post_target_start: usize,
    post_len: usize,
) -> bool {
    let Some(forbidden) = forbidden else {
        return false;
    };
    (0..pre_len)
        .any(|offset| forbidden.contains(&(query_start + offset, pre_target_start + offset)))
        || (0..post_len).any(|offset| {
            forbidden.contains(&(query_start + pre_len + offset, post_target_start + offset))
        })
}

fn align_suboptimal_pair(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut forbidden = HashSet::new();
    let mut out = Vec::new();
    loop {
        let alignment = align_with_scorer_forbidden(
            query,
            target,
            Model::Local,
            scoring,
            Strand::Forward,
            dna_score,
            Some(&forbidden),
        );
        if alignment.score < threshold || alignment.score <= 0 {
            break;
        }
        let added = forbid_alignment_pairs(
            &alignment,
            query.bases.len(),
            target.bases.len(),
            &mut forbidden,
        );
        out.push(alignment);
        if added == 0 {
            break;
        }
    }
    out
}

pub fn align_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_suboptimal_pair(query, target, scoring, threshold));
            if both_strands {
                let reverse_query = Sequence {
                    id: query.id.clone(),
                    bases: reverse_complement(&query.bases),
                };
                for mut alignment in
                    align_suboptimal_pair(&reverse_query, target, scoring, threshold)
                {
                    alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                    alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                    alignment.query_strand = Strand::Reverse;
                    out.push(alignment);
                }
            }
        }
    }
    out
}

pub fn align_protein_database(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.push(align_protein(query, target, model, scoring));
        }
    }
    out
}

/// Waterman--Eggert local protein alignments using the same pair exclusion
/// rule as DNA affine alignment.
pub fn align_protein_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut forbidden = HashSet::new();
    let mut out = Vec::new();
    loop {
        let mut alignment = align_with_scorer_forbidden(
            query,
            target,
            Model::Local,
            scoring,
            Strand::Forward,
            protein_score,
            Some(&forbidden),
        );
        alignment.query_strand = Strand::Unknown;
        alignment.target_strand = Strand::Unknown;
        if alignment.score < threshold || alignment.score <= 0 {
            break;
        }
        let added = forbid_alignment_pairs(
            &alignment,
            query.bases.len(),
            target.bases.len(),
            &mut forbidden,
        );
        out.push(alignment);
        if added == 0 {
            break;
        }
    }
    out
}

pub fn align_protein_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_suboptimal(query, target, scoring, threshold));
        }
    }
    out
}

/// Exhaustive local coding-DNA alignment with codon-affine gaps and frameshifts.
fn align_coding2coding_forbidden(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let (n, m, cols) = (
        query.bases.len(),
        target.bases.len(),
        target.bases.len() + 1,
    );
    let size = (n + 1) * (m + 1);
    let mut mm = vec![0; size];
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut pm: Vec<Option<P2Parent>> = vec![None; size];
    let mut pi: Vec<Option<P2Parent>> = vec![None; size];
    let mut pd: Vec<Option<P2Parent>> = vec![None; size];
    let mut end = (0usize, 0usize, State::M, 0);
    for i in 0..=n {
        for j in 0..=m {
            let k = idx(i, j, cols);
            if i >= 3 {
                let p = idx(i - 3, j, cols);
                let (from_m, from_i) = (
                    add(mm[p], scoring.codon_gap_open),
                    add(ii[p], scoring.codon_gap_extend),
                );
                if from_m >= from_i && from_m >= 0 {
                    ii[k] = from_m;
                    pi[k] = Some(P2Parent {
                        state: State::M,
                        op: Op::Insert,
                        transition_id: 2,
                        query_advance: 3,
                        target_advance: 0,
                    });
                } else if from_i >= 0 {
                    ii[k] = from_i;
                    pi[k] = Some(P2Parent {
                        state: State::I,
                        op: Op::Insert,
                        transition_id: 4,
                        query_advance: 3,
                        target_advance: 0,
                    });
                }
            }
            if j >= 3 {
                let p = idx(i, j - 3, cols);
                let (from_m, from_d) = (
                    add(mm[p], scoring.codon_gap_open),
                    add(dd[p], scoring.codon_gap_extend),
                );
                if from_m >= from_d && from_m >= 0 {
                    dd[k] = from_m;
                    pd[k] = Some(P2Parent {
                        state: State::M,
                        op: Op::Delete,
                        transition_id: 3,
                        query_advance: 0,
                        target_advance: 3,
                    });
                } else if from_d >= 0 {
                    dd[k] = from_d;
                    pd[k] = Some(P2Parent {
                        state: State::D,
                        op: Op::Delete,
                        transition_id: 5,
                        query_advance: 0,
                        target_advance: 3,
                    });
                }
            }
            for advance in [1_usize, 2, 4, 5] {
                if i >= advance {
                    let p = idx(i - advance, j, cols);
                    let (value, state) =
                        best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                    let value = add(value, scoring.frameshift);
                    if value > mm[k] {
                        mm[k] = value;
                        pm[k] = Some(P2Parent {
                            state,
                            op: Op::Frameshift,
                            transition_id: 6,
                            query_advance: advance as u32,
                            target_advance: 0,
                        });
                    }
                }
                if j >= advance {
                    let p = idx(i, j - advance, cols);
                    let (value, state) =
                        best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                    let value = add(value, scoring.frameshift);
                    if value > mm[k] {
                        mm[k] = value;
                        pm[k] = Some(P2Parent {
                            state,
                            op: Op::Frameshift,
                            transition_id: 7,
                            query_advance: 0,
                            target_advance: advance as u32,
                        });
                    }
                }
            }
            if i >= 3 && j >= 3 && !forbidden_match_transition(forbidden, i - 3, j - 3, 3, 3) {
                let p = idx(i - 3, j - 3, cols);
                let (value, state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let qa = translate_dna(&query.bases[i - 3..i], 0)[0];
                let ta = translate_dna(&target.bases[j - 3..j], 0)[0];
                let value = add(value, protein_score(qa, ta, scoring));
                if value > mm[k] {
                    mm[k] = value;
                    pm[k] = Some(P2Parent {
                        state,
                        op: Op::Match,
                        transition_id: 1,
                        query_advance: 3,
                        target_advance: 3,
                    });
                }
            }
            for (value, state) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if (value, i, j, state.rank()) > (end.3, end.0, end.1, end.2.rank()) {
                    end = (i, j, state, value);
                }
            }
        }
    }
    let (mut i, mut j, mut state, score) = end;
    let final_state = state;
    let (query_end, target_end) = (i as u64, j as u64);
    let mut trace = Vec::new();
    let mut reversed_raw_fragments: Vec<Vec<RawStep>> = Vec::new();
    while state != State::Stop {
        let parent = match state {
            State::M => pm[idx(i, j, cols)],
            State::I => pi[idx(i, j, cols)],
            State::D => pd[idx(i, j, cols)],
            State::Stop => None,
        };
        let Some(parent) = parent else { break };
        let mut raw_fragment = Vec::new();
        if state == State::M {
            let epsilon = match parent.state {
                State::I => Some(4),
                State::D => Some(7),
                _ => None,
            };
            if let Some(transition_id) = epsilon {
                raw_fragment.push(RawStep {
                    transition_id,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                });
            }
        }
        match parent.op {
            Op::Match => raw_fragment.push(RawStep {
                transition_id: 1,
                query_advance: 3,
                target_advance: 3,
                score: protein_score(
                    translate_dna(&query.bases[i - 3..i], 0)[0],
                    translate_dna(&target.bases[j - 3..j], 0)[0],
                    scoring,
                ),
            }),
            Op::Insert => raw_fragment.push(RawStep {
                transition_id: if parent.transition_id == 2 { 2 } else { 3 },
                query_advance: 3,
                target_advance: 0,
                score: if parent.transition_id == 2 {
                    scoring.codon_gap_open
                } else {
                    scoring.codon_gap_extend
                },
            }),
            Op::Delete => raw_fragment.push(RawStep {
                transition_id: if parent.transition_id == 3 { 5 } else { 6 },
                query_advance: 0,
                target_advance: 3,
                score: if parent.transition_id == 3 {
                    scoring.codon_gap_open
                } else {
                    scoring.codon_gap_extend
                },
            }),
            Op::Frameshift => {
                if parent.query_advance > 0 {
                    let short = parent.query_advance % 3;
                    raw_fragment.push(RawStep {
                        transition_id: if short == 1 { 8 } else { 9 },
                        query_advance: short,
                        target_advance: 0,
                        score: scoring.frameshift,
                    });
                    raw_fragment.push(RawStep {
                        transition_id: if parent.query_advance > 3 { 11 } else { 10 },
                        query_advance: if parent.query_advance > 3 { 3 } else { 0 },
                        target_advance: 0,
                        score: 0,
                    });
                } else {
                    let short = parent.target_advance % 3;
                    raw_fragment.push(RawStep {
                        transition_id: if short == 1 { 12 } else { 13 },
                        query_advance: 0,
                        target_advance: short,
                        score: scoring.frameshift,
                    });
                    raw_fragment.push(RawStep {
                        transition_id: if parent.target_advance > 3 { 15 } else { 14 },
                        query_advance: 0,
                        target_advance: if parent.target_advance > 3 { 3 } else { 0 },
                        score: 0,
                    });
                }
            }
            _ => unreachable!("coding2coding specialized traceback operation"),
        }
        reversed_raw_fragments.push(raw_fragment);
        i -= parent.query_advance as usize;
        j -= parent.target_advance as usize;
        trace.push(TraceRun {
            transition_id: parent.transition_id,
            op: parent.op,
            query_advance: parent.query_advance,
            target_advance: parent.target_advance,
            repeats: 1,
        });
        state = parent.state;
    }
    trace.reverse();
    reversed_raw_fragments.reverse();
    let mut raw_trace = vec![RawStep {
        transition_id: 0,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for fragment in reversed_raw_fragments {
        raw_trace.extend(fragment);
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
    let mut canonical = Vec::new();
    for run in trace {
        TraceFragment::one(run).append_to(&mut canonical);
    }
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
        score,
        raw_trace,
        trace: canonical,
    }
}

/// Exhaustive local coding-DNA alignment with codon-affine gaps and frameshifts.
pub fn align_coding2coding(query: &Sequence, target: &Sequence, scoring: Scoring) -> Alignment {
    align_coding2coding_forbidden(query, target, scoring, None)
}

fn enumerate_suboptimal_pair<F>(
    query: &Sequence,
    target: &Sequence,
    threshold: Score,
    mut aligner: F,
) -> Vec<Alignment>
where
    F: FnMut(&HashSet<(usize, usize)>) -> Alignment,
{
    let mut forbidden = HashSet::new();
    let mut out = Vec::new();
    loop {
        let alignment = aligner(&forbidden);
        if alignment.score < threshold || alignment.score <= 0 {
            break;
        }
        let added = forbid_alignment_pairs(
            &alignment,
            query.bases.len(),
            target.bases.len(),
            &mut forbidden,
        );
        out.push(alignment);
        if added == 0 {
            break;
        }
    }
    out
}

/// Waterman-Eggert-style pair-disjoint suboptimal coding alignments.
pub fn align_coding2coding_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    threshold: Score,
) -> Vec<Alignment> {
    enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_coding2coding_forbidden(query, target, scoring, Some(forbidden))
    })
}

pub fn align_coding2coding_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_coding2coding_suboptimal(
                query, target, scoring, threshold,
            ));
        }
    }
    out
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CdnaState {
    U5M,
    U5I,
    U5D,
    CodingM,
    CodingI,
    CodingD,
    U3M,
    U3I,
    U3D,
}

#[derive(Clone, Copy)]
struct CdnaParent {
    state: CdnaState,
    fragment: TraceFragment,
}

fn cdna_relax(
    score: &mut Score,
    parent: &mut Option<CdnaParent>,
    candidate: Score,
    state: CdnaState,
    fragment: TraceFragment,
) {
    if candidate > *score || (candidate == *score && fragment.atoms.iter().any(Option::is_some)) {
        *score = candidate;
        *parent = Some(CdnaParent { state, fragment });
    }
}

struct CdnaDp {
    cols: usize,
    u5m: Vec<Score>,
    u5i: Vec<Score>,
    u5d: Vec<Score>,
    cm: Vec<Score>,
    ci: Vec<Score>,
    cd: Vec<Score>,
    u3m: Vec<Score>,
    u3i: Vec<Score>,
    u3d: Vec<Score>,
    pu5m: Vec<Option<CdnaParent>>,
    pu5i: Vec<Option<CdnaParent>>,
    pu5d: Vec<Option<CdnaParent>>,
    pcm: Vec<Option<CdnaParent>>,
    pci: Vec<Option<CdnaParent>>,
    pcd: Vec<Option<CdnaParent>>,
    pu3m: Vec<Option<CdnaParent>>,
    pu3i: Vec<Option<CdnaParent>>,
    pu3d: Vec<Option<CdnaParent>>,
    end: (usize, usize, CdnaState, Score),
}
impl CdnaDp {
    fn value(&self, state: CdnaState, k: usize) -> Score {
        match state {
            CdnaState::U5M => self.u5m[k],
            CdnaState::U5I => self.u5i[k],
            CdnaState::U5D => self.u5d[k],
            CdnaState::CodingM => self.cm[k],
            CdnaState::CodingI => self.ci[k],
            CdnaState::CodingD => self.cd[k],
            CdnaState::U3M => self.u3m[k],
            CdnaState::U3I => self.u3i[k],
            CdnaState::U3D => self.u3d[k],
        }
    }
    fn parent(&self, state: CdnaState, k: usize) -> Option<CdnaParent> {
        match state {
            CdnaState::U5M => self.pu5m[k],
            CdnaState::U5I => self.pu5i[k],
            CdnaState::U5D => self.pu5d[k],
            CdnaState::CodingM => self.pcm[k],
            CdnaState::CodingI => self.pci[k],
            CdnaState::CodingD => self.pcd[k],
            CdnaState::U3M => self.pu3m[k],
            CdnaState::U3I => self.pu3i[k],
            CdnaState::U3D => self.pu3d[k],
        }
    }
}

/// Score-only cDNA-to-genome composition: 5-prime UTR, coding region, and 3-prime UTR.
///
/// This executes the upstream graph in one DP lattice rather than concatenating
/// independent alignments. UTR states use DNA affine/intron transitions; coding
/// states use codon affine gaps, bilateral frameshifts, and target phased introns.
fn cdna2genome_dp(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> CdnaDp {
    let intron = canonical_intron_scoring(intron);
    let (n, m, cols) = (
        query.bases.len(),
        target.bases.len(),
        target.bases.len() + 1,
    );
    let size = (n + 1) * (m + 1);
    let mut u5m = vec![0; size];
    let mut pu5m: Vec<Option<CdnaParent>> = vec![None; size];
    let mut u5i = vec![NEG_INF; size];
    let mut pu5i: Vec<Option<CdnaParent>> = vec![None; size];
    let mut u5d = vec![NEG_INF; size];
    let mut pu5d: Vec<Option<CdnaParent>> = vec![None; size];
    let mut cm = vec![0; size];
    let mut pcm: Vec<Option<CdnaParent>> = vec![None; size];
    let mut ci = vec![NEG_INF; size];
    let mut pci: Vec<Option<CdnaParent>> = vec![None; size];
    let mut cd = vec![NEG_INF; size];
    let mut pcd: Vec<Option<CdnaParent>> = vec![None; size];
    let mut u3m = vec![NEG_INF; size];
    let mut pu3m: Vec<Option<CdnaParent>> = vec![None; size];
    let mut u3i = vec![NEG_INF; size];
    let mut pu3i: Vec<Option<CdnaParent>> = vec![None; size];
    let mut u3d = vec![NEG_INF; size];
    let mut pu3d: Vec<Option<CdnaParent>> = vec![None; size];
    let mut end = (0usize, 0usize, CdnaState::CodingM, 0);
    for i in 0..=n {
        let mut w5 = IntronCandidateWindow::default();
        let mut wc = [
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
        ];
        let mut w3 = IntronCandidateWindow::default();
        for j in 0..=m {
            let k = idx(i, j, cols);
            // 5-prime UTR: local EST-to-genome subgraph.
            if j >= intron.min_len as usize {
                let s = j - intron.min_len as usize;
                let p = idx(i, s, cols);
                let (v, st) = best(&[(u5m[p], State::M), (u5i[p], State::I), (u5d[p], State::D)]);
                if v > 0 {
                    if let Some(d) = splice_score(
                        &target.bases,
                        s,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        w5.insert(IntronCandidate {
                            start: s,
                            score: add(v, intron.open_penalty.saturating_add(d)),
                            state_rank: st.rank(),
                        });
                    }
                }
            }
            if j > intron.max_len as usize {
                w5.expire_before(j - intron.max_len as usize);
            }
            if j >= 2 {
                if let (Some(c), Some(a)) = (
                    w5.best(),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    cdna_relax(
                        &mut u5m[k],
                        &mut pu5m[k],
                        add(c.score, a),
                        match state_from_rank(c.state_rank) {
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
                                target_advance: (j - c.start - 4) as u32,
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
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let p = idx(i - 1, j - 1, cols);
                let (v, st) = best(&[(u5m[p], State::M), (u5i[p], State::I), (u5d[p], State::D)]);
                cdna_relax(
                    &mut u5m[k],
                    &mut pu5m[k],
                    add(
                        v,
                        dna_score(query.bases[i - 1], target.bases[j - 1], scoring),
                    ),
                    match st {
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
            if i > 0 {
                let p = idx(i - 1, j, cols);
                let (candidate, state, transition_id) =
                    if add(u5m[p], scoring.gap_open) >= add(u5i[p], scoring.gap_extend) {
                        (add(u5m[p], scoring.gap_open), CdnaState::U5M, 2)
                    } else {
                        (add(u5i[p], scoring.gap_extend), CdnaState::U5I, 4)
                    };
                cdna_relax(
                    &mut u5i[k],
                    &mut pu5i[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j > 0 {
                let p = idx(i, j - 1, cols);
                let (candidate, state, transition_id) =
                    if add(u5m[p], scoring.gap_open) >= add(u5d[p], scoring.gap_extend) {
                        (add(u5m[p], scoring.gap_open), CdnaState::U5M, 3)
                    } else {
                        (add(u5d[p], scoring.gap_extend), CdnaState::U5D, 5)
                    };
                cdna_relax(
                    &mut u5d[k],
                    &mut pu5d[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            // epsilon transition: UTR into the coding match state.
            let (candidate, state) =
                best(&[(u5m[k], State::M), (u5i[k], State::I), (u5d[k], State::D)]);
            cdna_relax(
                &mut cm[k],
                &mut pcm[k],
                candidate,
                match state {
                    State::M => CdnaState::U5M,
                    State::I => CdnaState::U5I,
                    State::D => CdnaState::U5D,
                    State::Stop => unreachable!(),
                },
                TraceFragment::empty(),
            );
            // CDS target phase-0/1/2 introns.
            if i >= 3 {
                for (phase, window) in wc.iter_mut().enumerate() {
                    let post = 3 - phase;
                    if j >= post + intron.min_len as usize + phase {
                        let post_start = j - post;
                        let donor = post_start - intron.min_len as usize;
                        let source_j = donor - phase;
                        let p = idx(i - 3, source_j, cols);
                        let (v, st) =
                            best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                        if v > 0 {
                            if let Some(d) = splice_score(
                                &target.bases,
                                donor,
                                SpliceType::DonorForward,
                                intron.force_gtag,
                            ) {
                                window.insert(IntronCandidate {
                                    start: donor,
                                    score: add(v, intron.open_penalty.saturating_add(d)),
                                    state_rank: st.rank(),
                                });
                            }
                        }
                        if post_start > intron.max_len as usize {
                            window.expire_before(post_start - intron.max_len as usize);
                        }
                        if post_start >= 2 {
                            if let (Some(c), Some(a)) = (
                                window.best(),
                                splice_score(
                                    &target.bases,
                                    post_start - 2,
                                    SpliceType::AcceptorForward,
                                    intron.force_gtag,
                                ),
                            ) {
                                let pre = c.start - phase;
                                if forbidden_split_codon_transition(
                                    forbidden,
                                    i - 3,
                                    pre,
                                    phase,
                                    post_start,
                                    post,
                                ) {
                                    continue;
                                }
                                let mut codon = Vec::with_capacity(3);
                                codon.extend_from_slice(&target.bases[pre..c.start]);
                                codon.extend_from_slice(&target.bases[post_start..j]);
                                if codon.len() == 3 {
                                    let qaa = translate_dna(&query.bases[i - 3..i], 0)[0];
                                    let taa = translate_dna(&codon, 0)[0];
                                    let pre_q = phase as u32;
                                    let post_q = (3 - phase) as u32;
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
                                                target_advance: (post_start - c.start - 4) as u32,
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
                                                query_advance: pre_q,
                                                target_advance: pre_q,
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
                                                target_advance: (post_start - c.start - 4) as u32,
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
                                                query_advance: post_q,
                                                target_advance: post_q,
                                                repeats: 1,
                                            },
                                        )
                                    };
                                    cdna_relax(
                                        &mut cm[k],
                                        &mut pcm[k],
                                        add(add(c.score, a), protein_score(qaa, taa, scoring)),
                                        match state_from_rank(c.state_rank) {
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
                    let p = idx(i, j - 3, cols);
                    let (candidate, state, transition_id) = if add(cm[p], scoring.codon_gap_open)
                        >= add(cd[p], scoring.codon_gap_extend)
                    {
                        (add(cm[p], scoring.codon_gap_open), CdnaState::CodingM, 3)
                    } else {
                        (add(cd[p], scoring.codon_gap_extend), CdnaState::CodingD, 5)
                    };
                    cdna_relax(
                        &mut cd[k],
                        &mut pcd[k],
                        candidate,
                        state,
                        TraceFragment::one(TraceRun {
                            transition_id,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    );
                }
                let p = idx(i - 3, j, cols);
                let (candidate, state, transition_id) =
                    if add(cm[p], scoring.codon_gap_open) >= add(ci[p], scoring.codon_gap_extend) {
                        (add(cm[p], scoring.codon_gap_open), CdnaState::CodingM, 2)
                    } else {
                        (add(ci[p], scoring.codon_gap_extend), CdnaState::CodingI, 4)
                    };
                cdna_relax(
                    &mut ci[k],
                    &mut pci[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 3,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
                if j >= 3 && !forbidden_match_transition(forbidden, i - 3, j - 3, 3, 3) {
                    let p = idx(i - 3, j - 3, cols);
                    let (v, st) = best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                    cdna_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(
                            v,
                            protein_score(
                                translate_dna(&query.bases[i - 3..i], 0)[0],
                                translate_dna(&target.bases[j - 3..j], 0)[0],
                                scoring,
                            ),
                        ),
                        match st {
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
            for a in [1_usize, 2, 4, 5] {
                if i >= a {
                    let p = idx(i - a, j, cols);
                    let (v, st) = best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                    cdna_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(v, scoring.frameshift),
                        match st {
                            State::M => CdnaState::CodingM,
                            State::I => CdnaState::CodingI,
                            State::D => CdnaState::CodingD,
                            State::Stop => unreachable!(),
                        },
                        TraceFragment::one(TraceRun {
                            transition_id: 6,
                            op: Op::Frameshift,
                            query_advance: a as u32,
                            target_advance: 0,
                            repeats: 1,
                        }),
                    );
                }
                if j >= a {
                    let p = idx(i, j - a, cols);
                    let (v, st) = best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                    cdna_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(v, scoring.frameshift),
                        match st {
                            State::M => CdnaState::CodingM,
                            State::I => CdnaState::CodingI,
                            State::D => CdnaState::CodingD,
                            State::Stop => unreachable!(),
                        },
                        TraceFragment::one(TraceRun {
                            transition_id: 7,
                            op: Op::Frameshift,
                            query_advance: 0,
                            target_advance: a as u32,
                            repeats: 1,
                        }),
                    );
                }
            }
            if j >= intron.min_len as usize {
                let s = j - intron.min_len as usize;
                let p = idx(i, s, cols);
                let (v, st) = best(&[(u3m[p], State::M), (u3i[p], State::I), (u3d[p], State::D)]);
                if v > 0 {
                    if let Some(d) = splice_score(
                        &target.bases,
                        s,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        w3.insert(IntronCandidate {
                            start: s,
                            score: add(v, intron.open_penalty.saturating_add(d)),
                            state_rank: st.rank(),
                        });
                    }
                }
            }
            if j > intron.max_len as usize {
                w3.expire_before(j - intron.max_len as usize);
            }
            if j >= 2 {
                if let (Some(c), Some(a)) = (
                    w3.best(),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    cdna_relax(
                        &mut u3m[k],
                        &mut pu3m[k],
                        add(c.score, a),
                        match state_from_rank(c.state_rank) {
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
                                target_advance: (j - c.start - 4) as u32,
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
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let p = idx(i - 1, j - 1, cols);
                let (v, st) = best(&[(u3m[p], State::M), (u3i[p], State::I), (u3d[p], State::D)]);
                cdna_relax(
                    &mut u3m[k],
                    &mut pu3m[k],
                    add(
                        v,
                        dna_score(query.bases[i - 1], target.bases[j - 1], scoring),
                    ),
                    match st {
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
            if i > 0 {
                let p = idx(i - 1, j, cols);
                let (candidate, state, transition_id) =
                    if add(u3m[p], scoring.gap_open) >= add(u3i[p], scoring.gap_extend) {
                        (add(u3m[p], scoring.gap_open), CdnaState::U3M, 2)
                    } else {
                        (add(u3i[p], scoring.gap_extend), CdnaState::U3I, 4)
                    };
                cdna_relax(
                    &mut u3i[k],
                    &mut pu3i[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j > 0 {
                let p = idx(i, j - 1, cols);
                let (candidate, state, transition_id) =
                    if add(u3m[p], scoring.gap_open) >= add(u3d[p], scoring.gap_extend) {
                        (add(u3m[p], scoring.gap_open), CdnaState::U3M, 3)
                    } else {
                        (add(u3d[p], scoring.gap_extend), CdnaState::U3D, 5)
                    };
                cdna_relax(
                    &mut u3d[k],
                    &mut pu3d[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            // epsilon transition: coding into the 3-prime UTR match state.
            let (candidate, state) =
                best(&[(cm[k], State::M), (ci[k], State::I), (cd[k], State::D)]);
            cdna_relax(
                &mut u3m[k],
                &mut pu3m[k],
                candidate,
                match state {
                    State::M => CdnaState::CodingM,
                    State::I => CdnaState::CodingI,
                    State::D => CdnaState::CodingD,
                    State::Stop => unreachable!(),
                },
                TraceFragment::empty(),
            );
            for (value, state) in [
                (u3m[k], CdnaState::U3M),
                (u3i[k], CdnaState::U3I),
                (u3d[k], CdnaState::U3D),
                (cm[k], CdnaState::CodingM),
                (ci[k], CdnaState::CodingI),
                (cd[k], CdnaState::CodingD),
            ] {
                if value > end.3 {
                    end = (i, j, state, value);
                }
            }
        }
    }
    CdnaDp {
        cols,
        u5m,
        u5i,
        u5d,
        cm,
        ci,
        cd,
        u3m,
        u3i,
        u3d,
        pu5m,
        pu5i,
        pu5d,
        pcm,
        pci,
        pcd,
        pu3m,
        pu3i,
        pu3d,
        end,
    }
}

/// Score-only cDNA-to-genome composition wrapper.
pub fn cdna2genome_score(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Score {
    let dp = cdna2genome_dp(query, target, scoring, intron, None);
    dp.value(dp.end.2, idx(dp.end.0, dp.end.1, dp.cols))
}

fn cdna_coding_state(state: CdnaState) -> Option<State> {
    match state {
        CdnaState::CodingM => Some(State::M),
        CdnaState::CodingI => Some(State::I),
        CdnaState::CodingD => Some(State::D),
        _ => None,
    }
}

#[allow(clippy::too_many_arguments)]
fn cdna_raw_fragment(
    parent: CdnaParent,
    destination: CdnaState,
    query_end: usize,
    target_end: usize,
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<RawStep> {
    if parent.fragment.atoms.iter().all(Option::is_none) {
        let transition_id = match (parent.state, destination) {
            (CdnaState::U5M, CdnaState::CodingM) => 320,
            (CdnaState::CodingM, CdnaState::U3M) => 321,
            _ => 322,
        };
        return if transition_id == 321 {
            vec![
                RawStep {
                    transition_id,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                },
                RawStep {
                    transition_id: 400,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                },
            ]
        } else {
            vec![RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            }]
        };
    }
    if let (Some(source), Some(destination)) = (
        cdna_coding_state(parent.state),
        cdna_coding_state(destination),
    ) {
        return coding2genome_raw_fragment(
            P2GenomeParent {
                state: source,
                fragment: parent.fragment,
            },
            destination,
            query_end,
            target_end,
            query,
            &target.bases,
            scoring,
            intron,
        );
    }
    let (query_advance, target_advance) = parent.fragment.advances();
    let query_start = query_end - query_advance;
    let target_start = target_end - target_advance;
    let is_u5 = matches!(
        destination,
        CdnaState::U5M | CdnaState::U5I | CdnaState::U5D
    );
    let base = if is_u5 { 300 } else { 400 };
    let mut raw = Vec::new();
    if matches!(destination, CdnaState::U5M | CdnaState::U3M) {
        let epsilon = match parent.state {
            CdnaState::U5I | CdnaState::U3I => Some(base + 4),
            CdnaState::U5D | CdnaState::U3D => Some(base + 7),
            _ => None,
        };
        if let Some(transition_id) = epsilon {
            raw.push(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
        }
    }
    let atoms = parent
        .fragment
        .atoms
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if atoms.iter().any(|run| run.op == Op::Intron) {
        raw.push(RawStep {
            transition_id: base + 10,
            query_advance: 0,
            target_advance: 2,
            score: intron.open_penalty.saturating_add(
                splice_score(
                    &target.bases,
                    target_start,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                )
                .unwrap_or(0),
            ),
        });
        for _ in 0..(target_advance - 4) {
            raw.push(RawStep {
                transition_id: base + 11,
                query_advance: 0,
                target_advance: 1,
                score: 0,
            });
        }
        raw.push(RawStep {
            transition_id: base + 12,
            query_advance: 0,
            target_advance: 2,
            score: splice_score(
                &target.bases,
                target_end - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            )
            .unwrap_or(0),
        });
        return raw;
    }
    let run = atoms[0];
    match run.op {
        Op::Match => raw.push(RawStep {
            transition_id: base + 1,
            query_advance: 1,
            target_advance: 1,
            score: dna_score(
                query.bases[query_start],
                target.bases[target_start],
                scoring,
            ),
        }),
        Op::Insert => raw.push(RawStep {
            transition_id: if run.transition_id == 2 {
                base + 2
            } else {
                base + 3
            },
            query_advance: 1,
            target_advance: 0,
            score: if run.transition_id == 2 {
                scoring.gap_open
            } else {
                scoring.gap_extend
            },
        }),
        Op::Delete => raw.push(RawStep {
            transition_id: if run.transition_id == 3 {
                base + 5
            } else {
                base + 6
            },
            query_advance: 0,
            target_advance: 1,
            score: if run.transition_id == 3 {
                scoring.gap_open
            } else {
                scoring.gap_extend
            },
        }),
        _ => unreachable!("cDNA UTR parent fragment"),
    }
    raw
}

/// Full cDNA-to-genome traceback over the composed UTR/CDS/UTR DP lattice.
fn align_cdna_to_genome_forbidden(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let dp = cdna2genome_dp(query, target, scoring, intron, forbidden);
    let (mut i, mut j, mut state, score) = dp.end;
    let final_state = state;
    let (query_end, target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while let Some(parent) = dp.parent(state, idx(i, j, dp.cols)) {
        let raw_fragment = cdna_raw_fragment(parent, state, i, j, query, target, scoring, intron);
        let (query_advance, target_advance) = parent.fragment.advances();
        debug_assert!(i >= query_advance && j >= target_advance);
        i -= query_advance;
        j -= target_advance;
        fragments.push((parent.fragment, raw_fragment));
        state = parent.state;
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
        score,
        raw_trace,
        trace,
    }
}

pub fn align_cdna_to_genome(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Alignment {
    align_cdna_to_genome_forbidden(query, target, scoring, intron, None)
}

pub fn align_cdna_to_genome_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
) -> Vec<Alignment> {
    enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_cdna_to_genome_forbidden(query, target, scoring, intron, Some(forbidden))
    })
}

pub fn align_cdna_to_genome_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_cdna_to_genome_suboptimal(
                query, target, scoring, intron, threshold,
            ));
        }
    }
    out
}

pub fn align_cdna_to_genome_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets
                .iter()
                .map(move |target| align_cdna_to_genome(query, target, scoring, intron))
        })
        .collect()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GenomeState {
    UtrM,
    UtrI,
    UtrD,
    CdsM,
    CdsI,
    CdsD,
    U3M,
    U3I,
    U3D,
}

/// Score storage used by the genome-to-genome recurrence.  The full solver
/// keeps every row for now; checkpointed execution uses only the predecessor
/// history required by phase introns and frameshifts while retaining the same
/// indexed recurrence below.
enum GenomeScoreMatrix {
    Full {
        cols: usize,
        scores: Vec<Score>,
    },
    Rolling {
        cols: usize,
        rows: Vec<Vec<Score>>,
        row_numbers: Vec<usize>,
    },
}

impl GenomeScoreMatrix {
    fn full(cols: usize, size: usize, initial: Score) -> Self {
        Self::Full {
            cols,
            scores: vec![initial; size],
        }
    }

    fn rolling(cols: usize, history_rows: usize) -> Self {
        Self::Rolling {
            cols,
            rows: (0..history_rows.max(1))
                .map(|_| vec![NEG_INF; cols])
                .collect(),
            row_numbers: vec![usize::MAX; history_rows.max(1)],
        }
    }

    fn start_row(&mut self, row: usize, initial: Score) {
        let Self::Rolling {
            rows, row_numbers, ..
        } = self
        else {
            return;
        };
        let slot = row % rows.len();
        rows[slot].fill(initial);
        row_numbers[slot] = row;
    }

    fn slot(&self, index: usize) -> (usize, usize) {
        match self {
            Self::Full { .. } => unreachable!("full score matrices index directly"),
            Self::Rolling {
                cols,
                rows,
                row_numbers,
            } => {
                let row = index / *cols;
                let column = index % *cols;
                let slot = row % rows.len();
                assert_eq!(
                    row_numbers[slot], row,
                    "genome checkpoint history omitted required row {row}"
                );
                (slot, column)
            }
        }
    }

    fn row_clone(&self, row: usize) -> Vec<Score> {
        match self {
            Self::Full { cols, scores } => scores[row * *cols..(row + 1) * *cols].to_vec(),
            Self::Rolling {
                rows, row_numbers, ..
            } => {
                let slot = row % rows.len();
                assert_eq!(row_numbers[slot], row);
                rows[slot].clone()
            }
        }
    }

    #[allow(dead_code)] // consumed by genome2genome checkpoint block replay
    fn restore_row(&mut self, row: usize, values: &[Score]) {
        match self {
            Self::Full { cols, scores } => {
                assert_eq!(values.len(), *cols);
                scores[row * *cols..(row + 1) * *cols].copy_from_slice(values);
            }
            Self::Rolling {
                cols,
                rows,
                row_numbers,
            } => {
                assert_eq!(values.len(), *cols);
                let slot = row % rows.len();
                rows[slot].copy_from_slice(values);
                row_numbers[slot] = row;
            }
        }
    }
}

impl std::ops::Index<usize> for GenomeScoreMatrix {
    type Output = Score;

    fn index(&self, index: usize) -> &Self::Output {
        match self {
            Self::Full { scores, .. } => &scores[index],
            Self::Rolling { rows, .. } => {
                let (slot, column) = self.slot(index);
                &rows[slot][column]
            }
        }
    }
}

impl std::ops::IndexMut<usize> for GenomeScoreMatrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match self {
            Self::Full { scores, .. } => &mut scores[index],
            Self::Rolling {
                cols,
                rows,
                row_numbers,
            } => {
                let row = index / *cols;
                let column = index % *cols;
                let slot = row % rows.len();
                assert_eq!(
                    row_numbers[slot], row,
                    "genome checkpoint history omitted required row {row}"
                );
                &mut rows[slot][column]
            }
        }
    }
}

/// Parent-pointer storage for genome-to-genome traceback.  Checkpoint block
/// reconstruction writes only its active rows; updates outside that block are
/// routed to a scratch slot because their scores remain available in the
/// checkpoint history.
enum GenomeParentMatrix {
    Full(Vec<Option<GenomeParent>>),
    Block {
        cols: usize,
        first_row: usize,
        last_row: usize,
        parents: Vec<Option<GenomeParent>>,
        empty: Option<GenomeParent>,
        discarded: Option<GenomeParent>,
    },
}

impl GenomeParentMatrix {
    fn full(size: usize) -> Self {
        Self::Full(vec![None; size])
    }

    #[allow(dead_code)] // consumed by genome2genome block reconstruction
    fn block(cols: usize, first_row: usize, last_row: usize) -> Self {
        Self::Block {
            cols,
            first_row,
            last_row,
            parents: vec![None; (last_row - first_row + 1) * cols],
            empty: None,
            discarded: None,
        }
    }

    fn block_index(&self, index: usize) -> Option<usize> {
        let Self::Block {
            cols,
            first_row,
            last_row,
            ..
        } = self
        else {
            return None;
        };
        let row = index / *cols;
        if row >= *first_row && row <= *last_row {
            Some((row - *first_row) * *cols + index % *cols)
        } else {
            None
        }
    }
}

impl std::ops::Index<usize> for GenomeParentMatrix {
    type Output = Option<GenomeParent>;

    fn index(&self, index: usize) -> &Self::Output {
        match self {
            Self::Full(parents) => &parents[index],
            Self::Block { parents, empty, .. } => self
                .block_index(index)
                .map(|offset| &parents[offset])
                .unwrap_or(empty),
        }
    }
}

impl std::ops::IndexMut<usize> for GenomeParentMatrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let offset = self.block_index(index);
        match self {
            Self::Full(parents) => &mut parents[index],
            Self::Block {
                parents, discarded, ..
            } => offset
                .map(|offset| &mut parents[offset])
                .unwrap_or(discarded),
        }
    }
}

#[allow(clippy::too_many_arguments, dead_code)] // consumed by block replay
fn restore_genome_checkpoint_scores(
    checkpoint: &GenomeCheckpoint,
    u5m: &mut GenomeScoreMatrix,
    u5i: &mut GenomeScoreMatrix,
    u5d: &mut GenomeScoreMatrix,
    cm: &mut GenomeScoreMatrix,
    ci: &mut GenomeScoreMatrix,
    cd: &mut GenomeScoreMatrix,
    u3m: &mut GenomeScoreMatrix,
    u3i: &mut GenomeScoreMatrix,
    u3d: &mut GenomeScoreMatrix,
) {
    for row in &checkpoint.history {
        u5m.restore_row(row.row, &row.u5m);
        u5i.restore_row(row.row, &row.u5i);
        u5d.restore_row(row.row, &row.u5d);
        cm.restore_row(row.row, &row.cm);
        ci.restore_row(row.row, &row.ci);
        cd.restore_row(row.row, &row.cd);
        u3m.restore_row(row.row, &row.u3m);
        u3i.restore_row(row.row, &row.u3i);
        u3d.restore_row(row.row, &row.u3d);
    }
}

#[allow(dead_code)] // consumed by genome2genome checkpoint block replay
fn restore_genome_checkpoint_queues(checkpoint: &GenomeCheckpoint) -> GenomeIntronQueues {
    checkpoint.queues.clone()
}
impl From<State> for GenomeState {
    fn from(value: State) -> Self {
        match value {
            State::M => Self::UtrM,
            State::I => Self::UtrI,
            State::D => Self::UtrD,
            State::Stop => unreachable!("traceback parent cannot be stop"),
        }
    }
}

#[derive(Clone, Copy)]
struct GenomeParent {
    state: GenomeState,
    fragment: CompactGenomeFragment,
}

/// A genome-to-genome parent never needs the `repeats` field: every atomic
/// transition is materialized individually by the recurrence.  Keeping the
/// five atoms in this packed form halves the resident parent-matrix payload
/// while preserving enough information to recreate the public trace exactly.
#[derive(Clone, Copy)]
struct CompactGenomeRun {
    transition_id: u16,
    op: Op,
    query_advance: u32,
    target_advance: u32,
}

#[derive(Clone, Copy)]
struct CompactGenomeFragment {
    atoms: [Option<CompactGenomeRun>; 5],
}

impl CompactGenomeFragment {
    fn from_trace_fragment(fragment: TraceFragment) -> Self {
        Self {
            atoms: fragment.atoms.map(|atom| {
                atom.map(|run| {
                    debug_assert_eq!(run.repeats, 1);
                    CompactGenomeRun {
                        transition_id: run.transition_id,
                        op: run.op,
                        query_advance: run.query_advance,
                        target_advance: run.target_advance,
                    }
                })
            }),
        }
    }

    fn expand(self) -> TraceFragment {
        TraceFragment {
            atoms: self.atoms.map(|atom| {
                atom.map(|run| TraceRun {
                    transition_id: run.transition_id,
                    op: run.op,
                    query_advance: run.query_advance,
                    target_advance: run.target_advance,
                    repeats: 1,
                })
            }),
        }
    }

    fn advances(self) -> (usize, usize) {
        self.atoms
            .into_iter()
            .flatten()
            .fold((0, 0), |(q, t), run| {
                (
                    q + run.query_advance as usize,
                    t + run.target_advance as usize,
                )
            })
    }
}

#[derive(Clone, Copy)]
struct JointIntronCandidate {
    query_start: usize,
    target_start: usize,
    score: Score,
    state_rank: u8,
}

#[derive(Clone, Default)]
struct JointIntronWindow {
    entries: VecDeque<JointIntronCandidate>,
}
impl JointIntronWindow {
    fn insert(&mut self, candidate: JointIntronCandidate) {
        while self.entries.back().is_some_and(|entry| {
            entry.score < candidate.score
                || (entry.score == candidate.score && entry.state_rank <= candidate.state_rank)
        }) {
            self.entries.pop_back();
        }
        self.entries.push_back(candidate);
    }
    fn expire_before(&mut self, start: usize) {
        while self
            .entries
            .front()
            .is_some_and(|entry| entry.query_start < start)
        {
            self.entries.pop_front();
        }
    }
    fn best(&self) -> Option<JointIntronCandidate> {
        self.entries.front().copied()
    }

    #[allow(dead_code)] // consumed when genome2genome checkpoint snapshots are wired in
    fn checkpoint_bytes(&self) -> usize {
        size_of::<Self>() + self.entries.len() * size_of::<JointIntronCandidate>()
    }
}

/// All cross-query-row intron candidates needed by `genome2genome`.
///
/// Keeping these together makes a checkpoint self-contained: score rows alone
/// are insufficient because query and joint introns retain monotonic candidate
/// queues for every target column and coding phase.
#[derive(Clone)]
struct GenomeIntronQueues {
    query: Vec<IntronCandidateWindow>,
    joint: Vec<JointIntronWindow>,
    u3_query: Vec<IntronCandidateWindow>,
    u3_joint: Vec<JointIntronWindow>,
    query_coding: [Vec<IntronCandidateWindow>; 3],
    joint_coding: [Vec<JointIntronWindow>; 3],
}

impl GenomeIntronQueues {
    fn new(columns: usize) -> Self {
        Self {
            query: (0..columns)
                .map(|_| IntronCandidateWindow::default())
                .collect(),
            joint: (0..columns).map(|_| JointIntronWindow::default()).collect(),
            u3_query: (0..columns)
                .map(|_| IntronCandidateWindow::default())
                .collect(),
            u3_joint: (0..columns).map(|_| JointIntronWindow::default()).collect(),
            query_coding: std::array::from_fn(|_| {
                (0..columns)
                    .map(|_| IntronCandidateWindow::default())
                    .collect()
            }),
            joint_coding: std::array::from_fn(|_| {
                (0..columns).map(|_| JointIntronWindow::default()).collect()
            }),
        }
    }

    /// Exact payload size of the queue state that a genome2genome checkpoint
    /// must preserve.  Candidate deques span arbitrarily many query rows, so
    /// treating this as a fixed number of score rows would under-estimate
    /// checkpoint memory and make an exact reconstruction impossible.
    #[allow(dead_code)] // consumed when genome2genome checkpoint snapshots are wired in
    fn checkpoint_bytes(&self) -> usize {
        self.query
            .iter()
            .map(IntronCandidateWindow::checkpoint_bytes)
            .chain(
                self.query_coding
                    .iter()
                    .flatten()
                    .map(IntronCandidateWindow::checkpoint_bytes),
            )
            .chain(
                self.u3_query
                    .iter()
                    .map(IntronCandidateWindow::checkpoint_bytes),
            )
            .chain(self.joint.iter().map(JointIntronWindow::checkpoint_bytes))
            .chain(
                self.u3_joint
                    .iter()
                    .map(JointIntronWindow::checkpoint_bytes),
            )
            .chain(
                self.joint_coding
                    .iter()
                    .flatten()
                    .map(JointIntronWindow::checkpoint_bytes),
            )
            .sum()
    }
}

/// Number of full score rows retained by a genome2genome checkpoint.
///
/// Phase-aware query and joint intron candidates read as far back as
/// `min_intron + 3` rows; keeping the source row as well makes that
/// `min_intron + 4`.  The four/five-base frameshift edges independently need
/// six rows.  This must remain dynamic rather than adopting the five-row
/// coding2genome history.
fn genome_checkpoint_history_rows(intron: IntronScoring) -> usize {
    (intron.min_len as usize + 4).max(6)
}

/// Conservative payload bound for the cross-row intron queues at one
/// `genome2genome` checkpoint.
///
/// Every family retains candidates whose donor coordinate lies in the
/// inclusive `min_len..=max_len` intron window.  There are five
/// single-coordinate families (two UTRs plus three coding phases) and five
/// joint-coordinate families, each for every target column.  The bound is
/// intentionally expressed in live candidate payload rather than deque
/// capacity: cloned checkpoints allocate only semantic entries.
#[allow(dead_code)] // consumed when genome2genome checkpoint planning is wired in
fn genome_checkpoint_queue_upper_bound_bytes(columns: usize, intron: IntronScoring) -> u128 {
    let span = (intron.max_len as usize)
        .saturating_sub(intron.min_len as usize)
        .saturating_add(1) as u128;
    let columns = columns as u128;
    let single =
        size_of::<IntronCandidateWindow>() as u128 + span * size_of::<IntronCandidate>() as u128;
    let joint =
        size_of::<JointIntronWindow>() as u128 + span * size_of::<JointIntronCandidate>() as u128;
    columns * 5 * (single + joint)
}

fn genome_relax<S: Into<GenomeState>>(
    score: &mut Score,
    parent: &mut Option<GenomeParent>,
    candidate: Score,
    state: S,
    fragment: TraceFragment,
) {
    // Equal-scoring consuming edges replace a prior edge, matching C4's
    // transition-order tie resolution while preserving an existing path over
    // a later epsilon edge.
    if candidate > *score || (candidate == *score && fragment.atoms.iter().any(Option::is_some)) {
        *score = candidate;
        *parent = Some(GenomeParent {
            state: state.into(),
            fragment: CompactGenomeFragment::from_trace_fragment(fragment),
        });
    }
}

fn genome_state_is_coding(state: GenomeState) -> bool {
    matches!(
        state,
        GenomeState::CdsM | GenomeState::CdsI | GenomeState::CdsD
    )
}

#[allow(clippy::too_many_arguments)]
fn genome_raw_fragment(
    parent: GenomeParent,
    destination: GenomeState,
    query_end: usize,
    target_end: usize,
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<RawStep> {
    let (query_advance, target_advance) = parent.fragment.advances();
    let query_start = query_end - query_advance;
    let target_start = target_end - target_advance;
    let atoms = parent
        .fragment
        .expand()
        .atoms
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if atoms.is_empty() {
        return vec![RawStep {
            transition_id: if genome_state_is_coding(destination) {
                50
            } else {
                58
            },
            query_advance: 0,
            target_advance: 0,
            score: 0,
        }];
    }

    let first_id = atoms
        .iter()
        .map(|run| run.transition_id)
        .find(|id| matches!(id, 40 | 43 | 46 | 60 | 65 | 70))
        .unwrap_or(atoms[0].transition_id);
    if atoms.iter().any(|run| run.op == Op::Intron) {
        let phase = atoms
            .iter()
            .find(|run| run.op == Op::SplitCodon)
            .map_or(0, |run| run.query_advance as usize);
        let coding = first_id >= 60;
        let post = if coding { 3 - phase } else { 0 };
        let on_query = matches!(first_id, 43 | 46 | 65 | 70);
        let on_target = matches!(first_id, 40 | 46 | 60 | 70);
        let query_donor = query_start + phase;
        let target_donor = target_start + phase;
        let query_post = query_end - post;
        let target_post = target_end - post;
        let (pre_id, query_loop_id, target_loop_id, post_id, pre_split_id, post_split_id) =
            match first_id {
                40 => (40, 0, 41, 42, 0, 0),
                43 => (43, 44, 0, 45, 0, 0),
                46 => (46, 47, 49, 48, 0, 0),
                60 => (
                    110 + phase as u16,
                    0,
                    120 + phase as u16,
                    130 + phase as u16,
                    100 + phase as u16,
                    140 + phase as u16,
                ),
                65 => (
                    160 + phase as u16,
                    170 + phase as u16,
                    0,
                    180 + phase as u16,
                    150 + phase as u16,
                    190 + phase as u16,
                ),
                70 => (
                    210 + phase as u16,
                    220 + phase as u16,
                    230 + phase as u16,
                    240 + phase as u16,
                    200 + phase as u16,
                    250 + phase as u16,
                ),
                _ => unreachable!("genome intron fragment"),
            };
        let mut raw = Vec::new();
        if coding && phase > 0 {
            raw.push(RawStep {
                transition_id: pre_split_id,
                query_advance: phase as u32,
                target_advance: phase as u32,
                score: 0,
            });
        }
        let donor_score = intron
            .open_penalty
            .saturating_add(if on_query {
                splice_score(
                    &query.bases,
                    query_donor,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                )
                .unwrap_or(0)
            } else {
                0
            })
            .saturating_add(if on_target {
                splice_score(
                    &target.bases,
                    target_donor,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                )
                .unwrap_or(0)
            } else {
                0
            });
        raw.push(RawStep {
            transition_id: pre_id,
            query_advance: if on_query { 2 } else { 0 },
            target_advance: if on_target { 2 } else { 0 },
            score: donor_score,
        });
        if on_target {
            for _ in 0..target_post.saturating_sub(target_donor + 4) {
                raw.push(RawStep {
                    transition_id: target_loop_id,
                    query_advance: 0,
                    target_advance: 1,
                    score: 0,
                });
            }
        }
        if on_query {
            for _ in 0..query_post.saturating_sub(query_donor + 4) {
                raw.push(RawStep {
                    transition_id: query_loop_id,
                    query_advance: 1,
                    target_advance: 0,
                    score: 0,
                });
            }
        }

        let acceptor_score = (if on_query {
            splice_score(
                &query.bases,
                query_post - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            )
            .unwrap_or(0)
        } else {
            0
        })
        .saturating_add(if on_target {
            splice_score(
                &target.bases,
                target_post - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            )
            .unwrap_or(0)
        } else {
            0
        });
        raw.push(RawStep {
            transition_id: post_id,
            query_advance: if on_query { 2 } else { 0 },
            target_advance: if on_target { 2 } else { 0 },
            score: acceptor_score,
        });
        if coding {
            let mut query_codon = Vec::with_capacity(3);
            query_codon.extend_from_slice(&query.bases[query_start..query_donor]);
            query_codon.extend_from_slice(&query.bases[query_post..query_end]);
            let mut target_codon = Vec::with_capacity(3);
            target_codon.extend_from_slice(&target.bases[target_start..target_donor]);
            target_codon.extend_from_slice(&target.bases[target_post..target_end]);
            raw.push(RawStep {
                transition_id: post_split_id,
                query_advance: post as u32,
                target_advance: post as u32,
                score: protein_score(
                    translate_dna(&query_codon, 0)[0],
                    translate_dna(&target_codon, 0)[0],
                    scoring,
                ),
            });
        }
        return raw;
    }

    let run = atoms[0];
    let mut raw = Vec::new();
    if matches!(destination, GenomeState::UtrM | GenomeState::CdsM) {
        let close_id = match parent.state {
            GenomeState::UtrI => Some(35),
            GenomeState::UtrD => Some(36),
            GenomeState::CdsI => Some(59),
            GenomeState::CdsD => Some(79),
            _ => None,
        };
        if let Some(transition_id) = close_id {
            raw.push(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
        }
    }
    match run.op {
        Op::Match if run.query_advance == 1 => raw.push(RawStep {
            transition_id: 30,
            query_advance: 1,
            target_advance: 1,
            score: dna_score(
                query.bases[query_start],
                target.bases[target_start],
                scoring,
            ),
        }),
        Op::Match => raw.push(RawStep {
            transition_id: 51,
            query_advance: 3,
            target_advance: 3,
            score: protein_score(
                translate_dna(&query.bases[query_start..query_end], 0)[0],
                translate_dna(&target.bases[target_start..target_end], 0)[0],
                scoring,
            ),
        }),
        Op::Insert | Op::Delete => raw.push(RawStep {
            transition_id: run.transition_id,
            query_advance: run.query_advance,
            target_advance: run.target_advance,
            score: match run.transition_id {
                31 | 33 => scoring.gap_open,
                32 | 34 => scoring.gap_extend,
                52 | 53 => scoring.codon_gap_open,
                54 | 55 => scoring.codon_gap_extend,
                _ => unreachable!("genome gap transition"),
            },
        }),
        Op::Frameshift => {
            let advance = run.query_advance.max(run.target_advance);
            let short = advance % 3;
            let on_query = run.query_advance > 0;
            raw.push(RawStep {
                transition_id: match (on_query, short) {
                    (true, 1) => 80,
                    (true, _) => 81,
                    (false, 1) => 84,
                    (false, _) => 85,
                },
                query_advance: if on_query { short } else { 0 },
                target_advance: if on_query { 0 } else { short },
                score: scoring.frameshift,
            });
            raw.push(RawStep {
                transition_id: match (on_query, advance > 3) {
                    (true, false) => 82,
                    (true, true) => 83,
                    (false, false) => 86,
                    (false, true) => 87,
                },
                query_advance: if on_query && advance > 3 { 3 } else { 0 },
                target_advance: if !on_query && advance > 3 { 3 } else { 0 },
                score: 0,
            });
        }
        _ => unreachable!("genome parent fragment"),
    }
    raw
}

#[allow(clippy::too_many_arguments)]
fn genome_traceback_fragments(
    mut i: usize,
    mut j: usize,
    mut state: GenomeState,
    cols: usize,
    pmat: &GenomeParentMatrix,
    pins: &GenomeParentMatrix,
    pdel: &GenomeParentMatrix,
    pcm: &GenomeParentMatrix,
    pci: &GenomeParentMatrix,
    pcd: &GenomeParentMatrix,
    pu3m: &GenomeParentMatrix,
    pu3i: &GenomeParentMatrix,
    pu3d: &GenomeParentMatrix,
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> (
    usize,
    usize,
    GenomeState,
    Vec<(TraceFragment, Vec<RawStep>)>,
) {
    let mut fragments = Vec::new();
    loop {
        let parent = match state {
            GenomeState::UtrM => pmat[idx(i, j, cols)],
            GenomeState::UtrI => pins[idx(i, j, cols)],
            GenomeState::UtrD => pdel[idx(i, j, cols)],
            GenomeState::CdsM => pcm[idx(i, j, cols)],
            GenomeState::CdsI => pci[idx(i, j, cols)],
            GenomeState::CdsD => pcd[idx(i, j, cols)],
            GenomeState::U3M => pu3m[idx(i, j, cols)],
            GenomeState::U3I => pu3i[idx(i, j, cols)],
            GenomeState::U3D => pu3d[idx(i, j, cols)],
        };
        let Some(parent) = parent else { break };
        let raw_fragment = genome_raw_fragment(parent, state, i, j, query, target, scoring, intron);
        let (qa, ta) = parent.fragment.advances();
        debug_assert!(i >= qa && j >= ta);
        i -= qa;
        j -= ta;
        fragments.push((parent.fragment.expand(), raw_fragment));
        state = parent.state;
    }
    fragments.reverse();
    (i, j, state, fragments)
}

#[allow(clippy::needless_range_loop, clippy::too_many_arguments)]
/// Local genome-to-genome Viterbi with affine DNA gaps plus query-only,
/// target-only, and synchronised (joint) introns.  A joint intron has one
/// opening penalty and splice scores on both sequences, exactly as upstream
/// `Intron_create("joint", TRUE, TRUE, TRUE)` specifies.
fn align_genome_to_genome_forbidden(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    forbidden: Option<&HashSet<(usize, usize)>>,
    rolling_scores: bool,
    parent_rows: Option<(usize, usize)>,
    forced_end: Option<(usize, usize, GenomeState)>,
    traceback_start_sink: Option<&mut (usize, usize, GenomeState)>,
    mut checkpoint_sink: Option<&mut Vec<GenomeCheckpoint>>,
    checkpoint_stride: usize,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let (n, m, cols) = (
        query.bases.len(),
        target.bases.len(),
        target.bases.len() + 1,
    );
    let size = (n + 1) * (m + 1);
    let mut mat = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, 0)
    };
    let mut ins = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut del = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let parent_matrix = || match parent_rows {
        Some((first_row, last_row)) => GenomeParentMatrix::block(cols, first_row, last_row),
        None => GenomeParentMatrix::full(size),
    };
    let mut pmat = parent_matrix();
    let mut pins = parent_matrix();
    let mut pdel = parent_matrix();
    let mut cm = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut ci = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut cd = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut pcm = parent_matrix();
    let mut pci = parent_matrix();
    let mut pcd = parent_matrix();
    // `genome2genome` composes the pre- and post-CDS UTR graphs separately.
    // Keeping the post-CDS states distinct prevents an epsilon CDS exit from
    // re-entering the coding graph, which would create paths absent upstream.
    let mut u3m = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut u3i = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut u3d = if rolling_scores {
        GenomeScoreMatrix::rolling(cols, genome_checkpoint_history_rows(intron))
    } else {
        GenomeScoreMatrix::full(cols, size, NEG_INF)
    };
    let mut pu3m = parent_matrix();
    let mut pu3i = parent_matrix();
    let mut pu3d = parent_matrix();
    let mut queues = GenomeIntronQueues::new(m + 1);
    let mut end = (0usize, 0usize, GenomeState::UtrM, 0);
    let min_len = intron.min_len as usize;
    let max_len = intron.max_len as usize;
    debug_assert!(genome_checkpoint_history_rows(intron) > min_len + 3);

    for i in 0..=n {
        mat.start_row(i, 0);
        ins.start_row(i, NEG_INF);
        del.start_row(i, NEG_INF);
        cm.start_row(i, NEG_INF);
        ci.start_row(i, NEG_INF);
        cd.start_row(i, NEG_INF);
        u3m.start_row(i, NEG_INF);
        u3i.start_row(i, NEG_INF);
        u3d.start_row(i, NEG_INF);
        // Make query-side and joint-intron starts eligible once their minimum
        // query span has been reached.  Their source row is fully resolved.
        if i >= min_len {
            let source_i = i - min_len;
            for j in 0..=m {
                let source = idx(source_i, j, cols);
                let (value, state) = best(&[
                    (mat[source], State::M),
                    (ins[source], State::I),
                    (del[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor) = splice_score(
                        &query.bases,
                        source_i,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        queues.query[j].insert(IntronCandidate {
                            start: source_i,
                            score: add(value, intron.open_penalty.saturating_add(donor)),
                            state_rank: state.rank(),
                        });
                        if let Some(target_donor) = splice_score(
                            &target.bases,
                            j,
                            SpliceType::DonorForward,
                            intron.force_gtag,
                        ) {
                            queues.joint[j].insert(JointIntronCandidate {
                                query_start: source_i,
                                target_start: j,
                                score: add(
                                    value,
                                    intron
                                        .open_penalty
                                        .saturating_add(donor)
                                        .saturating_add(target_donor),
                                ),
                                state_rank: state.rank(),
                            });
                        }
                    }
                }
                if i > max_len {
                    queues.query[j].expire_before(i - max_len);
                    queues.joint[j].expire_before(i - max_len);
                }
                let (value, state) = best(&[
                    (u3m[source], State::M),
                    (u3i[source], State::I),
                    (u3d[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor) = splice_score(
                        &query.bases,
                        source_i,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        queues.u3_query[j].insert(IntronCandidate {
                            start: source_i,
                            score: add(value, intron.open_penalty.saturating_add(donor)),
                            state_rank: state.rank(),
                        });
                        if let Some(target_donor) = splice_score(
                            &target.bases,
                            j,
                            SpliceType::DonorForward,
                            intron.force_gtag,
                        ) {
                            queues.u3_joint[j].insert(JointIntronCandidate {
                                query_start: source_i,
                                target_start: j,
                                score: add(
                                    value,
                                    intron
                                        .open_penalty
                                        .saturating_add(donor)
                                        .saturating_add(target_donor),
                                ),
                                state_rank: state.rank(),
                            });
                        }
                    }
                }
                if i > max_len {
                    queues.u3_query[j].expire_before(i - max_len);
                    queues.u3_joint[j].expire_before(i - max_len);
                }
            }
        }

        // Make coding joint-intron starts eligible for each phase. The source
        // row precedes the current row, so all source CDS states are final.
        for phase in 0..3 {
            let post = 3 - phase;
            if i >= min_len + 3 {
                let post_start = i - post;
                let donor = post_start - min_len;
                let source_i = donor - phase;
                if let Some(query_donor) = splice_score(
                    &query.bases,
                    donor,
                    SpliceType::DonorForward,
                    intron.force_gtag,
                ) {
                    for target_donor in phase..=m {
                        let source = idx(source_i, target_donor - phase, cols);
                        let (value, state) = best(&[
                            (cm[source], State::M),
                            (ci[source], State::I),
                            (cd[source], State::D),
                        ]);
                        if value > 0 {
                            if let Some(target_score) = splice_score(
                                &target.bases,
                                target_donor,
                                SpliceType::DonorForward,
                                intron.force_gtag,
                            ) {
                                queues.joint_coding[phase][target_donor].insert(
                                    JointIntronCandidate {
                                        query_start: donor,
                                        target_start: target_donor,
                                        score: add(
                                            value,
                                            intron
                                                .open_penalty
                                                .saturating_add(query_donor)
                                                .saturating_add(target_score),
                                        ),
                                        state_rank: state.rank(),
                                    },
                                );
                            }
                        }
                    }
                }
            }
            if i >= min_len + 3 {
                let post_start = i - (3 - phase);
                if post_start > max_len {
                    for column in &mut queues.joint_coding[phase] {
                        column.expire_before(post_start - max_len);
                    }
                }
            }
        }

        // At this query coordinate, a monotonic target window turns the
        // per-column best joint starts into the required 2-D rectangle max.
        let mut joint_targets = JointIntronWindow::default();
        let mut target_window = IntronCandidateWindow::default();
        let mut u3_joint_targets = JointIntronWindow::default();
        let mut u3_target_window = IntronCandidateWindow::default();
        let mut coding_phase0 = IntronCandidateWindow::default();
        let mut coding_phase1 = IntronCandidateWindow::default();
        let mut coding_phase2 = IntronCandidateWindow::default();
        let mut joint_coding_targets: [JointIntronWindow; 3] =
            std::array::from_fn(|_| JointIntronWindow::default());
        for j in 0..=m {
            let k = idx(i, j, cols);
            if j >= min_len {
                let source_j = j - min_len;
                let source = idx(i, source_j, cols);
                let (value, state) = best(&[
                    (mat[source], State::M),
                    (ins[source], State::I),
                    (del[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor) = splice_score(
                        &target.bases,
                        source_j,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        target_window.insert(IntronCandidate {
                            start: source_j,
                            score: add(value, intron.open_penalty.saturating_add(donor)),
                            state_rank: state.rank(),
                        });
                    }
                }
            }
            if j > max_len {
                target_window.expire_before(j - max_len);
            }
            if j >= min_len {
                let source_j = j - min_len;
                if let Some(candidate) = queues.joint[source_j].best() {
                    joint_targets.insert(candidate);
                }
            }
            if j > max_len {
                while joint_targets
                    .entries
                    .front()
                    .is_some_and(|entry| entry.target_start < j - max_len)
                {
                    joint_targets.entries.pop_front();
                }
            }
            // The post-CDS UTR has the same UTR intron graph as the pre-CDS
            // UTR, but it is a distinct one-way subgraph.
            if j >= min_len {
                let source_j = j - min_len;
                let source = idx(i, source_j, cols);
                let (value, state) = best(&[
                    (u3m[source], State::M),
                    (u3i[source], State::I),
                    (u3d[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor) = splice_score(
                        &target.bases,
                        source_j,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        u3_target_window.insert(IntronCandidate {
                            start: source_j,
                            score: add(value, intron.open_penalty.saturating_add(donor)),
                            state_rank: state.rank(),
                        });
                    }
                }
                if let Some(candidate) = queues.u3_joint[source_j].best() {
                    u3_joint_targets.insert(candidate);
                }
            }
            if j > max_len {
                u3_target_window.expire_before(j - max_len);
                while u3_joint_targets
                    .entries
                    .front()
                    .is_some_and(|entry| entry.target_start < j - max_len)
                {
                    u3_joint_targets.entries.pop_front();
                }
            }

            let u3_state = |rank| match state_from_rank(rank) {
                State::M => GenomeState::U3M,
                State::I => GenomeState::U3I,
                State::D => GenomeState::U3D,
                State::Stop => unreachable!(),
            };
            if j >= 2 {
                if let (Some(candidate), Some(acceptor)) = (
                    u3_target_window.best(),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut u3m[k],
                        &mut pu3m[k],
                        add(candidate.score, acceptor),
                        u3_state(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 40,
                                op: Op::Splice5,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 41,
                                op: Op::Intron,
                                query_advance: 0,
                                target_advance: (j - candidate.start - 4) as u32,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 42,
                                op: Op::Splice3,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            if i >= 2 {
                if let (Some(candidate), Some(acceptor)) = (
                    queues.u3_query[j].best(),
                    splice_score(
                        &query.bases,
                        i - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut u3m[k],
                        &mut pu3m[k],
                        add(candidate.score, acceptor),
                        u3_state(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 43,
                                op: Op::Splice5,
                                query_advance: 2,
                                target_advance: 0,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 44,
                                op: Op::Intron,
                                query_advance: (i - candidate.start - 4) as u32,
                                target_advance: 0,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 45,
                                op: Op::Splice3,
                                query_advance: 2,
                                target_advance: 0,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            if i >= 2 && j >= 2 {
                if let (Some(candidate), Some(query_acceptor), Some(target_acceptor)) = (
                    u3_joint_targets.best(),
                    splice_score(
                        &query.bases,
                        i - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut u3m[k],
                        &mut pu3m[k],
                        add(add(candidate.score, query_acceptor), target_acceptor),
                        u3_state(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 46,
                                op: Op::Splice5,
                                query_advance: 2,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 47,
                                op: Op::Intron,
                                query_advance: (i - candidate.query_start - 4) as u32,
                                target_advance: (j - candidate.target_start - 4) as u32,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 48,
                                op: Op::Splice3,
                                query_advance: 2,
                                target_advance: 2,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }

            if j >= 2 {
                if let (Some(candidate), Some(acceptor)) = (
                    target_window.best(),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut mat[k],
                        &mut pmat[k],
                        add(candidate.score, acceptor),
                        state_from_rank(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 40,
                                op: Op::Splice5,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 41,
                                op: Op::Intron,
                                query_advance: 0,
                                target_advance: (j - candidate.start - 4) as u32,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 42,
                                op: Op::Splice3,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            if i >= 2 {
                if let (Some(candidate), Some(acceptor)) = (
                    queues.query[j].best(),
                    splice_score(
                        &query.bases,
                        i - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut mat[k],
                        &mut pmat[k],
                        add(candidate.score, acceptor),
                        state_from_rank(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 43,
                                op: Op::Splice5,
                                query_advance: 2,
                                target_advance: 0,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 44,
                                op: Op::Intron,
                                query_advance: (i - candidate.start - 4) as u32,
                                target_advance: 0,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 45,
                                op: Op::Splice3,
                                query_advance: 2,
                                target_advance: 0,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            if i >= 2 && j >= 2 {
                if let (Some(candidate), Some(query_acceptor), Some(target_acceptor)) = (
                    joint_targets.best(),
                    splice_score(
                        &query.bases,
                        i - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                    splice_score(
                        &target.bases,
                        j - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    genome_relax(
                        &mut mat[k],
                        &mut pmat[k],
                        add(add(candidate.score, query_acceptor), target_acceptor),
                        state_from_rank(candidate.state_rank),
                        TraceFragment::intron(
                            TraceRun {
                                transition_id: 46,
                                op: Op::Splice5,
                                query_advance: 2,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 47,
                                op: Op::Intron,
                                query_advance: (i - candidate.query_start - 4) as u32,
                                target_advance: (j - candidate.target_start - 4) as u32,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 48,
                                op: Op::Splice3,
                                query_advance: 2,
                                target_advance: 2,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let source = idx(i - 1, j - 1, cols);
                let (value, state) = best(&[
                    (mat[source], State::M),
                    (ins[source], State::I),
                    (del[source], State::D),
                ]);
                genome_relax(
                    &mut mat[k],
                    &mut pmat[k],
                    add(
                        value,
                        dna_score(query.bases[i - 1], target.bases[j - 1], scoring),
                    ),
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id: 30,
                        op: Op::Match,
                        query_advance: 1,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            if i > 0 {
                let source = idx(i - 1, j, cols);
                let (candidate, state, transition_id) =
                    if add(mat[source], scoring.gap_open) >= add(ins[source], scoring.gap_extend) {
                        (add(mat[source], scoring.gap_open), State::M, 31)
                    } else {
                        (add(ins[source], scoring.gap_extend), State::I, 32)
                    };
                genome_relax(
                    &mut ins[k],
                    &mut pins[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j > 0 {
                let source = idx(i, j - 1, cols);
                let (candidate, state, transition_id) =
                    if add(mat[source], scoring.gap_open) >= add(del[source], scoring.gap_extend) {
                        (add(mat[source], scoring.gap_open), State::M, 33)
                    } else {
                        (add(del[source], scoring.gap_extend), State::D, 34)
                    };
                genome_relax(
                    &mut del[k],
                    &mut pdel[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            // Enter the coding subgraph from a UTR match/gap state.
            let (candidate, state) =
                best(&[(mat[k], State::M), (ins[k], State::I), (del[k], State::D)]);
            genome_relax(
                &mut cm[k],
                &mut pcm[k],
                candidate,
                state,
                TraceFragment::empty(),
            );
            // Query-only coding introns in phases 0, 1, and 2.
            for phase in 0..3 {
                let post = 3 - phase;
                if i >= min_len + 3 && j >= 3 {
                    let post_start = i - post;
                    let donor = post_start - min_len;
                    let source = idx(donor - phase, j - 3, cols);
                    let (value, source_state) = best(&[
                        (cm[source], State::M),
                        (ci[source], State::I),
                        (cd[source], State::D),
                    ]);
                    if value > 0 {
                        if let Some(donor_score) = splice_score(
                            &query.bases,
                            donor,
                            SpliceType::DonorForward,
                            intron.force_gtag,
                        ) {
                            queues.query_coding[phase][j].insert(IntronCandidate {
                                start: donor,
                                score: add(value, intron.open_penalty.saturating_add(donor_score)),
                                state_rank: source_state.rank(),
                            });
                        }
                    }
                    if post_start > max_len {
                        queues.query_coding[phase][j].expire_before(post_start - max_len);
                    }
                    if let (Some(candidate), Some(acceptor)) = (
                        queues.query_coding[phase][j].best(),
                        splice_score(
                            &query.bases,
                            post_start - 2,
                            SpliceType::AcceptorForward,
                            intron.force_gtag,
                        ),
                    ) {
                        if forbidden.is_some_and(|pairs| {
                            pairs.contains(&(candidate.start - phase, j - 3))
                                || pairs.contains(&(post_start, j - post))
                        }) {
                            continue;
                        }
                        let mut query_codon = Vec::with_capacity(3);
                        query_codon.extend_from_slice(
                            &query.bases[candidate.start - phase..candidate.start],
                        );
                        query_codon.extend_from_slice(&query.bases[post_start..i]);
                        let state = match state_from_rank(candidate.state_rank) {
                            State::M => GenomeState::CdsM,
                            State::I => GenomeState::CdsI,
                            State::D => GenomeState::CdsD,
                            State::Stop => unreachable!(),
                        };
                        let splice5 = TraceRun {
                            transition_id: 65,
                            op: Op::Splice5,
                            query_advance: 2,
                            target_advance: 0,
                            repeats: 1,
                        };
                        let loop_run = TraceRun {
                            transition_id: 66,
                            op: Op::Intron,
                            query_advance: (post_start - candidate.start - 4) as u32,
                            target_advance: 0,
                            repeats: 1,
                        };
                        let splice3 = TraceRun {
                            transition_id: 67,
                            op: Op::Splice3,
                            query_advance: 2,
                            target_advance: 0,
                            repeats: 1,
                        };
                        let fragment = if phase == 0 {
                            TraceFragment::intron_match(
                                splice5,
                                loop_run,
                                splice3,
                                TraceRun {
                                    transition_id: 51,
                                    op: Op::Match,
                                    query_advance: 3,
                                    target_advance: 3,
                                    repeats: 1,
                                },
                            )
                        } else {
                            TraceFragment::phase_intron(
                                TraceRun {
                                    transition_id: 68,
                                    op: Op::SplitCodon,
                                    query_advance: phase as u32,
                                    target_advance: phase as u32,
                                    repeats: 1,
                                },
                                splice5,
                                loop_run,
                                splice3,
                                TraceRun {
                                    transition_id: 69,
                                    op: Op::SplitCodon,
                                    query_advance: post as u32,
                                    target_advance: post as u32,
                                    repeats: 1,
                                },
                            )
                        };
                        genome_relax(
                            &mut cm[k],
                            &mut pcm[k],
                            add(
                                add(candidate.score, acceptor),
                                protein_score(
                                    translate_dna(&query_codon, 0)[0],
                                    translate_dna(&target.bases[j - 3..j], 0)[0],
                                    scoring,
                                ),
                            ),
                            state,
                            fragment,
                        );
                    }
                }
            }

            // Joint coding introns in phases 0, 1, and 2.
            for phase in 0..3 {
                let post = 3 - phase;
                if i >= min_len + 3 && j >= min_len + 3 {
                    let query_post = i - post;
                    let target_post = j - post;
                    let target_donor = target_post - min_len;
                    if let Some(candidate) = queues.joint_coding[phase][target_donor].best() {
                        joint_coding_targets[phase].insert(candidate);
                    }
                    if target_post > max_len {
                        while joint_coding_targets[phase]
                            .entries
                            .front()
                            .is_some_and(|entry| entry.target_start < target_post - max_len)
                        {
                            joint_coding_targets[phase].entries.pop_front();
                        }
                    }
                    if let (Some(candidate), Some(query_acceptor), Some(target_acceptor)) = (
                        joint_coding_targets[phase].best(),
                        splice_score(
                            &query.bases,
                            query_post - 2,
                            SpliceType::AcceptorForward,
                            intron.force_gtag,
                        ),
                        splice_score(
                            &target.bases,
                            target_post - 2,
                            SpliceType::AcceptorForward,
                            intron.force_gtag,
                        ),
                    ) {
                        if forbidden.is_some_and(|pairs| {
                            pairs.contains(&(
                                candidate.query_start - phase,
                                candidate.target_start - phase,
                            )) || pairs.contains(&(query_post, target_post))
                        }) {
                            continue;
                        }
                        let mut query_codon = Vec::with_capacity(3);
                        query_codon.extend_from_slice(
                            &query.bases[candidate.query_start - phase..candidate.query_start],
                        );
                        query_codon.extend_from_slice(&query.bases[query_post..i]);
                        let mut target_codon = Vec::with_capacity(3);
                        target_codon.extend_from_slice(
                            &target.bases[candidate.target_start - phase..candidate.target_start],
                        );
                        target_codon.extend_from_slice(&target.bases[target_post..j]);
                        let state = match state_from_rank(candidate.state_rank) {
                            State::M => GenomeState::CdsM,
                            State::I => GenomeState::CdsI,
                            State::D => GenomeState::CdsD,
                            State::Stop => unreachable!(),
                        };
                        let splice5 = TraceRun {
                            transition_id: 70,
                            op: Op::Splice5,
                            query_advance: 2,
                            target_advance: 2,
                            repeats: 1,
                        };
                        let loop_run = TraceRun {
                            transition_id: 71,
                            op: Op::Intron,
                            query_advance: (query_post - candidate.query_start - 4) as u32,
                            target_advance: (target_post - candidate.target_start - 4) as u32,
                            repeats: 1,
                        };
                        let splice3 = TraceRun {
                            transition_id: 72,
                            op: Op::Splice3,
                            query_advance: 2,
                            target_advance: 2,
                            repeats: 1,
                        };
                        let fragment = if phase == 0 {
                            TraceFragment::intron_match(
                                splice5,
                                loop_run,
                                splice3,
                                TraceRun {
                                    transition_id: 51,
                                    op: Op::Match,
                                    query_advance: 3,
                                    target_advance: 3,
                                    repeats: 1,
                                },
                            )
                        } else {
                            TraceFragment::phase_intron(
                                TraceRun {
                                    transition_id: 73,
                                    op: Op::SplitCodon,
                                    query_advance: phase as u32,
                                    target_advance: phase as u32,
                                    repeats: 1,
                                },
                                splice5,
                                loop_run,
                                splice3,
                                TraceRun {
                                    transition_id: 74,
                                    op: Op::SplitCodon,
                                    query_advance: post as u32,
                                    target_advance: post as u32,
                                    repeats: 1,
                                },
                            )
                        };
                        genome_relax(
                            &mut cm[k],
                            &mut pcm[k],
                            add(
                                add(add(candidate.score, query_acceptor), target_acceptor),
                                protein_score(
                                    translate_dna(&query_codon, 0)[0],
                                    translate_dna(&target_codon, 0)[0],
                                    scoring,
                                ),
                            ),
                            state,
                            fragment,
                        );
                    }
                }
            }

            // Phase-0 target intron: whole codon after the intron.
            if i >= 3 && j >= 3 + min_len {
                let post_start = j - 3;
                let donor = post_start - min_len;
                let source = idx(i - 3, donor, cols);
                let (value, source_state) = best(&[
                    (cm[source], State::M),
                    (ci[source], State::I),
                    (cd[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor_score) = splice_score(
                        &target.bases,
                        donor,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        coding_phase0.insert(IntronCandidate {
                            start: donor,
                            score: add(value, intron.open_penalty.saturating_add(donor_score)),
                            state_rank: source_state.rank(),
                        });
                    }
                }
                if post_start > max_len {
                    coding_phase0.expire_before(post_start - max_len);
                }
                if let (Some(candidate), Some(acceptor)) = (
                    coding_phase0.best(),
                    splice_score(
                        &target.bases,
                        post_start - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    let state = match state_from_rank(candidate.state_rank) {
                        State::M => GenomeState::CdsM,
                        State::I => GenomeState::CdsI,
                        State::D => GenomeState::CdsD,
                        State::Stop => unreachable!(),
                    };
                    genome_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(
                            add(candidate.score, acceptor),
                            protein_score(
                                translate_dna(&query.bases[i - 3..i], 0)[0],
                                translate_dna(&target.bases[post_start..j], 0)[0],
                                scoring,
                            ),
                        ),
                        state,
                        TraceFragment::intron_match(
                            TraceRun {
                                transition_id: 60,
                                op: Op::Splice5,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 61,
                                op: Op::Intron,
                                query_advance: 0,
                                target_advance: (post_start - candidate.start - 4) as u32,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 62,
                                op: Op::Splice3,
                                query_advance: 0,
                                target_advance: 2,
                                repeats: 1,
                            },
                            TraceRun {
                                transition_id: 51,
                                op: Op::Match,
                                query_advance: 3,
                                target_advance: 3,
                                repeats: 1,
                            },
                        ),
                    );
                }
            }
            // Phase-1 target intron: one base before and two bases after the intron.
            if i >= 3 && j > 2 + min_len {
                let post_start = j - 2;
                let donor = post_start - min_len;
                let source = idx(i - 3, donor - 1, cols);
                let (value, source_state) = best(&[
                    (cm[source], State::M),
                    (ci[source], State::I),
                    (cd[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor_score) = splice_score(
                        &target.bases,
                        donor,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        coding_phase1.insert(IntronCandidate {
                            start: donor,
                            score: add(value, intron.open_penalty.saturating_add(donor_score)),
                            state_rank: source_state.rank(),
                        });
                    }
                }
                if post_start > max_len {
                    coding_phase1.expire_before(post_start - max_len);
                }
                if let (Some(candidate), Some(acceptor)) = (
                    coding_phase1.best(),
                    splice_score(
                        &target.bases,
                        post_start - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    let mut codon = Vec::with_capacity(3);
                    codon.extend_from_slice(&target.bases[candidate.start - 1..candidate.start]);
                    codon.extend_from_slice(&target.bases[post_start..j]);
                    if codon.len() == 3 {
                        let state = match state_from_rank(candidate.state_rank) {
                            State::M => GenomeState::CdsM,
                            State::I => GenomeState::CdsI,
                            State::D => GenomeState::CdsD,
                            State::Stop => unreachable!(),
                        };
                        genome_relax(
                            &mut cm[k],
                            &mut pcm[k],
                            add(
                                add(candidate.score, acceptor),
                                protein_score(
                                    translate_dna(&query.bases[i - 3..i], 0)[0],
                                    translate_dna(&codon, 0)[0],
                                    scoring,
                                ),
                            ),
                            state,
                            TraceFragment::phase_intron(
                                TraceRun {
                                    transition_id: 63,
                                    op: Op::SplitCodon,
                                    query_advance: 1,
                                    target_advance: 1,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 60,
                                    op: Op::Splice5,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 61,
                                    op: Op::Intron,
                                    query_advance: 0,
                                    target_advance: (post_start - candidate.start - 4) as u32,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 62,
                                    op: Op::Splice3,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 64,
                                    op: Op::SplitCodon,
                                    query_advance: 2,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                            ),
                        );
                    }
                }
            }
            // Phase-2 target intron: two bases before and one base after the intron.
            if i >= 3 && j > 1 + min_len + 2 {
                let post_start = j - 1;
                let donor = post_start - min_len;
                let source = idx(i - 3, donor - 2, cols);
                let (value, source_state) = best(&[
                    (cm[source], State::M),
                    (ci[source], State::I),
                    (cd[source], State::D),
                ]);
                if value > 0 {
                    if let Some(donor_score) = splice_score(
                        &target.bases,
                        donor,
                        SpliceType::DonorForward,
                        intron.force_gtag,
                    ) {
                        coding_phase2.insert(IntronCandidate {
                            start: donor,
                            score: add(value, intron.open_penalty.saturating_add(donor_score)),
                            state_rank: source_state.rank(),
                        });
                    }
                }
                if post_start > max_len {
                    coding_phase2.expire_before(post_start - max_len);
                }
                if let (Some(candidate), Some(acceptor)) = (
                    coding_phase2.best(),
                    splice_score(
                        &target.bases,
                        post_start - 2,
                        SpliceType::AcceptorForward,
                        intron.force_gtag,
                    ),
                ) {
                    let mut codon = Vec::with_capacity(3);
                    codon.extend_from_slice(&target.bases[candidate.start - 2..candidate.start]);
                    codon.extend_from_slice(&target.bases[post_start..j]);
                    if codon.len() == 3 {
                        let state = match state_from_rank(candidate.state_rank) {
                            State::M => GenomeState::CdsM,
                            State::I => GenomeState::CdsI,
                            State::D => GenomeState::CdsD,
                            State::Stop => unreachable!(),
                        };
                        genome_relax(
                            &mut cm[k],
                            &mut pcm[k],
                            add(
                                add(candidate.score, acceptor),
                                protein_score(
                                    translate_dna(&query.bases[i - 3..i], 0)[0],
                                    translate_dna(&codon, 0)[0],
                                    scoring,
                                ),
                            ),
                            state,
                            TraceFragment::phase_intron(
                                TraceRun {
                                    transition_id: 63,
                                    op: Op::SplitCodon,
                                    query_advance: 2,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 60,
                                    op: Op::Splice5,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 61,
                                    op: Op::Intron,
                                    query_advance: 0,
                                    target_advance: (post_start - candidate.start - 4) as u32,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 62,
                                    op: Op::Splice3,
                                    query_advance: 0,
                                    target_advance: 2,
                                    repeats: 1,
                                },
                                TraceRun {
                                    transition_id: 64,
                                    op: Op::SplitCodon,
                                    query_advance: 1,
                                    target_advance: 1,
                                    repeats: 1,
                                },
                            ),
                        );
                    }
                }
            }
            if i >= 3 {
                let p = idx(i - 3, j, cols);
                let (candidate, state, transition_id) =
                    if add(cm[p], scoring.codon_gap_open) >= add(ci[p], scoring.codon_gap_extend) {
                        (add(cm[p], scoring.codon_gap_open), GenomeState::CdsM, 52)
                    } else {
                        (add(ci[p], scoring.codon_gap_extend), GenomeState::CdsI, 54)
                    };
                genome_relax(
                    &mut ci[k],
                    &mut pci[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 3,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j >= 3 {
                let p = idx(i, j - 3, cols);
                let (candidate, state, transition_id) =
                    if add(cm[p], scoring.codon_gap_open) >= add(cd[p], scoring.codon_gap_extend) {
                        (add(cm[p], scoring.codon_gap_open), GenomeState::CdsM, 53)
                    } else {
                        (add(cd[p], scoring.codon_gap_extend), GenomeState::CdsD, 55)
                    };
                genome_relax(
                    &mut cd[k],
                    &mut pcd[k],
                    candidate,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 3,
                        repeats: 1,
                    }),
                );
            }
            if i >= 3 && j >= 3 && !forbidden_match_transition(forbidden, i - 3, j - 3, 3, 3) {
                let p = idx(i - 3, j - 3, cols);
                let (candidate, state) =
                    best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                let state = match state {
                    State::M => GenomeState::CdsM,
                    State::I => GenomeState::CdsI,
                    State::D => GenomeState::CdsD,
                    State::Stop => unreachable!(),
                };
                genome_relax(
                    &mut cm[k],
                    &mut pcm[k],
                    add(
                        candidate,
                        protein_score(
                            translate_dna(&query.bases[i - 3..i], 0)[0],
                            translate_dna(&target.bases[j - 3..j], 0)[0],
                            scoring,
                        ),
                    ),
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id: 51,
                        op: Op::Match,
                        query_advance: 3,
                        target_advance: 3,
                        repeats: 1,
                    }),
                );
            }
            for advance in [1_usize, 2, 4, 5] {
                if i >= advance {
                    let p = idx(i - advance, j, cols);
                    let (candidate, state) =
                        best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                    let state = match state {
                        State::M => GenomeState::CdsM,
                        State::I => GenomeState::CdsI,
                        State::D => GenomeState::CdsD,
                        State::Stop => unreachable!(),
                    };
                    genome_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(candidate, scoring.frameshift),
                        state,
                        TraceFragment::one(TraceRun {
                            transition_id: 56,
                            op: Op::Frameshift,
                            query_advance: advance as u32,
                            target_advance: 0,
                            repeats: 1,
                        }),
                    );
                }
                if j >= advance {
                    let p = idx(i, j - advance, cols);
                    let (candidate, state) =
                        best(&[(cm[p], State::M), (ci[p], State::I), (cd[p], State::D)]);
                    let state = match state {
                        State::M => GenomeState::CdsM,
                        State::I => GenomeState::CdsI,
                        State::D => GenomeState::CdsD,
                        State::Stop => unreachable!(),
                    };
                    genome_relax(
                        &mut cm[k],
                        &mut pcm[k],
                        add(candidate, scoring.frameshift),
                        state,
                        TraceFragment::one(TraceRun {
                            transition_id: 57,
                            op: Op::Frameshift,
                            query_advance: 0,
                            target_advance: advance as u32,
                            repeats: 1,
                        }),
                    );
                }
            }
            // The composed cDNA graph permits its CDS epsilon exit from every
            // coding state.  Keeping the predecessor state is significant
            // after a codon gap or frameshift at the boundary.
            let (candidate, state) =
                best(&[(cm[k], State::M), (ci[k], State::I), (cd[k], State::D)]);
            genome_relax(
                &mut u3m[k],
                &mut pu3m[k],
                candidate,
                match state {
                    State::M => GenomeState::CdsM,
                    State::I => GenomeState::CdsI,
                    State::D => GenomeState::CdsD,
                    State::Stop => unreachable!(),
                },
                TraceFragment::empty(),
            );
            if i > 0 && j > 0 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 1))) {
                let source = idx(i - 1, j - 1, cols);
                let (value, state) = best(&[
                    (u3m[source], State::M),
                    (u3i[source], State::I),
                    (u3d[source], State::D),
                ]);
                genome_relax(
                    &mut u3m[k],
                    &mut pu3m[k],
                    add(
                        value,
                        dna_score(query.bases[i - 1], target.bases[j - 1], scoring),
                    ),
                    match state {
                        State::M => GenomeState::U3M,
                        State::I => GenomeState::U3I,
                        State::D => GenomeState::U3D,
                        State::Stop => unreachable!(),
                    },
                    TraceFragment::one(TraceRun {
                        transition_id: 30,
                        op: Op::Match,
                        query_advance: 1,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            if i > 0 {
                let source = idx(i - 1, j, cols);
                let (value, state, transition_id) =
                    if add(u3m[source], scoring.gap_open) >= add(u3i[source], scoring.gap_extend) {
                        (add(u3m[source], scoring.gap_open), GenomeState::U3M, 31)
                    } else {
                        (add(u3i[source], scoring.gap_extend), GenomeState::U3I, 32)
                    };
                genome_relax(
                    &mut u3i[k],
                    &mut pu3i[k],
                    value,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                );
            }
            if j > 0 {
                let source = idx(i, j - 1, cols);
                let (value, state, transition_id) =
                    if add(u3m[source], scoring.gap_open) >= add(u3d[source], scoring.gap_extend) {
                        (add(u3m[source], scoring.gap_open), GenomeState::U3M, 33)
                    } else {
                        (add(u3d[source], scoring.gap_extend), GenomeState::U3D, 34)
                    };
                genome_relax(
                    &mut u3d[k],
                    &mut pu3d[k],
                    value,
                    state,
                    TraceFragment::one(TraceRun {
                        transition_id,
                        op: Op::Delete,
                        query_advance: 0,
                        target_advance: 1,
                        repeats: 1,
                    }),
                );
            }
            for (value, state) in [
                (mat[k], GenomeState::UtrM),
                (ins[k], GenomeState::UtrI),
                (del[k], GenomeState::UtrD),
                (cm[k], GenomeState::CdsM),
                (ci[k], GenomeState::CdsI),
                (cd[k], GenomeState::CdsD),
                (u3m[k], GenomeState::U3M),
                (u3i[k], GenomeState::U3I),
                (u3d[k], GenomeState::U3D),
            ] {
                // Upstream retains the first endpoint at an equal score.  In
                // particular this resolves reverse-orientation local paths
                // before the later, shifted CDS/UTR equivalent.
                if value > end.3 {
                    end = (i, j, state, value);
                }
            }
        }
        if let Some(checkpoints) = checkpoint_sink.as_deref_mut() {
            if i % checkpoint_stride.max(1) == 0 || i == n {
                let first_row = i + 1 - genome_checkpoint_history_rows(intron).min(i + 1);
                let history = (first_row..=i)
                    .map(|row| GenomeScoreRow {
                        row,
                        u5m: mat.row_clone(row),
                        u5i: ins.row_clone(row),
                        u5d: del.row_clone(row),
                        cm: cm.row_clone(row),
                        ci: ci.row_clone(row),
                        cd: cd.row_clone(row),
                        u3m: u3m.row_clone(row),
                        u3i: u3i.row_clone(row),
                        u3d: u3d.row_clone(row),
                    })
                    .collect();
                checkpoints.push(GenomeCheckpoint {
                    row: i,
                    history,
                    queues: queues.clone(),
                });
            }
        }
    }
    let (end_i, end_j, final_state, score) = if let Some((i, j, state)) = forced_end {
        debug_assert!(i <= n && j <= m);
        let k = idx(i, j, cols);
        let score = match state {
            GenomeState::UtrM => mat[k],
            GenomeState::UtrI => ins[k],
            GenomeState::UtrD => del[k],
            GenomeState::CdsM => cm[k],
            GenomeState::CdsI => ci[k],
            GenomeState::CdsD => cd[k],
            GenomeState::U3M => u3m[k],
            GenomeState::U3I => u3i[k],
            GenomeState::U3D => u3d[k],
        };
        (i, j, state, score)
    } else {
        end
    };
    let (query_end, target_end) = (end_i as u64, end_j as u64);
    let (i, j, start_state, fragments) = genome_traceback_fragments(
        end_i,
        end_j,
        final_state,
        cols,
        &pmat,
        &pins,
        &pdel,
        &pcm,
        &pci,
        &pcd,
        &pu3m,
        &pu3i,
        &pu3d,
        query,
        target,
        scoring,
        intron,
    );
    if let Some(sink) = traceback_start_sink {
        *sink = (i, j, start_state);
    }
    let mut trace = Vec::new();
    let mut raw_trace = vec![RawStep {
        transition_id: 29,
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
        GenomeState::UtrI => Some(35),
        GenomeState::UtrD => Some(36),
        GenomeState::CdsI => Some(59),
        GenomeState::CdsD => Some(79),
        GenomeState::U3I => Some(35),
        GenomeState::U3D => Some(36),
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
    raw_trace.push(RawStep {
        transition_id: if genome_state_is_coding(final_state) {
            89
        } else {
            39
        },
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    raw_trace.push(RawStep {
        transition_id: 90,
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
        score,
        raw_trace,
        trace,
    }
}

pub fn align_genome_to_genome(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Alignment {
    align_genome_to_genome_forbidden(
        query, target, scoring, intron, None, false, None, None, None, None, 1,
    )
}

/// Exact low-memory traceback for genome-to-genome.  The forward endpoint is
/// found with rolling scores and no resident parent matrix; traceback then
/// recomputes one parent-row block at a time.  Query/joint intron queues are
/// rebuilt from the prefix for every block, which is slower than restoring
/// the saved queue checkpoints but keeps memory bounded and preserves the
/// exact recurrence while that faster restoration path is completed.
fn align_genome_to_genome_checkpointed(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    budget_bytes: Option<u128>,
) -> Alignment {
    let mut endpoint_state = (0, 0, GenomeState::UtrM);
    let endpoint = align_genome_to_genome_forbidden(
        query,
        target,
        scoring,
        intron,
        None,
        true,
        Some((0, 0)),
        None,
        Some(&mut endpoint_state),
        None,
        1,
    );
    let final_state = endpoint_state.2;
    let mut cursor = (
        endpoint.query_end as usize,
        endpoint.target_end as usize,
        final_state,
    );
    let block_rows = budget_bytes
        .map(|budget| {
            genome_checkpoint_plan(query.bases.len(), target.bases.len(), intron, Some(budget))
                .stride
        })
        .unwrap_or(64)
        .max(genome_checkpoint_history_rows(intron))
        .min(query.bases.len().max(1));
    let mut trace_segments = Vec::new();
    let mut raw_segments = Vec::new();
    loop {
        let first_row = cursor.0.saturating_sub(block_rows.saturating_sub(1));
        let prefix = Sequence {
            id: query.id.clone(),
            bases: query.bases[..cursor.0].to_vec(),
        };
        let mut start = cursor;
        let segment = align_genome_to_genome_forbidden(
            &prefix,
            target,
            scoring,
            intron,
            None,
            true,
            Some((first_row, cursor.0)),
            Some(cursor),
            Some(&mut start),
            None,
            1,
        );
        let terminal_len = 2 + usize::from(genome_terminal_transition(cursor.2).is_some());
        debug_assert!(segment.raw_trace.len() > terminal_len);
        trace_segments.push(segment.trace);
        raw_segments.push(segment.raw_trace[1..segment.raw_trace.len() - terminal_len].to_vec());
        if start == cursor {
            break;
        }
        cursor = start;
    }
    trace_segments.reverse();
    raw_segments.reverse();
    let mut trace = Vec::new();
    for segment in trace_segments {
        for run in segment {
            TraceFragment::one(run).append_to(&mut trace);
        }
    }
    let mut raw_trace = vec![RawStep {
        transition_id: 29,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    raw_trace.extend(raw_segments.into_iter().flatten());
    if let Some(transition_id) = genome_terminal_transition(final_state) {
        raw_trace.push(RawStep {
            transition_id,
            query_advance: 0,
            target_advance: 0,
            score: 0,
        });
    }
    raw_trace.push(RawStep {
        transition_id: if genome_state_is_coding(final_state) {
            89
        } else {
            39
        },
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    raw_trace.push(RawStep {
        transition_id: 90,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: cursor.0 as u64,
        query_end: endpoint.query_end,
        query_strand: Strand::Forward,
        target_start: cursor.1 as u64,
        target_end: endpoint.target_end,
        target_len: target.bases.len() as u64,
        target_strand: Strand::Forward,
        score: endpoint.score,
        raw_trace,
        trace,
    }
}

fn genome_terminal_transition(state: GenomeState) -> Option<u16> {
    match state {
        GenomeState::UtrI | GenomeState::U3I => Some(35),
        GenomeState::UtrD | GenomeState::U3D => Some(36),
        GenomeState::CdsI => Some(59),
        GenomeState::CdsD => Some(79),
        GenomeState::UtrM | GenomeState::CdsM | GenomeState::U3M => None,
    }
}

fn align_genome_to_genome_memory_aware(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    rolling_scores: bool,
    budget_bytes: Option<u128>,
) -> Alignment {
    if rolling_scores {
        align_genome_to_genome_checkpointed(query, target, scoring, intron, budget_bytes)
    } else {
        align_genome_to_genome_forbidden(
            query, target, scoring, intron, None, false, None, None, None, None, 1,
        )
    }
}

fn full_genome2genome_dp_bytes(query_len: usize, target_len: usize) -> u128 {
    let cells = (query_len as u128 + 1) * (target_len as u128 + 1);
    cells * 9 * (size_of::<Score>() as u128 + size_of::<Option<GenomeParent>>() as u128)
}

fn genome2genome_uses_low_memory(query_len: usize, target_len: usize, memory_mb: usize) -> bool {
    full_genome2genome_dp_bytes(query_len, target_len) > ((memory_mb as u128) << 20)
}

pub fn align_genome_to_genome_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
) -> Vec<Alignment> {
    enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_genome_to_genome_forbidden(
            query,
            target,
            scoring,
            intron,
            Some(forbidden),
            false,
            None,
            None,
            None,
            None,
            1,
        )
    })
}

pub fn align_genome_to_genome_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_genome_to_genome_suboptimal(
                query, target, scoring, intron, threshold,
            ));
        }
    }
    out
}

/// Suboptimal genome-to-genome paths over all upstream orientation pairs.
pub fn align_genome_to_genome_database_suboptimal_stranded(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            for alignment in
                align_genome_to_genome_suboptimal(query, target, scoring, intron, threshold)
            {
                out.push(alignment);
            }
            if !both_strands {
                continue;
            }
            let reverse_query = Sequence {
                id: query.id.clone(),
                bases: reverse_complement(&query.bases),
            };
            let reverse_target = Sequence {
                id: target.id.clone(),
                bases: reverse_complement(&target.bases),
            };
            for (oriented_query, query_reverse) in [(query, false), (&reverse_query, true)] {
                for (oriented_target, target_reverse) in [(target, false), (&reverse_target, true)]
                {
                    if !query_reverse && !target_reverse {
                        continue;
                    }
                    for mut alignment in align_genome_to_genome_suboptimal(
                        oriented_query,
                        oriented_target,
                        scoring,
                        intron,
                        threshold,
                    ) {
                        alignment.target_id = target.id.clone();
                        if target_reverse {
                            alignment.target_start =
                                target.bases.len() as u64 - alignment.target_start;
                            alignment.target_end = target.bases.len() as u64 - alignment.target_end;
                            alignment.target_strand = Strand::Reverse;
                        }
                        if query_reverse {
                            alignment.query_start =
                                query.bases.len() as u64 - alignment.query_start;
                            alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                            alignment.query_strand = Strand::Reverse;
                        }
                        out.push(alignment);
                    }
                }
            }
        }
    }
    out
}

pub fn align_genome_to_genome_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<Alignment> {
    queries
        .iter()
        .flat_map(|query| {
            targets
                .iter()
                .map(move |target| align_genome_to_genome(query, target, scoring, intron))
        })
        .collect()
}

/// Genome-to-genome database alignment over the upstream four orientation
/// combinations when `both_strands` is enabled.
pub fn align_genome_to_genome_database_stranded(
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
            let rolling_scores =
                genome2genome_uses_low_memory(query.bases.len(), target.bases.len(), memory_mb);
            out.push(align_genome_to_genome_memory_aware(
                query,
                target,
                scoring,
                intron,
                rolling_scores,
                rolling_scores.then_some((memory_mb as u128) << 20),
            ));
            if !both_strands {
                continue;
            }
            let reverse_query = Sequence {
                id: query.id.clone(),
                bases: reverse_complement(&query.bases),
            };
            let reverse_target = Sequence {
                id: target.id.clone(),
                bases: reverse_complement(&target.bases),
            };
            for (oriented_query, query_reverse) in [(query, false), (&reverse_query, true)] {
                for (oriented_target, target_reverse) in [(target, false), (&reverse_target, true)]
                {
                    if !query_reverse && !target_reverse {
                        continue;
                    }
                    let mut alignment = align_genome_to_genome_memory_aware(
                        oriented_query,
                        oriented_target,
                        scoring,
                        intron,
                        rolling_scores,
                        rolling_scores.then_some((memory_mb as u128) << 20),
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

/// Local coding-DNA to coding-DNA Viterbi with codon affine gaps and
/// independent query/target 1, 2, 4, and 5 nt frameshifts.
pub fn coding2coding_score(query: &Sequence, target: &Sequence, scoring: Scoring) -> Score {
    align_coding2coding(query, target, scoring).score
}

/// Score a local protein-to-genome alignment with target introns in phases
/// 0, 1, and 2.  The split-codon paths follow upstream `Phase_create`:
/// 0:0, 1:2, and 2:1.
pub fn protein2genome_score(query: &Sequence, target: &Sequence, intron: IntronScoring) -> Score {
    protein2genome_score_with_scoring(query, target, Scoring::default(), intron)
}

pub fn align_protein_to_genome(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = vec![protein2genome_alignment_oriented(
        query,
        target,
        &target.bases,
        scoring,
        intron,
        Strand::Forward,
        false,
        None,
    )];
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.push(protein2genome_alignment_oriented(
            query,
            target,
            &reverse,
            scoring,
            intron,
            Strand::Reverse,
            false,
            None,
        ));
    }
    out
}

pub fn align_protein_to_genome_bestfit(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = vec![protein2genome_alignment_oriented(
        query,
        target,
        &target.bases,
        scoring,
        intron,
        Strand::Forward,
        true,
        None,
    )];
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.push(protein2genome_alignment_oriented(
            query,
            target,
            &reverse,
            scoring,
            intron,
            Strand::Reverse,
            true,
            None,
        ));
    }
    out
}

/// Configurable local protein-to-genome score Viterbi.
pub fn protein2genome_score_with_scoring(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Score {
    align_protein_to_genome(query, target, scoring, intron, true)
        .into_iter()
        .map(|alignment| alignment.score)
        .max()
        .unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn protein2genome_raw_fragment(
    parent: P2GenomeParent,
    destination: State,
    query_end: usize,
    target_end: usize,
    query: &Sequence,
    target: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<RawStep> {
    let (_, target_advance) = parent.fragment.advances();
    let target_start = target_end - target_advance;
    let mut raw = Vec::new();
    if destination == State::M {
        let epsilon = match parent.state {
            State::I => Some(6),
            State::D => Some(7),
            _ => None,
        };
        if let Some(transition_id) = epsilon {
            raw.push(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
        }
    }
    let atoms = parent
        .fragment
        .atoms
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if atoms.iter().any(|run| run.op == Op::Intron) {
        let phase = atoms
            .iter()
            .find(|run| run.op == Op::SplitCodon)
            .map_or(0, |run| run.target_advance);
        let post_len = 3 - phase;
        let donor = target_start + phase as usize;
        let post_start = target_end - post_len as usize;
        if phase > 0 {
            raw.push(RawStep {
                transition_id: 100 + phase as u16,
                query_advance: 0,
                target_advance: phase,
                score: 0,
            });
        }
        raw.push(RawStep {
            transition_id: 110 + phase as u16,
            query_advance: 0,
            target_advance: 2,
            score: intron.open_penalty.saturating_add(
                splice_score(target, donor, SpliceType::DonorForward, intron.force_gtag)
                    .unwrap_or(0),
            ),
        });
        let loop_len = post_start - donor - 4;
        for _ in 0..loop_len {
            raw.push(RawStep {
                transition_id: 120 + phase as u16,
                query_advance: 0,
                target_advance: 1,
                score: 0,
            });
        }
        raw.push(RawStep {
            transition_id: 130 + phase as u16,
            query_advance: 0,
            target_advance: 2,
            score: splice_score(
                target,
                post_start - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            )
            .unwrap_or(0),
        });
        let mut codon = Vec::with_capacity(3);
        codon.extend_from_slice(&target[target_start..donor]);
        codon.extend_from_slice(&target[post_start..target_end]);
        raw.push(RawStep {
            transition_id: 140 + phase as u16,
            query_advance: 1,
            target_advance: post_len,
            score: protein_score(
                query.bases[query_end - 1],
                translate_dna(&codon, 0)[0],
                scoring,
            ),
        });
        return raw;
    }
    let run = atoms[0];
    match run.op {
        Op::Match => raw.push(RawStep {
            transition_id: 1,
            query_advance: 1,
            target_advance: 3,
            score: protein_score(
                query.bases[query_end - 1],
                translate_dna(&target[target_end - 3..target_end], 0)[0],
                scoring,
            ),
        }),
        Op::Insert | Op::Delete => raw.push(RawStep {
            transition_id: run.transition_id,
            query_advance: run.query_advance,
            target_advance: run.target_advance,
            score: if matches!(run.transition_id, 2 | 3) {
                scoring.codon_gap_open
            } else {
                scoring.codon_gap_extend
            },
        }),
        Op::Frameshift => {
            let short = run.target_advance % 3;
            raw.push(RawStep {
                transition_id: if short == 1 { 8 } else { 9 },
                query_advance: 0,
                target_advance: short,
                score: scoring.frameshift,
            });
            raw.push(RawStep {
                transition_id: if run.target_advance > 3 { 11 } else { 10 },
                query_advance: 0,
                target_advance: if run.target_advance > 3 { 3 } else { 0 },
                score: 0,
            });
        }
        _ => unreachable!("protein2genome parent fragment"),
    }
    raw
}

#[allow(clippy::too_many_arguments)]
fn coding2genome_raw_fragment(
    parent: P2GenomeParent,
    destination: State,
    query_end: usize,
    target_end: usize,
    query: &Sequence,
    target: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
) -> Vec<RawStep> {
    let (query_advance, target_advance) = parent.fragment.advances();
    let query_start = query_end - query_advance;
    let target_start = target_end - target_advance;
    let mut raw = Vec::new();
    if destination == State::M {
        let epsilon = match parent.state {
            State::I => Some(4),
            State::D => Some(7),
            _ => None,
        };
        if let Some(transition_id) = epsilon {
            raw.push(RawStep {
                transition_id,
                query_advance: 0,
                target_advance: 0,
                score: 0,
            });
        }
    }
    let atoms = parent
        .fragment
        .atoms
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if atoms.iter().any(|run| run.op == Op::Intron) {
        let phase = atoms
            .iter()
            .find(|run| run.op == Op::SplitCodon)
            .map_or(0, |run| run.target_advance);
        let post_len = 3 - phase;
        let donor = target_start + phase as usize;
        let post_start = target_end - post_len as usize;
        if phase > 0 {
            raw.push(RawStep {
                transition_id: 100 + phase as u16,
                query_advance: phase,
                target_advance: phase,
                score: 0,
            });
        }
        raw.push(RawStep {
            transition_id: 110 + phase as u16,
            query_advance: 0,
            target_advance: 2,
            score: intron.open_penalty.saturating_add(
                splice_score(target, donor, SpliceType::DonorForward, intron.force_gtag)
                    .unwrap_or(0),
            ),
        });
        let loop_len = post_start - donor - 4;
        for _ in 0..loop_len {
            raw.push(RawStep {
                transition_id: 120 + phase as u16,
                query_advance: 0,
                target_advance: 1,
                score: 0,
            });
        }
        raw.push(RawStep {
            transition_id: 130 + phase as u16,
            query_advance: 0,
            target_advance: 2,
            score: splice_score(
                target,
                post_start - 2,
                SpliceType::AcceptorForward,
                intron.force_gtag,
            )
            .unwrap_or(0),
        });
        let mut target_codon = Vec::with_capacity(3);
        target_codon.extend_from_slice(&target[target_start..donor]);
        target_codon.extend_from_slice(&target[post_start..target_end]);
        raw.push(RawStep {
            transition_id: 140 + phase as u16,
            query_advance: post_len,
            target_advance: post_len,
            score: protein_score(
                translate_dna(&query.bases[query_start..query_end], 0)[0],
                translate_dna(&target_codon, 0)[0],
                scoring,
            ),
        });
        return raw;
    }
    let run = atoms[0];
    match run.op {
        Op::Match => raw.push(RawStep {
            transition_id: 1,
            query_advance: 3,
            target_advance: 3,
            score: protein_score(
                translate_dna(&query.bases[query_end - 3..query_end], 0)[0],
                translate_dna(&target[target_end - 3..target_end], 0)[0],
                scoring,
            ),
        }),
        Op::Insert => raw.push(RawStep {
            transition_id: if run.transition_id == 2 { 2 } else { 3 },
            query_advance: 3,
            target_advance: 0,
            score: if run.transition_id == 2 {
                scoring.codon_gap_open
            } else {
                scoring.codon_gap_extend
            },
        }),
        Op::Delete => raw.push(RawStep {
            transition_id: if run.transition_id == 3 { 5 } else { 6 },
            query_advance: 0,
            target_advance: 3,
            score: if run.transition_id == 3 {
                scoring.codon_gap_open
            } else {
                scoring.codon_gap_extend
            },
        }),
        Op::Frameshift => {
            if run.query_advance > 0 {
                let short = run.query_advance % 3;
                raw.push(RawStep {
                    transition_id: if short == 1 { 8 } else { 9 },
                    query_advance: short,
                    target_advance: 0,
                    score: scoring.frameshift,
                });
                raw.push(RawStep {
                    transition_id: if run.query_advance > 3 { 11 } else { 10 },
                    query_advance: if run.query_advance > 3 { 3 } else { 0 },
                    target_advance: 0,
                    score: 0,
                });
            } else {
                let short = run.target_advance % 3;
                raw.push(RawStep {
                    transition_id: if short == 1 { 12 } else { 13 },
                    query_advance: 0,
                    target_advance: short,
                    score: scoring.frameshift,
                });
                raw.push(RawStep {
                    transition_id: if run.target_advance > 3 { 15 } else { 14 },
                    query_advance: 0,
                    target_advance: if run.target_advance > 3 { 3 } else { 0 },
                    score: 0,
                });
            }
        }
        _ => unreachable!("coding2genome parent fragment"),
    }
    raw
}

#[allow(clippy::too_many_arguments)]
fn protein2genome_alignment_oriented(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    bestfit: bool,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let (n, m, cols) = (query.bases.len(), bases.len(), bases.len() + 1);
    let size = (n + 1) * (m + 1);
    let mut mm = vec![0; size];
    if bestfit {
        mm.fill(NEG_INF);
        for j in 0..=m {
            mm[idx(0, j, cols)] = 0;
        }
    }
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut phase_parents: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pm: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pi: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pd: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut end = (0usize, 0usize, State::M, if bestfit { NEG_INF } else { 0 });
    for i in 1..=n {
        let mut phase_windows = [
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
        ];
        for j in 0..=m {
            let k = idx(i, j, cols);
            for (phase, window) in phase_windows.iter_mut().enumerate() {
                let post_len = 3 - phase;
                if j < post_len + intron.min_len as usize + phase {
                    continue;
                }
                let post_start = j - post_len;
                let donor = post_start - intron.min_len as usize;
                let source_j = donor - phase;
                let source_k = idx(i - 1, source_j, cols);
                let (source, state) = best(&[
                    (mm[source_k], State::M),
                    (ii[source_k], State::I),
                    (dd[source_k], State::D),
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
                if forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, post_start))) {
                    continue;
                }
                let mut codon = Vec::with_capacity(3);
                codon.extend_from_slice(&bases[pre_start..candidate.start]);
                codon.extend_from_slice(&bases[post_start..j]);
                if codon.len() != 3 {
                    continue;
                }
                let aa = translate_dna(&codon, 0)[0];
                let score = add(
                    add(candidate.score, acceptor),
                    protein_score(query.bases[i - 1], aa, scoring),
                );
                if score > mm[k] {
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
                    phase_parents[k] = Some(P2GenomeParent {
                        state: state_from_rank(candidate.state_rank),
                        fragment,
                    });
                    pm[k] = phase_parents[k];
                    mm[k] = score;
                }
            }
            if j >= 3 {
                let p = idx(i, j - 3, cols);
                let from_m = add(mm[p], scoring.codon_gap_open);
                let from_d = add(dd[p], scoring.codon_gap_extend);
                if from_m >= from_d && (from_m >= 0 || bestfit) {
                    dd[k] = from_m;
                    pd[k] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 3,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                } else if from_d >= 0 || bestfit {
                    dd[k] = from_d;
                    pd[k] = Some(P2GenomeParent {
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
            for advance in [1_usize, 2, 4, 5] {
                if j < advance {
                    continue;
                }
                let p = idx(i, j - advance, cols);
                let (source, state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let value = add(source, scoring.frameshift);
                if value > mm[k] {
                    mm[k] = value;
                    pm[k] = Some(P2GenomeParent {
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
            if j >= 3 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 1, j - 3))) {
                let p = idx(i - 1, j - 3, cols);
                let aa = translate_dna(&bases[j - 3..j], 0)[0];
                let sub = protein_score(query.bases[i - 1], aa, scoring);
                let (source, state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let value = add(source, sub);
                if value >= mm[k] {
                    mm[k] = value;
                    pm[k] = Some(P2GenomeParent {
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
            let p = idx(i - 1, j, cols);
            let from_m = add(mm[p], scoring.codon_gap_open);
            let from_i = add(ii[p], scoring.codon_gap_extend);
            if from_m >= from_i && (from_m >= 0 || bestfit) {
                ii[k] = from_m;
                pi[k] = Some(P2GenomeParent {
                    state: State::M,
                    fragment: TraceFragment::one(TraceRun {
                        transition_id: 2,
                        op: Op::Insert,
                        query_advance: 1,
                        target_advance: 0,
                        repeats: 1,
                    }),
                });
            } else if from_i >= 0 || bestfit {
                ii[k] = from_i;
                pi[k] = Some(P2GenomeParent {
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
            for (value, state) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if (!bestfit || i == n)
                    && (value, i, j, state.rank()) > (end.3, end.0, end.1, end.2.rank())
                {
                    end = (i, j, state, value);
                }
            }
        }
    }
    debug_assert!(
        phase_parents
            .iter()
            .chain(pm.iter())
            .chain(pi.iter())
            .chain(pd.iter())
            .flatten()
            .all(|parent| parent.state != State::Stop && parent.fragment.atoms[0].is_some())
    );
    let (mut i, mut j, mut state, score) = end;
    let final_state = state;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop {
        let parent = match state {
            State::M => pm[idx(i, j, cols)],
            State::I => pi[idx(i, j, cols)],
            State::D => pd[idx(i, j, cols)],
            State::Stop => None,
        };
        let Some(parent) = parent else { break };
        let raw_fragment =
            protein2genome_raw_fragment(parent, state, i, j, query, bases, scoring, intron);
        let (query_advance, target_advance) =
            parent
                .fragment
                .atoms
                .iter()
                .flatten()
                .fold((0usize, 0usize), |(q, t), run| {
                    (
                        q + run.query_advance as usize * run.repeats as usize,
                        t + run.target_advance as usize * run.repeats as usize,
                    )
                });
        debug_assert!(i >= query_advance && j >= target_advance);
        i -= query_advance;
        j -= target_advance;
        fragments.push((parent.fragment, raw_fragment));
        state = parent.state;
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
        score,
        raw_trace,
        trace,
    }
}

pub fn align_protein_to_genome_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
    bestfit: bool,
) -> Vec<Alignment> {
    let mut out = enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        protein2genome_alignment_oriented(
            query,
            target,
            &target.bases,
            scoring,
            intron,
            Strand::Forward,
            bestfit,
            Some(forbidden),
        )
    });
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.extend(enumerate_suboptimal_pair(
            query,
            target,
            threshold,
            |forbidden| {
                protein2genome_alignment_oriented(
                    query,
                    target,
                    &reverse,
                    scoring,
                    intron,
                    Strand::Reverse,
                    bestfit,
                    Some(forbidden),
                )
            },
        ));
    }
    out
}

pub fn align_protein_to_genome_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
    bestfit: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_to_genome_suboptimal(
                query,
                target,
                scoring,
                intron,
                threshold,
                both_strands,
                bestfit,
            ));
        }
    }
    out
}

pub fn align_protein_to_genome_bestfit_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_to_genome_bestfit(
                query,
                target,
                scoring,
                intron,
                both_strands,
            ));
        }
    }
    out
}

/// Align every protein/genome pair, emitting one alignment per requested target strand.
pub fn align_protein_to_genome_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_to_genome(
                query,
                target,
                scoring,
                intron,
                both_strands,
            ));
        }
    }
    out
}

/// Score a local coding-DNA-to-genome alignment with target introns in phases 0, 1, and 2.
pub fn coding2genome_score(query: &Sequence, target: &Sequence, intron: IntronScoring) -> Score {
    coding2genome_score_with_scoring(query, target, Scoring::default(), intron)
}

pub fn align_coding_to_genome(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = vec![coding2genome_alignment_oriented(
        query,
        target,
        &target.bases,
        scoring,
        intron,
        Strand::Forward,
        None,
    )];
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.push(coding2genome_alignment_oriented(
            query,
            target,
            &reverse,
            scoring,
            intron,
            Strand::Reverse,
            None,
        ));

        // Coding DNA is itself strand-specific.  Upstream searches both
        // query orientations as well as both target orientations; remap the
        // reverse-query coordinates back onto the caller's forward sequence.
        let reverse_query = Sequence {
            id: query.id.clone(),
            bases: reverse_complement(&query.bases),
        };
        for (bases, strand) in [
            (&target.bases, Strand::Forward),
            (&reverse, Strand::Reverse),
        ] {
            let mut alignment = coding2genome_alignment_oriented(
                &reverse_query,
                target,
                bases,
                scoring,
                intron,
                strand,
                None,
            );
            alignment.query_start = query.bases.len() as u64 - alignment.query_start;
            alignment.query_end = query.bases.len() as u64 - alignment.query_end;
            alignment.query_strand = Strand::Reverse;
            out.push(alignment);
        }
    }
    out
}

/// Configurable local coding-DNA-to-genome score Viterbi.
pub fn coding2genome_score_with_scoring(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
) -> Score {
    align_coding_to_genome(query, target, scoring, intron, true)
        .into_iter()
        .map(|alignment| alignment.score)
        .max()
        .unwrap_or(0)
}

fn coding2genome_alignment_oriented(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    scoring: Scoring,
    intron: IntronScoring,
    strand: Strand,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let intron = canonical_intron_scoring(intron);
    let (n, m, cols) = (query.bases.len(), bases.len(), bases.len() + 1);
    let size = (n + 1) * (m + 1);
    let mut mm = vec![0; size];
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut phase_parents: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pm: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pi: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut pd: Vec<Option<P2GenomeParent>> = vec![None; size];
    let mut end = (0usize, 0usize, State::M, 0);
    for i in 0..=n {
        let mut phase_windows = [
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
            IntronCandidateWindow::default(),
        ];
        for j in 0..=m {
            let k = idx(i, j, cols);
            if i >= 3 {
                for (phase, window) in phase_windows.iter_mut().enumerate() {
                    let post_len = 3 - phase;
                    if j < post_len + intron.min_len as usize + phase {
                        continue;
                    }
                    let post_start = j - post_len;
                    let donor = post_start - intron.min_len as usize;
                    let source_j = donor - phase;
                    let source_k = idx(i - 3, source_j, cols);
                    let (source, state) = best(&[
                        (mm[source_k], State::M),
                        (ii[source_k], State::I),
                        (dd[source_k], State::D),
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
                    if forbidden.is_some_and(|pairs| {
                        pairs.contains(&(i - 3, pre_start))
                            || pairs.contains(&(i - post_len, post_start))
                    }) {
                        continue;
                    }
                    let mut codon = Vec::with_capacity(3);
                    codon.extend_from_slice(&bases[pre_start..candidate.start]);
                    codon.extend_from_slice(&bases[post_start..j]);
                    if codon.len() != 3 {
                        continue;
                    }
                    let aa = translate_dna(&codon, 0)[0];
                    let score = add(
                        add(candidate.score, acceptor),
                        protein_score(translate_dna(&query.bases[i - 3..i], 0)[0], aa, scoring),
                    );
                    if score > mm[k] {
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
                        phase_parents[k] = Some(P2GenomeParent {
                            state: state_from_rank(candidate.state_rank),
                            fragment,
                        });
                        pm[k] = phase_parents[k];
                        mm[k] = score;
                    }
                }
            }
            if j >= 3 {
                let p = idx(i, j - 3, cols);
                let from_m = add(mm[p], scoring.codon_gap_open);
                let from_d = add(dd[p], scoring.codon_gap_extend);
                if from_m >= from_d && from_m >= 0 {
                    dd[k] = from_m;
                    pd[k] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 3,
                            op: Op::Delete,
                            query_advance: 0,
                            target_advance: 3,
                            repeats: 1,
                        }),
                    });
                } else if from_d >= 0 {
                    dd[k] = from_d;
                    pd[k] = Some(P2GenomeParent {
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
            for advance in [1_usize, 2, 4, 5] {
                if i >= advance {
                    let p = idx(i - advance, j, cols);
                    let (source, state) =
                        best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                    let value = add(source, scoring.frameshift);
                    if value > mm[k] {
                        mm[k] = value;
                        pm[k] = Some(P2GenomeParent {
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
                if j < advance {
                    continue;
                }
                let p = idx(i, j - advance, cols);
                let (source, state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let value = add(source, scoring.frameshift);
                if value > mm[k] {
                    mm[k] = value;
                    pm[k] = Some(P2GenomeParent {
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
            if i >= 3 && j >= 3 && !forbidden.is_some_and(|pairs| pairs.contains(&(i - 3, j - 3))) {
                let p = idx(i - 3, j - 3, cols);
                let aa = translate_dna(&bases[j - 3..j], 0)[0];
                let sub = protein_score(translate_dna(&query.bases[i - 3..i], 0)[0], aa, scoring);
                let (source, state) =
                    best(&[(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]);
                let value = add(source, sub);
                if value >= mm[k] {
                    mm[k] = value;
                    pm[k] = Some(P2GenomeParent {
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
            if i >= 3 {
                let p = idx(i - 3, j, cols);
                let from_m = add(mm[p], scoring.codon_gap_open);
                let from_i = add(ii[p], scoring.codon_gap_extend);
                if from_m >= from_i && from_m >= 0 {
                    ii[k] = from_m;
                    pi[k] = Some(P2GenomeParent {
                        state: State::M,
                        fragment: TraceFragment::one(TraceRun {
                            transition_id: 2,
                            op: Op::Insert,
                            query_advance: 3,
                            target_advance: 0,
                            repeats: 1,
                        }),
                    });
                } else if from_i >= 0 {
                    ii[k] = from_i;
                    pi[k] = Some(P2GenomeParent {
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
            for (value, state) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if (value, i, j, state.rank()) > (end.3, end.0, end.1, end.2.rank()) {
                    end = (i, j, state, value);
                }
            }
        }
    }
    debug_assert!(
        phase_parents
            .iter()
            .chain(pm.iter())
            .chain(pi.iter())
            .chain(pd.iter())
            .flatten()
            .all(|parent| parent.state != State::Stop && parent.fragment.atoms[0].is_some())
    );
    let (mut i, mut j, mut state, score) = end;
    let final_state = state;
    let (query_end, oriented_target_end) = (i as u64, j as u64);
    let mut fragments = Vec::new();
    while state != State::Stop {
        let parent = match state {
            State::M => pm[idx(i, j, cols)],
            State::I => pi[idx(i, j, cols)],
            State::D => pd[idx(i, j, cols)],
            State::Stop => None,
        };
        let Some(parent) = parent else { break };
        let raw_fragment =
            coding2genome_raw_fragment(parent, state, i, j, query, bases, scoring, intron);
        let (query_advance, target_advance) =
            parent
                .fragment
                .atoms
                .iter()
                .flatten()
                .fold((0usize, 0usize), |(q, t), run| {
                    (
                        q + run.query_advance as usize * run.repeats as usize,
                        t + run.target_advance as usize * run.repeats as usize,
                    )
                });
        debug_assert!(i >= query_advance && j >= target_advance);
        i -= query_advance;
        j -= target_advance;
        fragments.push((parent.fragment, raw_fragment));
        state = parent.state;
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
        score,
        raw_trace,
        trace,
    }
}

pub fn align_coding_to_genome_suboptimal(
    query: &Sequence,
    target: &Sequence,
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        coding2genome_alignment_oriented(
            query,
            target,
            &target.bases,
            scoring,
            intron,
            Strand::Forward,
            Some(forbidden),
        )
    });
    if both_strands {
        let reverse = reverse_complement(&target.bases);
        out.extend(enumerate_suboptimal_pair(
            query,
            target,
            threshold,
            |forbidden| {
                coding2genome_alignment_oriented(
                    query,
                    target,
                    &reverse,
                    scoring,
                    intron,
                    Strand::Reverse,
                    Some(forbidden),
                )
            },
        ));

        let reverse_query = Sequence {
            id: query.id.clone(),
            bases: reverse_complement(&query.bases),
        };
        for (bases, strand) in [
            (&target.bases, Strand::Forward),
            (&reverse, Strand::Reverse),
        ] {
            for mut alignment in
                enumerate_suboptimal_pair(&reverse_query, target, threshold, |forbidden| {
                    coding2genome_alignment_oriented(
                        &reverse_query,
                        target,
                        bases,
                        scoring,
                        intron,
                        strand,
                        Some(forbidden),
                    )
                })
            {
                alignment.query_start = query.bases.len() as u64 - alignment.query_start;
                alignment.query_end = query.bases.len() as u64 - alignment.query_end;
                alignment.query_strand = Strand::Reverse;
                out.push(alignment);
            }
        }
    }
    out
}

pub fn align_coding_to_genome_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_coding_to_genome_suboptimal(
                query,
                target,
                scoring,
                intron,
                threshold,
                both_strands,
            ));
        }
    }
    out
}

/// Align every coding-DNA/genome pair, emitting one alignment per requested target strand.
pub fn align_coding_to_genome_database(
    queries: &[Sequence],
    targets: &[Sequence],
    scoring: Scoring,
    intron: IntronScoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_coding_to_genome(
                query,
                target,
                scoring,
                intron,
                both_strands,
            ));
        }
    }
    out
}

fn align_protein_to_dna_direct(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    strand: Strand,
    forbidden: Option<&HashSet<(usize, usize)>>,
) -> Alignment {
    let oriented = if strand == Strand::Forward {
        target.bases.clone()
    } else {
        reverse_complement(&target.bases)
    };
    let (n, m, cols) = (query.bases.len(), oriented.len(), oriented.len() + 1);
    let size = (n + 1) * (m + 1);
    let mut mm = vec![NEG_INF; size];
    let mut ii = vec![NEG_INF; size];
    let mut dd = vec![NEG_INF; size];
    let mut pm: Vec<Option<P2Parent>> = vec![None; size];
    let mut pi: Vec<Option<P2Parent>> = vec![None; size];
    let mut pd: Vec<Option<P2Parent>> = vec![None; size];
    let local = model == Model::Local || model == Model::Ungapped;
    let bestfit = model == Model::BestFit;
    mm[0] = 0;
    for i in 0..=n {
        for j in 0..=m {
            let k = idx(i, j, cols);
            if i == 0 && j == 0 {
                continue;
            }
            if (local && (i == 0 || j == 0)) || (bestfit && i == 0) {
                mm[k] = 0;
            }
            if i > 0 {
                let p = idx(i - 1, j, cols);
                let mut candidate = (NEG_INF, State::Stop, None);
                p2_update(
                    &mut candidate,
                    add(mm[p], scoring.codon_gap_open),
                    State::M,
                    P2Parent {
                        state: State::M,
                        op: Op::Insert,
                        transition_id: 2,
                        query_advance: 1,
                        target_advance: 0,
                    },
                );
                p2_update(
                    &mut candidate,
                    add(ii[p], scoring.codon_gap_extend),
                    State::I,
                    P2Parent {
                        state: State::I,
                        op: Op::Insert,
                        transition_id: 4,
                        query_advance: 1,
                        target_advance: 0,
                    },
                );
                if !(local && candidate.0 < 0) {
                    ii[k] = candidate.0;
                    pi[k] = candidate.2;
                }
            }
            if j >= 3 {
                let p = idx(i, j - 3, cols);
                let mut candidate = (NEG_INF, State::Stop, None);
                p2_update(
                    &mut candidate,
                    add(mm[p], scoring.codon_gap_open),
                    State::M,
                    P2Parent {
                        state: State::M,
                        op: Op::Delete,
                        transition_id: 3,
                        query_advance: 0,
                        target_advance: 3,
                    },
                );
                p2_update(
                    &mut candidate,
                    add(dd[p], scoring.codon_gap_extend),
                    State::D,
                    P2Parent {
                        state: State::D,
                        op: Op::Delete,
                        transition_id: 5,
                        query_advance: 0,
                        target_advance: 3,
                    },
                );
                if !(local && candidate.0 < 0) {
                    dd[k] = candidate.0;
                    pd[k] = candidate.2;
                }
            }
            if i > 0 && j >= 3 {
                let p = idx(i - 1, j - 3, cols);
                let aa = translate_dna(&oriented[j - 3..j], 0)[0];
                let sub = protein_score(query.bases[i - 1], aa, scoring);
                let mut candidate = (NEG_INF, State::Stop, None);
                if !forbidden_match_transition(forbidden, i - 1, j - 3, 1, 3) {
                    for (value, state) in [(mm[p], State::M), (ii[p], State::I), (dd[p], State::D)]
                    {
                        p2_update(
                            &mut candidate,
                            add(value, sub),
                            state,
                            P2Parent {
                                state,
                                op: Op::Match,
                                transition_id: 1,
                                query_advance: 1,
                                target_advance: 3,
                            },
                        );
                    }
                }
                for advance in [1_usize, 2, 4, 5] {
                    if j < advance {
                        continue;
                    }
                    let fp = idx(i, j - advance, cols);
                    for (value, state) in
                        [(mm[fp], State::M), (ii[fp], State::I), (dd[fp], State::D)]
                    {
                        p2_update(
                            &mut candidate,
                            add(value, scoring.frameshift),
                            state,
                            P2Parent {
                                state,
                                op: Op::Frameshift,
                                transition_id: if advance < 3 { 6 } else { 7 },
                                query_advance: 0,
                                target_advance: advance as u32,
                            },
                        );
                    }
                }
                if local && candidate.0 < 0 {
                    mm[k] = 0;
                    pm[k] = None;
                } else {
                    mm[k] = candidate.0;
                    pm[k] = candidate.2;
                }
            }
        }
    }
    let mut end = (0usize, 0usize, State::M, NEG_INF);
    for i in 0..=n {
        for j in 0..=m {
            let valid = if bestfit { i == n } else { true };
            if !valid {
                continue;
            }
            let k = idx(i, j, cols);
            for (value, state) in [(mm[k], State::M), (ii[k], State::I), (dd[k], State::D)] {
                if value > end.3
                    || (value == end.3 && (i, j, state.rank()) > (end.0, end.1, end.2.rank()))
                {
                    end = (i, j, state, value);
                }
            }
        }
    }
    let (mut i, mut j, mut state, score) = end;
    let final_state = state;
    let (query_end, target_end) = (i as u64, j as u64);
    let mut trace: Vec<TraceRun> = Vec::new();
    let mut reversed_raw_fragments: Vec<Vec<RawStep>> = Vec::new();
    while state != State::Stop {
        let k = idx(i, j, cols);
        let parent = match state {
            State::M => pm[k],
            State::I => pi[k],
            State::D => pd[k],
            State::Stop => None,
        };
        let Some(parent) = parent else { break };
        let mut raw_fragment = Vec::new();
        if state == State::M {
            let epsilon = match parent.state {
                State::I => Some(6),
                State::D => Some(7),
                _ => None,
            };
            if let Some(transition_id) = epsilon {
                raw_fragment.push(RawStep {
                    transition_id,
                    query_advance: 0,
                    target_advance: 0,
                    score: 0,
                });
            }
        }
        match parent.op {
            Op::Match => raw_fragment.push(RawStep {
                transition_id: 1,
                query_advance: 1,
                target_advance: 3,
                score: protein_score(
                    query.bases[i - 1],
                    translate_dna(&oriented[j - 3..j], 0)[0],
                    scoring,
                ),
            }),
            Op::Insert | Op::Delete => raw_fragment.push(RawStep {
                transition_id: parent.transition_id,
                query_advance: parent.query_advance,
                target_advance: parent.target_advance,
                score: if matches!(parent.transition_id, 2 | 3) {
                    scoring.codon_gap_open
                } else {
                    scoring.codon_gap_extend
                },
            }),
            Op::Frameshift => {
                let short = parent.target_advance % 3;
                raw_fragment.push(RawStep {
                    transition_id: if short == 1 { 8 } else { 9 },
                    query_advance: 0,
                    target_advance: short,
                    score: scoring.frameshift,
                });
                if parent.target_advance > 3 {
                    raw_fragment.push(RawStep {
                        transition_id: 11,
                        query_advance: 0,
                        target_advance: 3,
                        score: 0,
                    });
                } else {
                    raw_fragment.push(RawStep {
                        transition_id: 10,
                        query_advance: 0,
                        target_advance: 0,
                        score: 0,
                    });
                }
            }
            _ => unreachable!("protein2dna specialized traceback operation"),
        }
        reversed_raw_fragments.push(raw_fragment);
        i -= parent.query_advance as usize;
        j -= parent.target_advance as usize;
        if let Some(last) = trace.last_mut()
            && last.op == parent.op
            && last.transition_id == parent.transition_id
            && last.query_advance == parent.query_advance
            && last.target_advance == parent.target_advance
        {
            last.repeats += 1;
        } else {
            trace.push(TraceRun {
                transition_id: parent.transition_id,
                op: parent.op,
                query_advance: parent.query_advance,
                target_advance: parent.target_advance,
                repeats: 1,
            });
        }
        state = parent.state;
    }
    trace.reverse();
    reversed_raw_fragments.reverse();
    let mut raw_trace = vec![RawStep {
        transition_id: 0,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    }];
    for fragment in reversed_raw_fragments {
        raw_trace.extend(fragment);
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
    let (target_start, target_end) = if strand == Strand::Forward {
        (j as u64, target_end)
    } else {
        (
            target.bases.len() as u64 - j as u64,
            target.bases.len() as u64 - target_end,
        )
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
        score,
        raw_trace,
        trace,
    }
}

/// Exhaustive protein-to-DNA alignment over the original nucleotide sequence.
/// It implements the upstream `protein2dna` codon-affine graph, including
/// target-side 1/2/4/5-nt frameshift transitions.
pub fn align_protein_to_dna(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::with_capacity(if both_strands { 2 } else { 1 });
    out.push(align_protein_to_dna_direct(
        query,
        target,
        model,
        scoring,
        Strand::Forward,
        None,
    ));
    if both_strands {
        out.push(align_protein_to_dna_direct(
            query,
            target,
            model,
            scoring,
            Strand::Reverse,
            None,
        ));
    }
    out
}

/// Pair-disjoint suboptimal protein-to-DNA paths, enumerated independently
/// on each requested target orientation.
pub fn align_protein_to_dna_suboptimal(
    query: &Sequence,
    target: &Sequence,
    model: Model,
    scoring: Scoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = enumerate_suboptimal_pair(query, target, threshold, |forbidden| {
        align_protein_to_dna_direct(
            query,
            target,
            model,
            scoring,
            Strand::Forward,
            Some(forbidden),
        )
    });
    if both_strands {
        out.extend(enumerate_suboptimal_pair(
            query,
            target,
            threshold,
            |forbidden| {
                align_protein_to_dna_direct(
                    query,
                    target,
                    model,
                    scoring,
                    Strand::Reverse,
                    Some(forbidden),
                )
            },
        ));
    }
    out
}

pub fn align_protein_to_dna_database_suboptimal(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    threshold: Score,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_to_dna_suboptimal(
                query,
                target,
                model,
                scoring,
                threshold,
                both_strands,
            ));
        }
    }
    out
}

pub fn align_protein_to_dna_database(
    queries: &[Sequence],
    targets: &[Sequence],
    model: Model,
    scoring: Scoring,
    both_strands: bool,
) -> Vec<Alignment> {
    let mut out = Vec::new();
    for query in queries {
        for target in targets {
            out.extend(align_protein_to_dna(
                query,
                target,
                model,
                scoring,
                both_strands,
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn s(id: &str, b: &str) -> Sequence {
        Sequence {
            id: id.into(),
            bases: b.bytes().collect(),
        }
    }
    #[test]
    fn heuristic_seed_refinement_matches_exhaustive_and_reduces_area() {
        let pseudo_dna = |length: usize, mut state: u64| {
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((state >> 32) & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let query_bases = pseudo_dna(600, 17);
        let mut target_bases = pseudo_dna(900, 91);
        target_bases[350..650].copy_from_slice(&query_bases[150..450]);
        let query = Sequence {
            id: "query".into(),
            bases: query_bases,
        };
        let target = Sequence {
            id: "target".into(),
            bases: target_bases,
        };
        let config = HeuristicConfig {
            word_len: 12,
            padding: 128,
            max_word_occurrences: 32,
        };
        let region = heuristic_region(&query.bases, &target.bases, config)
            .expect("planted exact segment must seed");
        assert!(
            (region.1 - region.0) * (region.3 - region.2)
                < query.bases.len() * target.bases.len() * 3 / 4
        );
        let exhaustive = align(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        let heuristic =
            align_heuristic_pair(&query, &target, Model::Local, Scoring::default(), config);
        assert_eq!(heuristic.score, exhaustive.score);
        assert_eq!(
            (
                heuristic.query_start,
                heuristic.query_end,
                heuristic.target_start,
                heuristic.target_end
            ),
            (
                exhaustive.query_start,
                exhaustive.query_end,
                exhaustive.target_start,
                exhaustive.target_end
            )
        );
        assert_eq!(heuristic.vulgar(), exhaustive.vulgar());
    }

    #[test]
    fn est2genome_spliced_heuristic_matches_exhaustive() {
        let pseudo_dna = |length: usize, mut state: u64| {
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((state >> 32) & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let exon1 = pseudo_dna(120, 17);
        let exon2 = pseudo_dna(120, 23);
        let mut query_bases = exon1.clone();
        query_bases.extend_from_slice(&exon2);
        let mut target_bases = pseudo_dna(3_000, 91);
        let mut intron_bases = pseudo_dna(200, 101);
        intron_bases[0..2].copy_from_slice(b"GT");
        intron_bases[198..200].copy_from_slice(b"AG");
        target_bases[1_200..1_320].copy_from_slice(&exon1);
        target_bases[1_320..1_520].copy_from_slice(&intron_bases);
        target_bases[1_520..1_640].copy_from_slice(&exon2);
        let query = Sequence {
            id: "est".into(),
            bases: query_bases,
        };
        let target = Sequence {
            id: "genome".into(),
            bases: target_bases,
        };
        let intron = IntronScoring {
            min_len: 200,
            max_len: 200,
            open_penalty: -30,
            force_gtag: true,
        };
        let config = HeuristicConfig {
            word_len: 12,
            padding: 64,
            max_word_occurrences: 32,
        };
        let region =
            spliced_target_region(&query.bases, &target.bases, config, intron.max_len as usize)
                .expect("planted exons must seed");
        assert!(region.1 - region.0 < target.bases.len() / 2);
        let exhaustive = align_est2genome_stranded(&query, &target, intron, false);
        let heuristic = align_est2genome_heuristic_pair(&query, &target, intron, false, config);
        assert_eq!(heuristic.score, exhaustive.score);
        assert_eq!(heuristic.vulgar(), exhaustive.vulgar());
        assert!(heuristic.trace.iter().any(|run| run.op == Op::Intron));

        let cdna_exhaustive = align_cdna_to_genome(&query, &target, Scoring::default(), intron);
        let cdna_heuristic = align_cdna_to_genome_database_heuristic(
            std::slice::from_ref(&query),
            std::slice::from_ref(&target),
            Scoring::default(),
            intron,
            config,
        )
        .pop()
        .unwrap();
        assert_eq!(cdna_heuristic.vulgar(), cdna_exhaustive.vulgar());

        let coding_exhaustive =
            align_coding_to_genome(&query, &target, Scoring::default(), intron, false)
                .pop()
                .unwrap();
        let coding_heuristic = align_coding_to_genome_database_heuristic(
            std::slice::from_ref(&query),
            std::slice::from_ref(&target),
            Scoring::default(),
            intron,
            false,
            config,
        )
        .pop()
        .unwrap();
        assert_eq!(coding_heuristic.vulgar(), coding_exhaustive.vulgar());

        let mut query_intron = pseudo_dna(200, 202);
        query_intron[0..2].copy_from_slice(b"GT");
        query_intron[198..200].copy_from_slice(b"AG");
        let mut genome_query_bases = exon1;
        genome_query_bases.extend_from_slice(&query_intron);
        genome_query_bases.extend_from_slice(&exon2);
        let genome_query = Sequence {
            id: "query_genome".into(),
            bases: genome_query_bases,
        };
        let genome_exhaustive =
            align_genome_to_genome(&genome_query, &target, Scoring::default(), intron);
        let genome_heuristic = align_genome_to_genome_database_heuristic(
            std::slice::from_ref(&genome_query),
            std::slice::from_ref(&target),
            Scoring::default(),
            intron,
            config,
        )
        .pop()
        .unwrap();
        assert_eq!(genome_heuristic.vulgar(), genome_exhaustive.vulgar());
    }

    #[test]
    fn protein2genome_translated_heuristic_matches_exhaustive() {
        let pseudo_dna = |length: usize, mut state: u64| {
            (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((state >> 32) & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let mut amino_state = 505_u64;
        let amino_acids = (0..80)
            .map(|_| {
                amino_state = amino_state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1);
                b"ACDEFGHIKLMNPQRSTVWY"[((amino_state >> 32) % 20) as usize]
            })
            .collect::<Vec<_>>();
        let codon = |amino_acid| match amino_acid {
            b'A' => b"GCT",
            b'C' => b"TGT",
            b'D' => b"GAT",
            b'E' => b"GAA",
            b'F' => b"TTT",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'K' => b"AAA",
            b'L' => b"CTG",
            b'M' => b"ATG",
            b'N' => b"AAT",
            b'P' => b"CCT",
            b'Q' => b"CAA",
            b'R' => b"CGT",
            b'S' => b"TCT",
            b'T' => b"ACT",
            b'V' => b"GTT",
            b'W' => b"TGG",
            b'Y' => b"TAT",
            _ => unreachable!(),
        };
        let mut coding = Vec::with_capacity(3 * amino_acids.len());
        for &amino_acid in &amino_acids {
            coding.extend_from_slice(codon(amino_acid));
        }
        let mut target_bases = pseudo_dna(5_000, 303);
        let mut intron_bases = pseudo_dna(200, 404);
        intron_bases[0..2].copy_from_slice(b"GT");
        intron_bases[198..200].copy_from_slice(b"AG");
        target_bases[2_000..2_120].copy_from_slice(&coding[..120]);
        target_bases[2_120..2_320].copy_from_slice(&intron_bases);
        target_bases[2_320..2_440].copy_from_slice(&coding[120..]);
        let query = Sequence {
            id: "protein".into(),
            bases: amino_acids,
        };
        let target = Sequence {
            id: "genome".into(),
            bases: target_bases,
        };
        let intron = IntronScoring {
            min_len: 200,
            max_len: 200,
            open_penalty: -30,
            force_gtag: true,
        };
        let config = HeuristicConfig {
            word_len: 6,
            padding: 64,
            max_word_occurrences: 32,
        };
        let region = protein_to_dna_target_region(
            &query.bases,
            &target.bases,
            config,
            intron.max_len as usize,
        )
        .expect("translated exon words must seed");
        assert!(region.1 - region.0 < target.bases.len() / 2);
        let exhaustive =
            align_protein_to_genome(&query, &target, Scoring::default(), intron, false)
                .pop()
                .unwrap();
        let heuristic = align_protein_to_genome_database_heuristic(
            std::slice::from_ref(&query),
            std::slice::from_ref(&target),
            Scoring::default(),
            intron,
            false,
            config,
        )
        .pop()
        .unwrap();
        assert_eq!(heuristic.vulgar(), exhaustive.vulgar());
        assert!(heuristic.trace.iter().any(|run| run.op == Op::Intron));

        let exhaustive_bestfit =
            align_protein_to_genome_bestfit(&query, &target, Scoring::default(), intron, false)
                .pop()
                .unwrap();
        let heuristic_bestfit = align_protein_to_genome_bestfit_database_heuristic(
            std::slice::from_ref(&query),
            std::slice::from_ref(&target),
            Scoring::default(),
            intron,
            false,
            config,
        )
        .pop()
        .unwrap();
        assert_eq!(heuristic_bestfit.vulgar(), exhaustive_bestfit.vulgar());

        let reverse_target = Sequence {
            id: "reverse_genome".into(),
            bases: reverse_complement(&target.bases),
        };
        let reverse_exhaustive =
            align_protein_to_genome(&query, &reverse_target, Scoring::default(), intron, true);
        let reverse_heuristic = align_protein_to_genome_database_heuristic(
            std::slice::from_ref(&query),
            std::slice::from_ref(&reverse_target),
            Scoring::default(),
            intron,
            true,
            config,
        );
        assert_eq!(
            reverse_heuristic
                .iter()
                .map(Alignment::vulgar)
                .collect::<Vec<_>>(),
            reverse_exhaustive
                .iter()
                .map(Alignment::vulgar)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn heuristic_without_valid_words_falls_back_to_exhaustive() {
        let query = s("query", "NNNNNNNNNNNNNNNN");
        let target = s("target", "ACGTACGTACGTACGT");
        let config = HeuristicConfig::default();
        assert!(heuristic_region(&query.bases, &target.bases, config).is_none());
        let expected = align(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        let actual =
            align_heuristic_pair(&query, &target, Model::Local, Scoring::default(), config);
        assert_eq!(actual.vulgar(), expected.vulgar());
    }

    #[test]
    fn suboptimal_coding_paths_are_score_ordered_and_pair_disjoint() {
        let query = s("query", "ATGGCTATGGCT");
        let target = s("target", "ATGGCTTTTATGGCT");
        let alignments = align_coding2coding_suboptimal(&query, &target, Scoring::default(), 5);
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(
            alignments
                .windows(2)
                .all(|pair| pair[0].score >= pair[1].score)
        );
        assert_eq!(
            alignments[0].vulgar(),
            align_coding2coding(&query, &target, Scoring::default()).vulgar()
        );
        let mut forbidden = HashSet::new();
        for alignment in &alignments {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn codon_match_suboptimal_exclusion_blocks_each_nucleotide_pair() {
        let alignment = Alignment {
            query_id: "query".into(),
            target_id: "target".into(),
            query_start: 3,
            query_end: 6,
            query_strand: Strand::Forward,
            target_start: 9,
            target_end: 12,
            target_len: 20,
            target_strand: Strand::Forward,
            score: 1,
            raw_trace: Vec::new(),
            trace: vec![TraceRun {
                transition_id: 1,
                op: Op::Match,
                query_advance: 3,
                target_advance: 3,
                repeats: 1,
            }],
        };
        let mut forbidden = HashSet::new();
        assert_eq!(
            forbid_alignment_pairs(&alignment, 20, 20, &mut forbidden),
            3
        );
        assert_eq!(forbidden, HashSet::from([(3, 9), (4, 10), (5, 11)]));
        assert!(forbidden_match_transition(Some(&forbidden), 3, 9, 3, 3));
        assert!(forbidden_match_transition(Some(&forbidden), 4, 10, 3, 3));
    }

    #[test]
    fn suboptimal_local_paths_are_score_ordered_and_pair_disjoint() {
        let query = s("query", "ACGTACGTACGT");
        let target = s("target", "ACGTGGGGACGT");
        let alignments = align_suboptimal_pair(&query, &target, Scoring::default(), 5);
        assert!(alignments.len() >= 2);
        assert!(
            alignments
                .windows(2)
                .all(|pair| pair[0].score >= pair[1].score)
        );
        let optimal = align(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(alignments[0].vulgar(), optimal.vulgar());
        let mut forbidden = HashSet::new();
        for alignment in &alignments {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn model_parser_accepts_all_upstream_affine_short_names() {
        for (name, expected) in [
            ("u", Model::Ungapped),
            ("a:g", Model::Global),
            ("a:b", Model::BestFit),
            ("a:l", Model::Local),
            ("a:o", Model::Overlap),
        ] {
            assert_eq!(
                name.parse::<Model>().expect("upstream short name"),
                expected
            );
        }
    }

    #[test]
    fn protein_affine_uses_upstream_blosum62() {
        let a = align_protein(
            &s("q", "ARN"),
            &s("t", "ADN"),
            Model::Global,
            Scoring::default(),
        );
        assert_eq!(a.score, 8);
        assert_eq!(a.sugar(), "sugar: q 0 3 . t 0 3 . 8");
        assert_eq!(a.cigar(), "cigar: q 0 3 . t 0 3 . 8  M 3");
    }

    #[test]
    fn protein2genome_bestfit_spans_the_complete_query() {
        let query = s("protein", "WMM");
        let target = s("genome", "ATGATG");
        let local = align_protein_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
        )
        .pop()
        .expect("local forward alignment");
        let bestfit = align_protein_to_genome_bestfit(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
        )
        .pop()
        .expect("best-fit forward alignment");
        assert_eq!((local.query_start, local.query_end), (1, 3));
        assert_eq!((bestfit.query_start, bestfit.query_end), (0, 3));
        assert_eq!(bestfit.score, -8);
        assert_eq!(
            bestfit.vulgar(),
            "vulgar: protein 0 3 . genome 0 6 + -8 G 1 0 M 2 6"
        );
        assert!(bestfit.score < local.score);
        assert_eq!(
            bestfit
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            bestfit.score
        );
        let query_span = bestfit
            .trace
            .iter()
            .map(|run| u64::from(run.query_advance) * run.repeats)
            .sum::<u64>();
        assert_eq!(query_span, query.bases.len() as u64);
    }

    #[test]
    fn protein2genome_honors_configured_codon_gap_penalty() {
        let query = s("protein", "MAM");
        let target = s("genome", "ATGGCTNNNATG");
        let relaxed = Scoring {
            codon_gap_open: -1,
            codon_gap_extend: -1,
            ..Scoring::default()
        };
        assert_eq!(
            protein2genome_score_with_scoring(&query, &target, relaxed, IntronScoring::default()),
            13
        );
        assert_eq!(
            protein2genome_score(&query, &target, IntronScoring::default()),
            9
        );
    }

    #[test]
    fn suboptimal_protein2genome_recovers_two_phase_intron_loci() {
        let query = Sequence {
            id: "protein".into(),
            bases: b"MADQLTEQIAEFKEAFSLFDKDGDGTITT".to_vec(),
        };
        let locus = b"ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGGAGGCCTTCTCCCTCTTTGACAAGGATGGAGATGGCACTATTACCACC";
        let mut bases = locus.to_vec();
        bases.extend_from_slice(&[b'N'; 20]);
        bases.extend_from_slice(locus);
        let target = Sequence {
            id: "genome".into(),
            bases,
        };
        let alignments = align_protein_to_genome_suboptimal(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            100,
            false,
            false,
        );
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(alignments[0].score >= alignments[1].score);
        assert!(
            alignments
                .iter()
                .take(2)
                .all(|alignment| alignment.trace.iter().any(|run| run.op == Op::Intron))
        );
        let mut forbidden = HashSet::new();
        for alignment in alignments.iter().take(2) {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn protein2genome_matches_upstream_phase_intron_oracle() {
        let query = Sequence {
            id: "protein".into(),
            bases: b"MADQLTEQIAEFKEAFSLFDKDGDGTITT".to_vec(),
        };
        let target = Sequence {
            id: "genome".into(),
            bases: b"ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGGAGGCCTTCTCCCTCTTTGACAAGGATGGAGATGGCACTATTACCACC".to_vec(),
        };
        let alignment = align_protein_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
        )
        .pop()
        .expect("one forward protein2genome alignment");
        assert_eq!(alignment.score, 125);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.query_advance))
                .sum::<u64>(),
            alignment.query_end - alignment.query_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.target_advance))
                .sum::<u64>(),
            alignment.target_end - alignment.target_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .filter(|step| step.transition_id == 122)
                .count(),
            42
        );
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end,
            ),
            (0, 29, 0, 133)
        );
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48"
        );
        assert!(alignment.trace.iter().any(|run| run.op == Op::SplitCodon));
        assert!(alignment.trace.iter().any(|run| run.op == Op::Intron));

        let reverse_target = Sequence {
            id: "reverse-genome".into(),
            bases: reverse_complement(&target.bases),
        };
        let reverse = align_protein_to_genome(
            &query,
            &reverse_target,
            Scoring::default(),
            IntronScoring::default(),
            true,
        )
        .into_iter()
        .find(|candidate| candidate.target_strand == Strand::Reverse)
        .expect("reverse-strand protein2genome alignment");
        assert_eq!(reverse.score, 125);
        assert_eq!((reverse.target_start, reverse.target_end), (133, 0));
    }

    #[test]
    fn suboptimal_genome2genome_recovers_two_spliced_loci() {
        let query = s(
            "query",
            concat!(
                "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "TGAAGCAGCTCCAGCTATGTGTGAAGAA",
                "GAGGACAGCACTGCCTTGGTGTGTGACAATG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "GCTCTGGGCTCTGTAAGGCCGGCTTTGCT"
            ),
        );
        let locus = concat!(
            "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
            "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
            "TGAAGCAGCTCCAGCTATGTGTGAAGAA",
            "GAGGACAGCACTGCCTTGGTGTGTGACAATGGC",
            "TCTGGGCTCTGTAAGGCCGGCTTTGCT"
        );
        let target = s("target", &format!("{locus}{}{}", "N".repeat(20), locus));
        let alignments = align_genome_to_genome_suboptimal(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            400,
        );
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(alignments[0].score >= alignments[1].score);
        assert!(
            alignments
                .iter()
                .take(2)
                .all(|alignment| alignment.trace.iter().any(|run| run.op == Op::Intron))
        );
        let mut forbidden = HashSet::new();
        for alignment in alignments.iter().take(2) {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn genome2genome_matches_upstream_joint_intron_oracle() {
        let query = s(
            "query",
            concat!(
                "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "TGAAGCAGCTCCAGCTATGTGTGAAGAA",
                "GAGGACAGCACTGCCTTGGTGTGTGACAATG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "GCTCTGGGCTCTGTAAGGCCGGCTTTGCT"
            ),
        );
        let target = s(
            "target",
            concat!(
                "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "TGAAGCAGCTCCAGCTATGTGTGAAGAA",
                "GAGGACAGCACTGCCTTGGTGTGTGACAATGGC",
                "TCTGGGCTCTGTAAGGCCGGCTTTGCT"
            ),
        );
        let alignment = align_genome_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        assert_eq!(alignment.score, 557);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.query_advance))
                .sum::<u64>(),
            alignment.query_end - alignment.query_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.target_advance))
                .sum::<u64>(),
            alignment.target_end - alignment.target_start
        );
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.target_start, 0);
        assert!(alignment.trace.iter().any(|run| run.op == Op::Intron));
        assert!(alignment.trace.iter().any(|run| {
            run.op == Op::Intron && run.query_advance == 0 && run.target_advance > 0
        }));
        assert!(alignment.trace.iter().any(|run| {
            run.op == Op::Intron && run.query_advance > 0 && run.target_advance == 0
        }));
        let (query_span, target_span) =
            alignment.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
                (
                    q + u64::from(run.query_advance) * run.repeats,
                    t + u64::from(run.target_advance) * run.repeats,
                )
            });
        assert_eq!(query_span, alignment.query_end - alignment.query_start);
        assert_eq!(target_span, alignment.target_end - alignment.target_start);
    }

    #[test]
    fn genome2genome_rolling_score_rows_preserve_full_traceback() {
        let query = s(
            "query",
            concat!(
                "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "TGAAGCAGCTCCAGCTATGTGTGAAGAA",
                "GAGGACAGCACTGCCTTGGTGTGTGACAATG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "GCTCTGGGCTCTGTAAGGCCGGCTTTGCT"
            ),
        );
        let target = s(
            "target",
            concat!(
                "AGCCCAGCCAAGCACTGTCAGGAATCCTG",
                "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG",
                "TGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGC",
                "TCTGGGCTCTGTAAGGCCGGCTTTGCT"
            ),
        );
        let full = align_genome_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        let rolling = align_genome_to_genome_forbidden(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            None,
            true,
            None,
            None,
            None,
            None,
            1,
        );
        assert_eq!(rolling.score, full.score);
        assert_eq!(rolling.sugar(), full.sugar());
        assert_eq!(rolling.cigar(), full.cigar());
        assert_eq!(rolling.vulgar(), full.vulgar());
        assert_eq!(rolling.raw_trace, full.raw_trace);
    }

    #[test]
    fn checkpointed_genome2genome_preserves_the_complete_traceback() {
        let query = s("query", &"ACGT".repeat(45));
        let target = s(
            "target",
            &format!("{}{}", "N".repeat(23), "ACGT".repeat(45)),
        );
        let full = align_genome_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        let checkpointed = align_genome_to_genome_checkpointed(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            None,
        );
        assert_eq!(checkpointed.score, full.score);
        assert_eq!(checkpointed.sugar(), full.sugar());
        assert_eq!(checkpointed.cigar(), full.cigar());
        assert_eq!(checkpointed.vulgar(), full.vulgar());
        assert_eq!(checkpointed.trace, full.trace);
        assert_eq!(checkpointed.raw_trace, full.raw_trace);
    }

    #[test]
    fn budgeted_genome2genome_checkpoint_replay_preserves_the_complete_traceback() {
        let query = s("query", &"ACGT".repeat(45));
        let target = s(
            "target",
            &format!("{}{}", "N".repeat(23), "ACGT".repeat(45)),
        );
        let full = align_genome_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        let checkpointed = align_genome_to_genome_checkpointed(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            Some(0),
        );
        assert_eq!(checkpointed.score, full.score);
        assert_eq!(checkpointed.sugar(), full.sugar());
        assert_eq!(checkpointed.cigar(), full.cigar());
        assert_eq!(checkpointed.vulgar(), full.vulgar());
        assert_eq!(checkpointed.trace, full.trace);
        assert_eq!(checkpointed.raw_trace, full.raw_trace);
    }

    #[test]
    fn genome2genome_parent_block_preserves_the_forward_endpoint() {
        let query = s("query", "ACGTACGTACGTACGTACGTACGT");
        let target = s("target", "NNNNACGTACGTACGTACGTACGTACGTNN");
        let full = align_genome_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        // A forward pass with only row zero of each parent lane retained is
        // the first half of checkpoint replay.  It cannot traceback yet, but
        // must select exactly the same local endpoint and score.
        let blocked = align_genome_to_genome_forbidden(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            None,
            true,
            Some((0, 0)),
            None,
            None,
            None,
            1,
        );
        assert_eq!(blocked.score, full.score);
        assert_eq!(blocked.query_end, full.query_end);
        assert_eq!(blocked.target_end, full.target_end);
        assert!(blocked.trace.is_empty());

        let replayed_block = align_genome_to_genome_forbidden(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            None,
            true,
            Some((12, 24)),
            Some((
                full.query_end as usize,
                full.target_end as usize,
                GenomeState::UtrM,
            )),
            None,
            None,
            1,
        );
        assert_eq!(replayed_block.score, full.score);
        assert_eq!(replayed_block.query_end, full.query_end);
        assert_eq!(replayed_block.target_end, full.target_end);
        // A retained row may contain a five-base frameshift whose source lies
        // just before the block boundary; the next replay begins there.
        assert!((7..=12).contains(&replayed_block.query_start));
        assert!(!replayed_block.trace.is_empty());
    }

    #[test]
    fn genome2genome_checkpoints_capture_score_history_and_live_queues() {
        let query = s("query", &"TGG".repeat(40));
        let target = s(
            "target",
            &format!(
                "{}GT{}AG{}",
                "TGG".repeat(12),
                "N".repeat(30),
                "TGG".repeat(28)
            ),
        );
        let intron = IntronScoring {
            min_len: 30,
            max_len: 80,
            ..IntronScoring::default()
        };
        let mut checkpoints = Vec::new();
        let rolling = align_genome_to_genome_forbidden(
            &query,
            &target,
            Scoring::default(),
            intron,
            None,
            true,
            None,
            None,
            None,
            Some(&mut checkpoints),
            17,
        );
        assert!(rolling.score > 0);
        assert_eq!(
            checkpoints.first().map(|checkpoint| checkpoint.row),
            Some(0)
        );
        assert_eq!(
            checkpoints.last().map(|checkpoint| checkpoint.row),
            Some(query.bases.len())
        );
        assert!(checkpoints.iter().all(|checkpoint| {
            checkpoint.history.last().map(|row| row.row) == Some(checkpoint.row)
                && checkpoint.history.len() <= genome_checkpoint_history_rows(intron)
        }));
        assert!(
            checkpoints
                .iter()
                .any(|checkpoint| checkpoint.queues.checkpoint_bytes() > 0)
        );

        let checkpoint = checkpoints.last().expect("final checkpoint");
        let cols = target.bases.len() + 1;
        let history_rows = genome_checkpoint_history_rows(intron);
        let mut u5m = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut u5i = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut u5d = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut cm = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut ci = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut cd = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut u3m = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut u3i = GenomeScoreMatrix::rolling(cols, history_rows);
        let mut u3d = GenomeScoreMatrix::rolling(cols, history_rows);
        restore_genome_checkpoint_scores(
            checkpoint, &mut u5m, &mut u5i, &mut u5d, &mut cm, &mut ci, &mut cd, &mut u3m,
            &mut u3i, &mut u3d,
        );
        let final_row = checkpoint.history.last().expect("checkpoint final row");
        assert_eq!(u5m.row_clone(final_row.row), final_row.u5m);
        assert_eq!(cm.row_clone(final_row.row), final_row.cm);
        assert_eq!(u3d.row_clone(final_row.row), final_row.u3d);
    }

    #[test]
    fn genome_parent_block_retains_only_its_reconstruction_rows() {
        let cols = 5;
        let mut parents = GenomeParentMatrix::block(cols, 3, 4);
        parents[idx(3, 2, cols)] = Some(GenomeParent {
            state: GenomeState::UtrM,
            fragment: CompactGenomeFragment::from_trace_fragment(TraceFragment::empty()),
        });
        parents[idx(2, 2, cols)] = Some(GenomeParent {
            state: GenomeState::CdsM,
            fragment: CompactGenomeFragment::from_trace_fragment(TraceFragment::empty()),
        });
        assert_eq!(
            parents[idx(3, 2, cols)].map(|parent| parent.state),
            Some(GenomeState::UtrM)
        );
    }

    #[test]
    fn genome_parent_fragment_is_compacter_than_public_trace_fragment() {
        assert!(
            size_of::<CompactGenomeFragment>() < size_of::<TraceFragment>(),
            "the resident g2g parent fragment must stay smaller than its public trace form"
        );
    }

    #[test]
    fn genome_rolling_scores_restore_checkpoint_history_by_row_number() {
        let mut scores = GenomeScoreMatrix::rolling(3, 4);
        scores.restore_row(8, &[10, 11, 12]);
        scores.restore_row(9, &[20, 21, 22]);
        assert_eq!(scores.row_clone(8), vec![10, 11, 12]);
        assert_eq!(scores.row_clone(9), vec![20, 21, 22]);
        scores.start_row(10, NEG_INF);
        assert_eq!(scores.row_clone(10), vec![NEG_INF; 3]);
    }

    #[test]
    fn genome2genome_cds_query_and_joint_introns_cover_all_phases() {
        let with_intron = |codon: &str, phase: usize| {
            format!(
                "{}{}GT{}AG{}{}",
                codon.repeat(10),
                &codon[..phase],
                "N".repeat(30),
                &codon[phase..],
                codon.repeat(9)
            )
        };
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        for phase in 0..3 {
            let query_only = align_genome_to_genome(
                &s("query", &with_intron("CGT", phase)),
                &s("target", &"AGA".repeat(20)),
                Scoring::default(),
                intron,
            );
            assert!(
                query_only.trace.iter().any(|run| {
                    run.op == Op::Intron && run.query_advance > 0 && run.target_advance == 0
                }),
                "missing query-only phase {phase}: {}",
                query_only.vulgar()
            );
            assert_eq!(
                query_only.trace.iter().any(|run| run.op == Op::SplitCodon),
                phase != 0,
                "query-only phase {phase}: {}",
                query_only.vulgar()
            );

            let target_only = align_genome_to_genome(
                &s("query", &"CGT".repeat(20)),
                &s("target", &with_intron("AGA", phase)),
                Scoring::default(),
                intron,
            );
            assert!(
                target_only.trace.iter().any(|run| {
                    run.op == Op::Intron && run.query_advance == 0 && run.target_advance > 0
                }),
                "missing target-only phase {phase}: {}",
                target_only.vulgar()
            );
            assert_eq!(
                target_only.trace.iter().any(|run| run.op == Op::SplitCodon),
                phase != 0,
                "target-only phase {phase}: {}",
                target_only.vulgar()
            );

            let joint = align_genome_to_genome(
                &s("query", &with_intron("CGT", phase)),
                &s("target", &with_intron("AGA", phase)),
                Scoring::default(),
                intron,
            );
            assert!(
                joint.trace.iter().any(|run| {
                    run.op == Op::Intron && run.query_advance == 0 && run.target_advance > 0
                }) && joint.trace.iter().any(|run| {
                    run.op == Op::Intron && run.query_advance > 0 && run.target_advance == 0
                }),
                "missing joint phase {phase}: {}",
                joint.vulgar()
            );
            assert_eq!(
                joint.trace.iter().any(|run| run.op == Op::SplitCodon),
                phase != 0,
                "joint phase {phase}: {}",
                joint.vulgar()
            );
            let (query_span, target_span) =
                joint.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
                    (
                        q + u64::from(run.query_advance) * run.repeats,
                        t + u64::from(run.target_advance) * run.repeats,
                    )
                });
            assert_eq!(query_span, joint.query_end - joint.query_start);
            assert_eq!(target_span, joint.target_end - joint.target_start);
        }
    }

    #[test]
    fn suboptimal_cdna2genome_recovers_two_spliced_loci() {
        let query = s("cdna", "ACGTACGT");
        let locus = format!("ACGTGT{}AGACGT", "N".repeat(30));
        let target = s("genome", &format!("{locus}{}{}", "N".repeat(20), locus));
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        let alignments =
            align_cdna_to_genome_suboptimal(&query, &target, Scoring::default(), intron, 5);
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(alignments[0].score >= alignments[1].score);
        assert!(
            alignments
                .iter()
                .take(2)
                .all(|alignment| alignment.trace.iter().any(|run| run.op == Op::Intron))
        );
        let mut forbidden = HashSet::new();
        for alignment in alignments.iter().take(2) {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn cdna2genome_matches_upstream_composed_oracle() {
        let source = include_str!("../../../upstream/src/model/cdna2genome.test.c");
        let extract = |name: &str| {
            source
                .split(name)
                .nth(1)
                .expect("fixture chain")
                .split("NULL,")
                .nth(1)
                .expect("fixture bases")
                .split("0, Sequence_Strand_UNKNOWN")
                .next()
                .expect("fixture terminator")
                .lines()
                .filter_map(|line| {
                    let text = line.trim();
                    text.split_once('"')
                        .and_then(|(_, text)| text.split_once('"'))
                        .map(|(bases, _)| bases)
                })
                .collect::<String>()
        };
        let query = s("cdna", &extract("*qy = Sequence_create"));
        let target = s("genome", &extract("*tg = Sequence_create"));
        assert_eq!(query.bases.len(), 270);
        assert_eq!(target.bases.len(), 540);
        let dp = cdna2genome_dp(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            None,
        );
        assert_eq!(
            dp.value(dp.end.2, idx(dp.end.0, dp.end.1, dp.cols)),
            dp.end.3
        );
        let alignment = align_cdna_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
        );
        assert_eq!(alignment.score, 1281);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.query_advance))
                .sum::<u64>(),
            alignment.query_end - alignment.query_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.target_advance))
                .sum::<u64>(),
            alignment.target_end - alignment.target_start
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 321)
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 400)
        );
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end
            ),
            (0, 270, 54, 486)
        );
        assert!(
            alignment
                .trace
                .iter()
                .filter(|run| run.op == Op::Intron)
                .count()
                >= 3
        );
        assert_eq!(
            cdna2genome_score(
                &query,
                &target,
                Scoring::default(),
                IntronScoring::default()
            ),
            1281
        );
    }

    #[test]
    fn suboptimal_coding2genome_recovers_two_phase_intron_loci() {
        let query = s(
            "coding",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let locus = b"AGCCCAGCCAAGCACTGTCAGGAATCCTGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT";
        let mut bases = locus.to_vec();
        bases.extend_from_slice(&[b'N'; 20]);
        bases.extend_from_slice(locus);
        let target = Sequence {
            id: "genome".into(),
            bases,
        };
        let alignments = align_coding_to_genome_suboptimal(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            150,
            false,
        );
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(alignments[0].score >= alignments[1].score);
        assert!(
            alignments
                .iter()
                .take(2)
                .all(|alignment| alignment.trace.iter().any(|run| run.op == Op::Intron))
        );
        let mut forbidden = HashSet::new();
        for alignment in alignments.iter().take(2) {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn coding2genome_matches_upstream_phase_intron_oracle() {
        let query = s(
            "coding",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let target = s(
            "genome",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let alignment = align_coding_to_genome(
            &query,
            &target,
            Scoring::default(),
            IntronScoring::default(),
            false,
        )
        .pop()
        .expect("one forward coding2genome alignment");
        assert_eq!(alignment.score, 194);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.query_advance))
                .sum::<u64>(),
            alignment.query_end - alignment.query_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.target_advance))
                .sum::<u64>(),
            alignment.target_end - alignment.target_start
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .filter(|step| step.transition_id == 122)
                .count(),
            30
        );
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end,
            ),
            (0, 117, 0, 151)
        );
        assert_eq!(
            alignment.vulgar(),
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87"
        );
        assert!(alignment.trace.iter().any(|run| run.op == Op::SplitCodon));
        assert!(alignment.trace.iter().any(|run| run.op == Op::Intron));

        let reverse_target = Sequence {
            id: "reverse-genome".into(),
            bases: reverse_complement(&target.bases),
        };
        let reverse = align_coding_to_genome(
            &query,
            &reverse_target,
            Scoring::default(),
            IntronScoring::default(),
            true,
        )
        .into_iter()
        .find(|candidate| candidate.target_strand == Strand::Reverse)
        .expect("reverse-strand coding2genome alignment");
        assert_eq!(reverse.score, 194);
        assert_eq!((reverse.target_start, reverse.target_end), (151, 0));
    }

    #[test]
    fn suboptimal_protein2dna_paths_are_pair_disjoint() {
        let query = s("protein", "MA");
        let target = s("target", "ATGGCTNNNNNNATGGCT");
        let alignments = align_protein_to_dna_suboptimal(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            5,
            false,
        );
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(
            alignments
                .windows(2)
                .all(|pair| pair[0].score >= pair[1].score)
        );
        assert_eq!(
            alignments[0].vulgar(),
            align_protein_to_dna(&query, &target, Model::Local, Scoring::default(), false,)[0]
                .vulgar()
        );
        let mut forbidden = HashSet::new();
        for alignment in &alignments {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn protein_to_dna_matches_upstream_model_test_score() {
        let alignment = align_protein_to_dna(
            &s("protein", "NNNNNNMADQLTEQIAEFKEAFSLFDKDGTVHNCXWYFSGRWDGTITT"),
            &s("dna", "ATGGCTGACCAGCTGACTGAGGAGCAGATTGCAGAGTTCNAAGGAGGCCTTCTCCCTCTTTGACAAGGATGGANNACTGTCCATAATTGCTGGTACTTCAGCGGTCGATGGGATGGCACTCTGACCACC"),
            Model::Local,
            Scoring::default(),
            false,
        )
        .pop()
        .expect("one forward alignment");
        // `upstream/src/model/protein2dna.test.c` asserts this same score.
        assert_eq!(alignment.score, 134);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end
            ),
            (6, 48, 0, 129)
        );
    }
    #[test]
    fn protein_to_dna_direct_viterbi_recovers_frameshift_path() {
        let scoring = Scoring {
            frameshift: -1,
            ..Scoring::default()
        };
        let alignment = align_protein_to_dna(
            &s("protein", "MA"),
            &s("dna", "ATGTGCT"),
            Model::Local,
            scoring,
            false,
        )
        .pop()
        .expect("one forward alignment");
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end,
                alignment.score
            ),
            (0, 2, 0, 7, 8)
        );
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 2 . dna 0 7 + 8 M 1 3 F 0 1 M 1 3"
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 8)
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 10)
        );
    }
    #[test]
    fn protein_to_dna_reverse_frame_uses_forward_reference_coordinates() {
        let alignments = align_protein_to_dna(
            &s("protein", "MA"),
            &s("dna", "TTTAGCCATGGG"),
            Model::Local,
            Scoring::default(),
            true,
        );
        let alignment = alignments
            .iter()
            .find(|a| a.score == 9 && a.target_strand == Strand::Reverse)
            .expect("reverse-complement frame must recover MA");
        assert_eq!(
            (
                alignment.target_start,
                alignment.target_end,
                alignment.target_strand
            ),
            (9, 3, Strand::Reverse)
        );
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 2 . dna 9 3 - 9 M 2 6"
        );
    }
    #[test]
    fn protein_to_dna_projects_trace_and_coordinates() {
        let alignments = align_protein_to_dna(
            &s("protein", "MA"),
            &s("dna", "ATGGCT"),
            Model::Local,
            Scoring::default(),
            false,
        );
        let alignment = alignments
            .iter()
            .find(|a| a.score == 9)
            .expect("frame one must recover MA");
        assert_eq!(
            (
                alignment.query_start,
                alignment.query_end,
                alignment.target_start,
                alignment.target_end,
                alignment.query_strand,
                alignment.target_strand,
            ),
            (0, 2, 0, 6, Strand::Unknown, Strand::Forward)
        );
        assert_eq!(alignment.cigar(), "cigar: protein 0 2 . dna 0 6 + 9  M 6");
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 2 . dna 0 6 + 9 M 2 6"
        );
    }

    #[test]
    fn split_codon_uses_protein2dna_translation_and_blosum() {
        assert_eq!(split_codon_score(b'M', b"AT", b"G"), Some(5));
        assert_eq!(split_codon_score(b'M', b"A", b"TG"), Some(5));
        assert_eq!(split_codon_score(b'M', b"AT", b""), None);
    }
    #[test]
    fn splice_predictor_enforces_directional_canonical_sites() {
        let forward_donor = b"AAAGTAAAA";
        assert!(splice_score(forward_donor, 3, SpliceType::DonorForward, true).is_some());
        assert!(splice_score(forward_donor, 3, SpliceType::AcceptorForward, true).is_none());
        let reverse_donor = b"AAAAACAAA";
        assert!(splice_score(reverse_donor, 4, SpliceType::DonorReverse, true).is_some());
        assert!(splice_score(reverse_donor, 4, SpliceType::DonorForward, true).is_none());
        // Values reported by upstream/src/sequence/splice.test.c after C rounding.
        assert_eq!(
            splice_score(
                b"TGGCGCTCCTGGCCTCTGCCCGTAAGCACTTGGTGGGACTGGG",
                21,
                SpliceType::DonorForward,
                false
            ),
            Some(6)
        );
        assert_eq!(
            splice_score(
                b"CAGCCCTGCTCTTTCCTCAGGAGCTTCAGAGGCCGAGGATGCC",
                18,
                SpliceType::AcceptorForward,
                false
            ),
            Some(12)
        );
    }
    #[test]
    fn standard_genetic_code_and_frames_are_correct() {
        assert_eq!(translate_dna(b"ATGGCTTAA", 0), b"MA*");
        assert_eq!(translate_dna(b"AATGGCT", 1), b"MA");
        assert_eq!(translate_dna(b"ATN", 0), b"X");
    }
    #[test]
    fn nucleic_matrix_matches_upstream_iupac_scores() {
        let scoring = Scoring::default();
        assert_eq!(dna_score(b'R', b'A', scoring), 1);
        assert_eq!(dna_score(b'N', b'A', scoring), -2);
        assert_eq!(dna_score(b'A', b'A', scoring), 5);
    }

    #[test]
    fn protein_to_dna_ir_matches_upstream_frameshift_graph() {
        let ir = model::protein_to_dna(model::Scope::Anywhere);
        assert!(ir.validate().is_ok());
        assert!(ir.transitions.iter().any(|transition| {
            transition.label == model::Label::Frameshift
                && transition.query_advance == 0
                && transition.target_advance == 1
        }));
        assert!(ir.transitions.iter().any(|transition| {
            transition.label == model::Label::Frameshift
                && transition.query_advance == 0
                && transition.target_advance == 3
        }));
    }

    #[test]
    fn intron_fragment_preserves_splice_loop_splice_order() {
        let atom = |transition_id, op, target_advance| TraceRun {
            transition_id,
            op,
            query_advance: 0,
            target_advance,
            repeats: 1,
        };
        let mut trace = Vec::new();
        TraceFragment::intron(
            atom(10, Op::Splice5, 2),
            atom(11, Op::Intron, 36),
            atom(12, Op::Splice3, 2),
        )
        .append_to(&mut trace);
        assert_eq!(
            trace.iter().map(|run| run.op).collect::<Vec<_>>(),
            vec![Op::Splice5, Op::Intron, Op::Splice3]
        );
    }

    #[test]
    fn gff3_preserves_reverse_target_strand_coordinates() {
        let alignment = align_protein_to_dna(
            &s("protein", "MA"),
            &s("dna", "TTTAGCCATGGG"),
            Model::Local,
            Scoring::default(),
            true,
        )
        .into_iter()
        .find(|a| a.target_strand == Strand::Reverse && a.score == 9)
        .expect("reverse hit");
        assert_eq!(
            alignment.gff3(),
            "dna\texonerate-rs\tmatch\t4\t9\t9\t-\t.\tID=match_protein_dna;Target=protein 1 2\ndna\texonerate-rs\tmatch_part\t4\t9\t.\t-\t.\tID=match_protein_dna.part1;Parent=match_protein_dna;Target=protein 1 2"
        );
    }

    #[test]
    fn gff3_renders_match_and_parts_from_trace() {
        let alignment =
            align_coding2coding(&s("q", "ATGGCT"), &s("t", "ATGGCT"), Scoring::default());
        assert_eq!(
            alignment.gff3(),
            "t\texonerate-rs\tmatch\t2\t4\t11\t+\t.\tID=match_q_t;Target=q 2 4\nt\texonerate-rs\tmatch_part\t2\t4\t.\t+\t.\tID=match_q_t.part1;Parent=match_q_t;Target=q 2 4"
        );
    }

    #[test]
    fn vulgar_reports_split_codon_as_s() {
        let alignment = Alignment {
            query_id: "protein".into(),
            target_id: "genome".into(),
            query_start: 0,
            query_end: 1,
            query_strand: Strand::Unknown,
            target_start: 0,
            target_end: 3,
            target_len: 3,
            target_strand: Strand::Forward,
            score: 5,
            raw_trace: Vec::new(),
            trace: vec![TraceRun {
                transition_id: 1,
                op: Op::SplitCodon,
                query_advance: 1,
                target_advance: 3,
                repeats: 1,
            }],
        };
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 1 . genome 0 3 + 5 S 1 3"
        );
    }

    #[test]
    fn phase_intron_fragment_preserves_split_codon_order() {
        let atom = |id, op, q, t| TraceRun {
            transition_id: id,
            op,
            query_advance: q,
            target_advance: t,
            repeats: 1,
        };
        let mut trace = Vec::new();
        TraceFragment::phase_intron(
            atom(20, Op::SplitCodon, 0, 2),
            atom(21, Op::Splice5, 0, 2),
            atom(22, Op::Intron, 0, 42),
            atom(23, Op::Splice3, 0, 2),
            atom(24, Op::SplitCodon, 1, 1),
        )
        .append_to(&mut trace);
        assert_eq!(
            trace.iter().map(|run| run.op).collect::<Vec<_>>(),
            vec![
                Op::SplitCodon,
                Op::Splice5,
                Op::Intron,
                Op::Splice3,
                Op::SplitCodon
            ]
        );
    }

    #[test]
    fn intron_candidate_window_uses_state_rank_for_score_ties() {
        let mut window = IntronCandidateWindow::default();
        window.insert(IntronCandidate {
            start: 2,
            score: 9,
            state_rank: 1,
        });
        window.insert(IntronCandidate {
            start: 3,
            score: 9,
            state_rank: 3,
        });
        assert_eq!(
            window
                .best()
                .map(|candidate| (candidate.start, candidate.state_rank)),
            Some((3, 3))
        );
    }

    #[test]
    fn intron_window_snapshots_preserve_candidates_independently() {
        let mut window = IntronCandidateWindow::default();
        window.insert(IntronCandidate {
            start: 2,
            score: 9,
            state_rank: 1,
        });
        let snapshot = window.clone();
        window.expire_before(3);
        assert!(window.best().is_none());
        assert_eq!(snapshot.best().map(|candidate| candidate.start), Some(2));

        let mut joint = JointIntronWindow::default();
        joint.insert(JointIntronCandidate {
            query_start: 4,
            target_start: 7,
            score: 12,
            state_rank: 2,
        });
        let snapshot = joint.clone();
        joint.expire_before(5);
        assert!(joint.best().is_none());
        assert_eq!(
            snapshot
                .best()
                .map(|candidate| (candidate.query_start, candidate.target_start)),
            Some((4, 7))
        );
    }

    #[test]
    fn genome_intron_queue_checkpoint_size_includes_live_candidates() {
        let mut queues = GenomeIntronQueues::new(1);
        let empty = 5 * size_of::<IntronCandidateWindow>() + 5 * size_of::<JointIntronWindow>();
        assert_eq!(queues.checkpoint_bytes(), empty);

        queues.query[0].insert(IntronCandidate {
            start: 3,
            score: 11,
            state_rank: 2,
        });
        queues.joint_coding[2][0].insert(JointIntronCandidate {
            query_start: 4,
            target_start: 6,
            score: 12,
            state_rank: 3,
        });
        assert_eq!(
            queues.checkpoint_bytes(),
            empty + size_of::<IntronCandidate>() + size_of::<JointIntronCandidate>()
        );
    }

    #[test]
    fn genome_intron_queue_snapshot_preserves_each_queue_family() {
        let mut queues = GenomeIntronQueues::new(1);
        queues.query[0].insert(IntronCandidate {
            start: 1,
            score: 10,
            state_rank: 1,
        });
        queues.query_coding[1][0].insert(IntronCandidate {
            start: 2,
            score: 11,
            state_rank: 2,
        });
        queues.u3_query[0].insert(IntronCandidate {
            start: 3,
            score: 12,
            state_rank: 3,
        });
        queues.joint[0].insert(JointIntronCandidate {
            query_start: 3,
            target_start: 4,
            score: 12,
            state_rank: 3,
        });
        queues.u3_joint[0].insert(JointIntronCandidate {
            query_start: 4,
            target_start: 5,
            score: 13,
            state_rank: 4,
        });
        queues.joint_coding[2][0].insert(JointIntronCandidate {
            query_start: 5,
            target_start: 6,
            score: 13,
            state_rank: 4,
        });

        let snapshot = queues.clone();
        queues.query[0].expire_before(2);
        queues.query_coding[1][0].expire_before(3);
        queues.u3_query[0].expire_before(4);
        queues.joint[0].expire_before(4);
        queues.u3_joint[0].expire_before(5);
        queues.joint_coding[2][0].expire_before(6);

        assert!(queues.query[0].best().is_none());
        assert!(queues.query_coding[1][0].best().is_none());
        assert!(queues.u3_query[0].best().is_none());
        assert!(queues.joint[0].best().is_none());
        assert!(queues.u3_joint[0].best().is_none());
        assert!(queues.joint_coding[2][0].best().is_none());
        assert_eq!(
            snapshot.query[0].best().map(|candidate| candidate.start),
            Some(1)
        );
        assert_eq!(
            snapshot.query_coding[1][0]
                .best()
                .map(|candidate| candidate.start),
            Some(2)
        );
        assert_eq!(
            snapshot.u3_query[0].best().map(|candidate| candidate.start),
            Some(3)
        );
        assert_eq!(
            snapshot.joint[0]
                .best()
                .map(|candidate| (candidate.query_start, candidate.target_start)),
            Some((3, 4))
        );
        assert_eq!(
            snapshot.u3_joint[0]
                .best()
                .map(|candidate| (candidate.query_start, candidate.target_start)),
            Some((4, 5))
        );
        assert_eq!(
            snapshot.joint_coding[2][0]
                .best()
                .map(|candidate| (candidate.query_start, candidate.target_start)),
            Some((5, 6))
        );
    }

    #[test]
    fn genome_checkpoint_history_covers_phase_introns_and_frameshifts() {
        assert_eq!(
            genome_checkpoint_history_rows(IntronScoring {
                min_len: 4,
                ..IntronScoring::default()
            }),
            8
        );
        assert_eq!(
            genome_checkpoint_history_rows(IntronScoring {
                min_len: 30,
                ..IntronScoring::default()
            }),
            34
        );
    }

    #[test]
    fn genome_checkpoint_queue_bound_covers_a_populated_snapshot() {
        let intron = IntronScoring {
            min_len: 2,
            max_len: 4,
            ..IntronScoring::default()
        };
        let mut queues = GenomeIntronQueues::new(2);
        for column in 0..2 {
            for start in 0..3 {
                queues.query[column].entries.push_back(IntronCandidate {
                    start,
                    score: start as Score,
                    state_rank: 0,
                });
                queues.query_coding[0][column]
                    .entries
                    .push_back(IntronCandidate {
                        start,
                        score: start as Score,
                        state_rank: 0,
                    });
                queues.joint[column]
                    .entries
                    .push_back(JointIntronCandidate {
                        query_start: start,
                        target_start: start,
                        score: start as Score,
                        state_rank: 0,
                    });
                queues.joint_coding[0][column]
                    .entries
                    .push_back(JointIntronCandidate {
                        query_start: start,
                        target_start: start,
                        score: start as Score,
                        state_rank: 0,
                    });
            }
        }
        assert!(
            genome_checkpoint_queue_upper_bound_bytes(2, intron)
                >= queues.checkpoint_bytes() as u128
        );
    }

    #[test]
    fn intron_window_keeps_best_eligible_start_in_linear_space() {
        let mut window = IntronWindow::default();
        window.insert(2, 10);
        window.insert(5, 7);
        window.insert(8, 12);
        assert_eq!(window.best(), Some((8, 12)));
        window.expire_before(9);
        assert_eq!(window.best(), None);
        window.insert(10, 4);
        window.insert(11, 6);
        assert_eq!(window.best(), Some((11, 6)));
    }

    fn c_concat(source: &str, marker: &str) -> String {
        let section = source.split_once(marker).expect("C test marker").1;
        let section = section.split_once(';').expect("C declaration terminator").0;
        let mut uncommented = String::new();
        let mut rest = section;
        while let Some((before, after_open)) = rest.split_once("/*") {
            uncommented.push_str(before);
            rest = after_open.split_once("*/").expect("closed C comment").1;
        }
        uncommented.push_str(rest);
        let mut out = String::new();
        let mut quoted = false;
        for chunk in uncommented.split('"') {
            if quoted {
                out.push_str(chunk);
            }
            quoted = !quoted;
        }
        out
    }
    #[test]
    fn suboptimal_est2genome_recovers_two_spliced_loci() {
        let source = include_str!("../../../upstream/src/model/est2genome.test.c");
        let query = Sequence {
            id: "query".into(),
            bases: c_concat(source, "*query_seq =").into_bytes(),
        };
        let locus = c_concat(source, "*target_seq =").into_bytes();
        let mut duplicated = locus.clone();
        duplicated.extend_from_slice(&[b'N'; 20]);
        duplicated.extend_from_slice(&locus);
        let target = Sequence {
            id: "target".into(),
            bases: duplicated,
        };
        let alignments =
            align_est2genome_suboptimal(&query, &target, IntronScoring::default(), 100, false);
        assert!(alignments.len() >= 2, "{alignments:?}");
        assert!(alignments[0].score >= 100);
        assert!(alignments[1].score >= 100);
        assert!(alignments[0].score >= alignments[1].score);
        assert!(
            alignments
                .iter()
                .take(2)
                .all(|alignment| alignment.trace.iter().any(|run| run.op == Op::Intron))
        );
        let mut forbidden = HashSet::new();
        for alignment in alignments.iter().take(2) {
            assert!(
                forbid_alignment_pairs(
                    alignment,
                    query.bases.len(),
                    target.bases.len(),
                    &mut forbidden,
                ) > 0
            );
        }
    }

    #[test]
    fn est2genome_gff3_splits_match_parts_at_intron() {
        let source = include_str!("../../../upstream/src/model/est2genome.test.c");
        let query = Sequence {
            id: "query".into(),
            bases: c_concat(source, "*query_seq =").into_bytes(),
        };
        let target = Sequence {
            id: "target".into(),
            bases: c_concat(source, "*target_seq =").into_bytes(),
        };
        let gff = align_est2genome(&query, &target, IntronScoring::default()).gff3();
        assert_eq!(
            gff.lines()
                .filter(|line| line.contains("\tmatch_part\t"))
                .count(),
            4
        );
        assert!(!gff.contains("\t11\t202\t"));
    }

    #[test]
    fn est2genome_traceback_runs_on_upstream_test_fixture() {
        let source = include_str!("../../../upstream/src/model/est2genome.test.c");
        let query = Sequence {
            id: "query".into(),
            bases: c_concat(source, "*query_seq =").into_bytes(),
        };
        let target = Sequence {
            id: "target".into(),
            bases: c_concat(source, "*target_seq =").into_bytes(),
        };
        let alignment = align_est2genome(&query, &target, IntronScoring::default());
        assert_eq!(alignment.score, 157);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 6)
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| step.transition_id == 8)
        );
        assert!(alignment.trace.iter().any(|run| run.op == Op::Splice5));
        assert!(alignment.trace.iter().any(|run| run.op == Op::Intron));
        assert!(alignment.trace.iter().any(|run| run.op == Op::Splice3));
        assert_eq!(alignment.query_end - alignment.query_start, 44);
        assert_eq!(alignment.target_end - alignment.target_start, 238);
    }

    #[test]
    fn est2genome_reports_reverse_target_strand() {
        let source = include_str!("../../../upstream/src/model/est2genome.test.c");
        let query = Sequence {
            id: "query".into(),
            bases: c_concat(source, "*query_seq =").into_bytes(),
        };
        let original = c_concat(source, "*target_seq =").into_bytes();
        let target = Sequence {
            id: "target".into(),
            bases: reverse_complement(&original),
        };
        let alignment = align_est2genome(&query, &target, IntronScoring::default());
        assert_eq!(alignment.score, 157);
        assert_eq!(alignment.target_strand, Strand::Reverse);
        assert_eq!((alignment.target_start, alignment.target_end), (238, 0));
    }

    #[test]
    fn est2genome_score_runs_on_upstream_test_fixture() {
        let source = include_str!("../../../upstream/src/model/est2genome.test.c");
        let query = Sequence {
            id: "query".into(),
            bases: c_concat(source, "*query_seq =").into_bytes(),
        };
        let target = Sequence {
            id: "target".into(),
            bases: c_concat(source, "*target_seq =").into_bytes(),
        };
        assert_eq!(
            est2genome_score(&query, &target, IntronScoring::default()),
            157
        );
    }

    #[test]
    fn coding2coding_trace_conserves_reported_coordinates() {
        let query = s("q", "ATGGCTATG");
        let target = s("t", "ATGGCTTATG");
        let alignment = align_coding2coding(
            &query,
            &target,
            Scoring {
                frameshift: -1,
                ..Scoring::default()
            },
        );
        assert_eq!(alignment.score, 16);
        let query_advance: u64 = alignment
            .trace
            .iter()
            .map(|run| u64::from(run.query_advance) * run.repeats)
            .sum();
        let target_advance: u64 = alignment
            .trace
            .iter()
            .map(|run| u64::from(run.target_advance) * run.repeats)
            .sum();
        assert_eq!(query_advance, alignment.query_end - alignment.query_start);
        assert_eq!(
            target_advance,
            alignment.target_end - alignment.target_start
        );
        assert!(alignment.trace.iter().any(|run| run.op == Op::Frameshift));
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.query_advance))
                .sum::<u64>(),
            query_advance
        );
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| u64::from(step.target_advance))
                .sum::<u64>(),
            target_advance
        );
        assert!(
            alignment
                .raw_trace
                .iter()
                .any(|step| matches!(step.transition_id, 8 | 9 | 12 | 13))
        );
    }

    #[test]
    fn coding2coding_matches_upstream_score_oracle() {
        let query = s(
            "qy",
            "AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        let target = s(
            "tg",
            "AGCCCAGCCAAACACTGTCAGGAATCCTGTNNNGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTGCCTTGGTGTGTGACAATGGCNNTCTGGGCTCTGTAAGGCCGGCTTTGCT",
        );
        assert_eq!(
            coding2coding_score(&query, &target, Scoring::default()),
            169
        );
    }

    #[test]
    fn coding_to_coding_ir_has_query_and_target_frameshifts() {
        let ir = model::coding_to_coding();
        assert!(ir.validate().is_ok());
        assert!(
            ir.transitions
                .iter()
                .any(|edge| edge.label == model::Label::Frameshift && edge.query_advance == 1)
        );
        assert!(
            ir.transitions
                .iter()
                .any(|edge| edge.label == model::Label::Frameshift && edge.target_advance == 1)
        );
    }
    #[test]
    fn target_phase_intron_ir_has_all_three_phase_paths() {
        let ir = model::target_phase_intron(30, 200_000);
        assert!(ir.validate().is_ok());
        assert_eq!(
            ir.transitions
                .iter()
                .filter(|edge| edge.label == model::Label::SplitCodon)
                .count(),
            4
        );
        let post_advances: Vec<_> = ir
            .transitions
            .iter()
            .filter(|edge| edge.label == model::Label::SplitCodon && edge.query_advance == 1)
            .map(|edge| edge.target_advance)
            .collect();
        assert_eq!(post_advances, vec![2, 1]);
    }
    #[test]
    fn est2genome_ir_has_bounded_stereo_introns() {
        let ir = model::est2genome(30, 200_000);
        assert!(ir.validate().is_ok());
        assert_eq!(
            ir.transitions
                .iter()
                .filter(|edge| edge.label == model::Label::Intron)
                .count(),
            2
        );
        assert!(ir.transitions.iter().any(|edge| matches!(
            edge.kernel,
            model::ScoreKernel::IntronClose {
                min_len: 30,
                max_len: 200_000
            }
        )));
    }
    #[test]
    fn generic_c4_ungapped_translated_matches_upstream_22_point_oracle() {
        let query = s("dna1", "CGATCAGCTAGCTAGCTACGATCGATCGAT");
        let target = s("dna2", "CGATACGATCGCTCTGAGATCTCGACTCAG");
        let alignment = align_model_ir(
            &query,
            &target,
            &model::ungapped_translated(model::Scope::Anywhere),
            Scoring::default(),
            Strand::Forward,
        )
        .expect("translated ungapped IR must execute");
        assert_eq!(alignment.score, 22);
        assert!(alignment.trace.iter().all(|run| {
            run.op == Op::Match && run.query_advance == 3 && run.target_advance == 3
        }));
    }

    #[test]
    fn generic_c4_ner_matches_upstream_208_point_oracle() {
        let query = s(
            "qy",
            concat!(
                "TTTTATCTTCCCAAGAGNCCCCATNNNGCGA",
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAA",
                "GTGATTGAAATGTGGATGAAACATTTC"
            ),
        );
        let target = s(
            "tg",
            concat!(
                "TTTTATCTTCCCAAGAGCCCCATGAGGCGA",
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                "GTGANTGAAATGTGGATGAACATTTC"
            ),
        );
        let alignment = align_model_ir(
            &query,
            &target,
            &model::ner(10, 50_000, -20),
            Scoring::default(),
            Strand::Forward,
        )
        .expect("NER span must compile");
        assert_eq!(alignment.score, 208);
        assert_eq!(
            alignment
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            alignment.score
        );
        let ner = alignment
            .trace
            .iter()
            .find(|run| run.op == Op::Ner)
            .expect("upstream fixture must use one NER span");
        assert!((10..=50_000).contains(&ner.query_advance));
        assert!((10..=50_000).contains(&ner.target_advance));
        let (query_span, target_span) =
            alignment.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
                (
                    q + u64::from(run.query_advance) * run.repeats,
                    t + u64::from(run.target_advance) * run.repeats,
                )
            });
        assert_eq!(query_span, alignment.query_end - alignment.query_start);
        assert_eq!(target_span, alignment.target_end - alignment.target_start);
        assert!(alignment.vulgar().contains(" N "));
    }

    #[test]
    fn generic_c4_executor_compiles_unequal_joint_intron_long_states() {
        let exon = "ACGT".repeat(10);
        let with_intron =
            |loop_len: usize| format!("{}GT{}AG{}", &exon[..20], "N".repeat(loop_len), &exon[20..]);
        let query = s("query", &with_intron(30));
        let target = s("target", &with_intron(40));
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        let expected = align_genome_to_genome(&query, &target, Scoring::default(), intron);
        let actual = align_model_ir_with_intron(
            &query,
            &target,
            &model::joint_intron(30, 200_000),
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .expect("joint intron long state must compile");
        assert_eq!(actual.score, expected.score);
        assert_eq!(
            actual
                .raw_trace
                .iter()
                .map(|step| step.score)
                .sum::<Score>(),
            actual.score
        );
        assert_eq!((actual.query_start, actual.query_end), (0, 74));
        assert_eq!((actual.target_start, actual.target_end), (0, 84));
        assert!(actual.trace.iter().any(|run| {
            run.op == Op::Intron && run.query_advance == 30 && run.target_advance == 40
        }));
        let (query_span, target_span) = actual.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
            (
                q + u64::from(run.query_advance) * run.repeats,
                t + u64::from(run.target_advance) * run.repeats,
            )
        });
        assert_eq!(query_span, 74);
        assert_eq!(target_span, 84);
    }

    #[test]
    fn generic_c4_executor_compiles_query_intron_long_states() {
        let exon = "ACGT".repeat(10);
        let intron_sequence = format!("{}GT{}AG{}", &exon[..20], "N".repeat(30), &exon[20..]);
        let query = s("query", &intron_sequence);
        let target = s("target", &exon);
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        let alignment = align_model_ir_with_intron(
            &query,
            &target,
            &model::query_intron(30, 200_000),
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .expect("query intron long state must compile");
        assert_eq!((alignment.query_start, alignment.query_end), (0, 74));
        assert_eq!((alignment.target_start, alignment.target_end), (0, 40));
        assert!(alignment.trace.iter().any(|run| {
            run.op == Op::Intron && run.query_advance == 30 && run.target_advance == 0
        }));
        let (query_span, target_span) =
            alignment.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
                (
                    q + u64::from(run.query_advance) * run.repeats,
                    t + u64::from(run.target_advance) * run.repeats,
                )
            });
        assert_eq!(query_span, 74);
        assert_eq!(target_span, 40);
    }

    #[test]
    fn generic_c4_executor_matches_specialized_affine_scopes() {
        let query = s("query", "ACGTTGCA");
        let target = s("target", "TTACGTCGCAAA");
        for model in [Model::Global, Model::BestFit, Model::Local, Model::Overlap] {
            let expected = align(&query, &target, model, Scoring::default(), Strand::Forward);
            let actual = align_model_ir(
                &query,
                &target,
                &model.ir(),
                Scoring::default(),
                Strand::Forward,
            )
            .expect("finite affine IR must execute");
            assert_eq!(actual.score, expected.score, "{model:?}");
            assert_eq!(
                (
                    actual.query_start,
                    actual.query_end,
                    actual.target_start,
                    actual.target_end
                ),
                (
                    expected.query_start,
                    expected.query_end,
                    expected.target_start,
                    expected.target_end
                ),
                "{model:?}"
            );
            assert_eq!(actual.cigar(), expected.cigar(), "{model:?}");
        }
    }

    #[test]
    fn generic_c4_executor_exhaustively_matches_short_dna_paths() {
        let mut sequences = Vec::new();
        for length in 1..=3 {
            let count = 4usize.pow(length);
            for mut code in 0..count {
                let mut bases = vec![b'A'; length as usize];
                for base in &mut bases {
                    *base = b"ACGT"[code % 4];
                    code /= 4;
                }
                sequences.push(Sequence {
                    id: "s".into(),
                    bases,
                });
            }
        }
        for model in [Model::Global, Model::BestFit, Model::Local, Model::Overlap] {
            for query in &sequences {
                for target in &sequences {
                    let expected = align(query, target, model, Scoring::default(), Strand::Forward);
                    let actual = align_model_ir(
                        query,
                        target,
                        &model.ir(),
                        Scoring::default(),
                        Strand::Forward,
                    )
                    .expect("finite DNA IR must execute");
                    assert_eq!(
                        (
                            actual.score,
                            actual.query_start,
                            actual.query_end,
                            actual.target_start,
                            actual.target_end
                        ),
                        (
                            expected.score,
                            expected.query_start,
                            expected.query_end,
                            expected.target_start,
                            expected.target_end
                        ),
                        "model={model:?} query={:?} target={:?} actual={} expected={}",
                        query.bases,
                        target.bases,
                        actual.vulgar(),
                        expected.vulgar()
                    );
                    assert_eq!(actual.cigar(), expected.cigar());
                }
            }
        }
    }

    #[test]
    fn generic_c4_executor_matches_ungapped_scores_for_all_short_sequences() {
        let alphabet = b"ACGT";
        for query_code in 0..64usize {
            for target_code in 0..64usize {
                let decode = |mut code: usize| {
                    let mut bases = vec![b'A'; 3];
                    for base in &mut bases {
                        *base = alphabet[code % 4];
                        code /= 4;
                    }
                    Sequence {
                        id: "s".into(),
                        bases,
                    }
                };
                let query = decode(query_code);
                let target = decode(target_code);
                let expected = align(
                    &query,
                    &target,
                    Model::Ungapped,
                    Scoring::default(),
                    Strand::Forward,
                );
                let actual = align_model_ir(
                    &query,
                    &target,
                    &Model::Ungapped.ir(),
                    Scoring::default(),
                    Strand::Forward,
                )
                .expect("ungapped IR must execute");
                assert_eq!(actual.score, expected.score);
            }
        }
    }

    #[test]
    fn generic_c4_executor_matches_protein2dna_oracle() {
        let query = s("protein", "MAMAPRTEINSEQWENCE");
        let target = s(
            "dna",
            "ATGGCTATGGCTCCTCGTACTGAAATTAATTCAGAATCTGAATGGGAAAATGCT",
        );
        let expected =
            align_protein_to_dna(&query, &target, Model::Local, Scoring::default(), false)
                .pop()
                .expect("forward alignment");
        let actual = align_model_ir(
            &query,
            &target,
            &model::protein_to_dna(model::Scope::Anywhere),
            Scoring::default(),
            Strand::Forward,
        )
        .expect("protein2dna IR must execute");
        assert_eq!(actual.score, expected.score);
        assert_eq!(actual.vulgar(), expected.vulgar());
    }

    #[test]
    fn generic_c4_executor_reconstructs_all_target_split_codon_phases() {
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        for phase in 0..3 {
            let codon = "AGA";
            let target = s(
                "target",
                &format!(
                    "{codon}{}GT{}AG{}",
                    &codon[..phase],
                    "N".repeat(30),
                    &codon[phase..]
                ),
            );
            let query = s("protein", "RR");
            let expected =
                align_protein_to_genome(&query, &target, Scoring::default(), intron, false)
                    .pop()
                    .expect("forward protein2genome alignment");
            let actual = align_model_ir_with_intron(
                &query,
                &target,
                &model::protein_phase_intron(30, 200_000),
                Scoring::default(),
                intron,
                Strand::Forward,
            )
            .expect("target phase long state must compile");
            assert_eq!(actual.score, expected.score, "phase {phase}");
            assert_eq!(
                (
                    actual.query_start,
                    actual.query_end,
                    actual.target_start,
                    actual.target_end
                ),
                (
                    expected.query_start,
                    expected.query_end,
                    expected.target_start,
                    expected.target_end
                ),
                "phase {phase}: actual={} expected={}",
                actual.vulgar(),
                expected.vulgar()
            );
            assert_eq!(
                actual.trace.iter().any(|run| run.op == Op::SplitCodon),
                phase != 0,
                "phase {phase}: {}",
                actual.vulgar()
            );
            assert!(actual.trace.iter().any(|run| run.op == Op::Intron));
        }
    }

    #[test]
    fn generic_c4_executor_compiles_target_intron_long_states() {
        let exon = "ACGT".repeat(10);
        let query = s("query", &exon);
        let target = s(
            "target",
            &format!("{}GT{}AG{}", &exon[..20], "N".repeat(30), &exon[20..]),
        );
        let intron = IntronScoring {
            open_penalty: -1,
            force_gtag: true,
            ..IntronScoring::default()
        };
        let expected = align_est2genome(&query, &target, intron);
        let actual = align_model_ir_with_intron(
            &query,
            &target,
            &model::est2genome(30, 200_000),
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .expect("target intron long state must compile");
        assert_eq!(actual.score, expected.score);
        assert_eq!(
            (
                actual.query_start,
                actual.query_end,
                actual.target_start,
                actual.target_end
            ),
            (
                expected.query_start,
                expected.query_end,
                expected.target_start,
                expected.target_end
            )
        );
        assert_eq!(actual.vulgar(), expected.vulgar());
        assert!(actual.trace.iter().any(|run| run.op == Op::Intron));
    }

    #[test]
    fn generic_c4_short_intron_minimum_keeps_splice_loop_nonnegative() {
        let exon = "ACGT".repeat(10);
        let query = s("query", &exon);
        let target = s("target", &format!("{}GTAG{}", &exon[..20], &exon[20..]));
        let intron = IntronScoring {
            min_len: 2,
            max_len: 20,
            open_penalty: -1,
            force_gtag: true,
        };
        let actual = align_model_ir_with_intron(
            &query,
            &target,
            &model::est2genome(2, 20),
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .expect("short intron IR must compile without underflow");
        assert!(actual.trace.iter().any(|run| {
            run.op == Op::Intron && run.query_advance == 0 && run.target_advance == 0
        }));
    }

    #[test]
    fn model_ir_rejects_invalid_graphs() {
        use crate::model::{ModelError, ScoreKernel};
        let mut cyclic = Model::Ungapped.ir();
        cyclic.transitions[1].query_advance = 0;
        cyclic.transitions[1].target_advance = 0;
        assert_eq!(cyclic.validate(), Err(ModelError::EpsilonCycle));
        let mut intron = Model::Local.ir();
        intron.transitions[0].kernel = ScoreKernel::IntronClose {
            min_len: 10,
            max_len: 2,
        };
        assert_eq!(intron.validate(), Err(ModelError::InvalidIntronBounds));
    }
    #[test]
    fn global_has_affine_trace() {
        let a = align(
            &s("q", "ACGT"),
            &s("t", "ACT"),
            Model::Global,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(a.score, 3);
        assert_eq!(a.cigar(), "cigar: q 0 4 + t 0 3 + 3  M 2 I 1 M 1");
        assert_eq!(
            a.trace
                .iter()
                .map(|run| run.transition_id)
                .collect::<Vec<_>>(),
            vec![1, 2, 1]
        );
    }

    #[test]
    fn frameshift_trace_uses_upstream_vulgar_label() {
        let alignment = Alignment {
            query_id: "protein".into(),
            target_id: "dna".into(),
            query_start: 0,
            query_end: 1,
            query_strand: Strand::Unknown,
            target_start: 0,
            target_end: 4,
            target_len: 4,
            target_strand: Strand::Forward,
            score: -19,
            raw_trace: Vec::new(),
            trace: vec![
                TraceRun {
                    transition_id: 6,
                    op: Op::Frameshift,
                    query_advance: 0,
                    target_advance: 1,
                    repeats: 1,
                },
                TraceRun {
                    transition_id: 1,
                    op: Op::Match,
                    query_advance: 1,
                    target_advance: 3,
                    repeats: 1,
                },
            ],
        };
        assert_eq!(
            alignment.cigar(),
            "cigar: protein 0 1 . dna 0 4 + -19 D 1 M 3"
        );
        assert_eq!(
            alignment.vulgar(),
            "vulgar: protein 0 1 . dna 0 4 + -19 F 0 1 M 1 3"
        );
    }
    #[test]
    fn report_layer_merges_gap_open_and_extend() {
        let a = align(
            &s("q", "ACGGGGT"),
            &s("t", "ACGT"),
            Model::Global,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(
            a.trace
                .iter()
                .map(|run| run.transition_id)
                .collect::<Vec<_>>(),
            vec![1, 2, 4, 1]
        );
        assert_eq!(a.cigar(), "cigar: q 0 7 + t 0 4 + 0  M 2 I 3 M 2");
        assert_eq!(a.vulgar(), "vulgar: q 0 7 + t 0 4 + 0 M 2 2 G 3 0 M 2 2");
    }
    #[test]
    fn local_clips() {
        let a = align(
            &s("q", "TTACGTAA"),
            &s("t", "GGACGTCC"),
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(
            (
                a.query_start,
                a.query_end,
                a.target_start,
                a.target_end,
                a.score
            ),
            (2, 6, 2, 6, 20)
        );
    }
    #[test]
    fn local_can_start_on_a_sequence_boundary() {
        let a = align(
            &s("q", "ACG"),
            &s("t", "CGT"),
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        assert_eq!(
            (
                a.query_start,
                a.query_end,
                a.target_start,
                a.target_end,
                a.score
            ),
            (1, 3, 0, 2, 10)
        );
    }
    #[test]
    fn reverse_coordinates() {
        let a = align(
            &s("q", "ACG"),
            &s("t", "CGT"),
            Model::Local,
            Scoring::default(),
            Strand::Reverse,
        );
        assert_eq!(
            (a.target_start, a.target_end, a.target_strand),
            (3, 0, Strand::Reverse)
        );
    }
    #[test]
    fn scopes_and_ungapped_have_expected_boundaries() {
        let scoring = Scoring::default();
        let q = s("q", "ACGT");
        let target = s("t", "TTACGTT");
        let bestfit = align(&q, &target, Model::BestFit, scoring, Strand::Forward);
        assert_eq!(
            (
                bestfit.query_start,
                bestfit.query_end,
                bestfit.target_start,
                bestfit.target_end,
                bestfit.score
            ),
            (0, 4, 2, 6, 20)
        );
        let overlap = align(
            &s("q", "ACGT"),
            &s("t", "GTAA"),
            Model::Overlap,
            scoring,
            Strand::Forward,
        );
        assert_eq!(
            (
                overlap.query_start,
                overlap.query_end,
                overlap.target_start,
                overlap.target_end,
                overlap.score
            ),
            (2, 4, 0, 2, 10)
        );
        let ungapped = align(
            &s("q", "TTACG"),
            &s("t", "GGACG"),
            Model::Ungapped,
            scoring,
            Strand::Forward,
        );
        assert_eq!(
            (
                ungapped.query_start,
                ungapped.query_end,
                ungapped.target_start,
                ungapped.target_end,
                ungapped.score
            ),
            (2, 5, 2, 5, 15)
        );
    }
    #[test]
    fn pretty_alignment_renders_trace_rows() {
        let query = s("query", "ACGT");
        let target = s("target", "TTACGTGG");
        let alignment = align(
            &query,
            &target,
            Model::Local,
            Scoring::default(),
            Strand::Forward,
        );
        let pretty = alignment.pretty(&query, &target);
        assert!(pretty.contains("C4 Alignment:"));
        assert!(pretty.contains("Raw score: 20"));
        assert!(pretty.contains("ACGT\n||||\nACGT"));
    }

    #[test]
    fn reports_and_reverse_query_are_canonical() {
        let hits = align_database(
            &[s("q", "ACG")],
            &[s("t", "CGT")],
            Model::Local,
            Scoring::default(),
            true,
        );
        let reverse = hits
            .iter()
            .find(|a| a.query_strand == Strand::Reverse)
            .unwrap();
        assert_eq!(
            (
                reverse.query_start,
                reverse.query_end,
                reverse.target_start,
                reverse.target_end
            ),
            (3, 0, 0, 3)
        );
        assert_eq!(reverse.vulgar(), "vulgar: q 3 0 - t 0 3 + 15 M 3 3");
        for a in hits {
            let (q, t) = a.trace.iter().fold((0_u64, 0_u64), |(q, t), run| {
                let (dq, dt) = (run.query_advance as u64, run.target_advance as u64);
                (q + dq * run.repeats, t + dt * run.repeats)
            });
            assert_eq!(q, a.query_start.abs_diff(a.query_end));
            assert_eq!(t, a.target_start.abs_diff(a.target_end));
        }
    }
    #[test]
    fn generic_c4_query_and_joint_phase_shadows_score_split_codons() {
        let query = s("query", "ATGGGTAAAGCT");
        let compact_target = s("target", "ATGGCT");
        let intron_target = s("target", "ATGGGTAAAGCT");
        let intron = IntronScoring {
            min_len: 6,
            max_len: 6,
            open_penalty: 100,
            force_gtag: true,
        };
        let donor = splice_score(&query.bases, 4, SpliceType::DonorForward, true).unwrap();
        let acceptor = splice_score(&query.bases, 8, SpliceType::AcceptorForward, true).unwrap();

        let mut query_ir = model::query_codon_phase_intron(6, 6);
        query_ir.scope = model::Scope::Query;
        let query_alignment = align_model_ir_with_intron(
            &query,
            &compact_target,
            &query_ir,
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .unwrap();
        assert_eq!(query_alignment.score, 9 + 100 + donor + acceptor);
        assert_eq!(
            (query_alignment.query_start, query_alignment.query_end),
            (0, 12)
        );
        assert!(
            query_alignment
                .trace
                .iter()
                .any(|run| run.op == Op::SplitCodon)
        );
        assert!(query_alignment.trace.iter().any(|run| run.op == Op::Intron));

        let mut joint_ir = model::joint_codon_phase_intron(6, 6);
        joint_ir.scope = model::Scope::Query;
        let joint_alignment = align_model_ir_with_intron(
            &query,
            &intron_target,
            &joint_ir,
            Scoring::default(),
            intron,
            Strand::Forward,
        )
        .unwrap();
        assert_eq!(joint_alignment.score, 9 + 100 + 2 * (donor + acceptor));
        assert_eq!(
            (joint_alignment.query_start, joint_alignment.query_end),
            (0, 12)
        );
        assert!(
            joint_alignment
                .trace
                .iter()
                .any(|run| run.op == Op::SplitCodon)
        );
        assert!(joint_alignment.trace.iter().any(|run| run.op == Op::Intron));
    }

    #[test]
    fn translated_self_score_uses_complete_codons() {
        let sequence = s("coding", "ATGATGT");
        assert_eq!(translated_self_score(&sequence), 10);
        assert_eq!(
            translated_self_score(&sequence),
            protein_self_score(&s("aa", "MM"))
        );
    }
}
