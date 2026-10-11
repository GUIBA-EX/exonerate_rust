//! Default HSP construction and sparse gapped extension for spliced models.
//!
//! Defaults and extension/tie rules follow comparison/hspset.c. In particular,
//! protein seeds include the BLOSUM word neighbourhood, not just exact words.
//! The default upstream engine is SDP (hub/gam.c), not BSDP. Following sdp/sdp.c,
//! a reverse pass starts at HSP score centres, skips span loops and records the
//! start/span boundary. Forward extension starts on that boundary and may thaw
//! frozen intron candidates only there. Both passes retain a path maximum for
//! the 50-point gapped dropoff; reversed shadow-dependent costs are zero.
//! Rolling score rows and active-cell bitsets avoid a full score matrix. Trace
//! nodes are stored only for accepted updates; dpmemory is not an RSS cap here.
use super::*;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Hsp {
    pub query_start: usize,
    pub target_start: usize,
    pub length: usize,
    pub score: Score,
    pub centre: usize,
}

pub(super) fn defaults(config: HeuristicConfig) -> bool {
    let default = HeuristicConfig::default();
    config.word_len == default.word_len
        && config.padding == default.padding
        && config.max_word_occurrences == default.max_word_occurrences
}

fn extend(query: &[u8], target: &[u8], qi: usize, ti: usize, length: usize, protein: bool) -> Hsp {
    let score_at = |i, j| {
        if protein {
            protein_score(query[i], target[j], Scoring::default())
        } else {
            dna_score(query[i], target[j], Scoring::default())
        }
    };
    let mut start = (qi, ti);
    let mut end = (qi + length, ti + length);
    // HSP_trim_ends removes unfavourable ends of neighbourhood seeds.
    while start.0 < end.0 && score_at(start.0, start.1) <= 0 {
        start.0 += 1;
        start.1 += 1;
    }
    while start.0 < end.0 && score_at(end.0 - 1, end.1 - 1) <= 0 {
        end.0 -= 1;
        end.1 -= 1;
    }
    let mut best: Score = (0..end.0 - start.0)
        .map(|k| score_at(start.0 + k, start.1 + k))
        .sum();
    let dropoff = if protein { 20 } else { 30 };
    let (mut i, mut j) = start;
    let mut score = best;
    while i > 0 && j > 0 {
        i -= 1;
        j -= 1;
        score += score_at(i, j);
        if score >= best {
            best = score;
            start = (i, j);
        } else if score < 0 || best - score >= dropoff {
            break;
        }
    }
    (i, j) = end;
    score = best;
    while i < query.len() && j < target.len() {
        score += score_at(i, j);
        i += 1;
        j += 1;
        if score >= best {
            best = score;
            end = (i, j);
        } else if score < 0 || best - score >= dropoff {
            break;
        }
    }
    let mut sum = 0;
    let mut centre = 0;
    for k in 0..end.0 - start.0 {
        sum += score_at(start.0 + k, start.1 + k);
        centre = k;
        if sum >= best >> 1 {
            break;
        }
    }
    Hsp {
        query_start: start.0,
        target_start: start.1,
        length: end.0 - start.0,
        score: best,
        centre,
    }
}

fn protein_words(word: &[u8], position: usize, budget: Score, key: u64, keys: &mut Vec<u64>) {
    const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    if position == word.len() {
        keys.push(key);
        return;
    }
    let base = word[position];
    let identity = protein_substitution_score(base, base);
    for &candidate in ALPHABET {
        let loss = identity - protein_substitution_score(base, candidate);
        if loss <= budget {
            protein_words(
                word,
                position + 1,
                budget - loss,
                (key << 8) | u64::from(candidate),
                keys,
            );
        }
    }
}

fn find(query: &[u8], target: &[u8], protein: bool, first_only: bool) -> Vec<Hsp> {
    let k = if protein { 6 } else { 12 };
    if query.len() < k || target.len() < k {
        return Vec::new();
    }
    let key_for = |word: &[u8]| {
        if protein {
            Some(word.iter().fold(0u64, |key, &base| {
                (key << 8) | u64::from(base.to_ascii_uppercase())
            }))
        } else {
            encoded_kmer(word)
        }
    };
    // Positions form linked lists in one allocation; the index contains no
    // separately allocated word or position vector per target word.
    let mut index: HashMap<u64, usize> = HashMap::with_capacity(target.len() - k + 1);
    let mut next = vec![usize::MAX; target.len() - k + 1];
    for (j, word) in target.windows(k).enumerate() {
        if let Some(key) = key_for(word) {
            next[j] = index.insert(key, j).unwrap_or(usize::MAX);
        }
    }
    let mut horizon: HashMap<isize, usize> = HashMap::new();
    let mut out = Vec::new();
    let mut words = Vec::new();
    let mut positions = Vec::new();
    for (i, word) in query.windows(k).enumerate() {
        words.clear();
        positions.clear();
        if protein {
            protein_words(word, 0, 4, 0, &mut words);
        } else if let Some(key) = key_for(word) {
            words.push(key);
        } else {
            continue;
        }
        for key in &words {
            let mut j = index.get(key).copied().unwrap_or(usize::MAX);
            while j != usize::MAX {
                positions.push(j);
                j = next[j];
            }
        }
        positions.sort_unstable();
        positions.dedup();
        for &j in &positions {
            let diagonal = j as isize - i as isize;
            if i < *horizon.get(&diagonal).unwrap_or(&0) {
                continue;
            }
            let hsp = extend(query, target, i, j, k, protein);
            horizon.insert(diagonal, hsp.query_start + hsp.length);
            if hsp.score >= if protein { 30 } else { 75 } {
                out.push(hsp);
                if first_only {
                    return out;
                }
            }
        }
    }
    out.sort_by_key(|h| {
        (
            h.target_start + h.centre,
            h.query_start + h.centre,
            h.length,
            h.score,
        )
    });
    out
}

/// None requests the fallback for unsupported query alphabets.
pub(super) fn has_protein_hsp(query: &[u8], target: &[u8]) -> Option<bool> {
    if query
        .iter()
        .any(|b| !b"ACDEFGHIKLMNPQRSTVWY".contains(&b.to_ascii_uppercase()))
    {
        return None;
    }
    for frame in 0..3 {
        let translated = translate_dna(target, frame);
        if !find(query, &translated, true, true).is_empty() {
            return Some(true);
        }
    }
    Some(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dna_hsps(q: &[u8], t: &[u8]) -> Vec<Hsp> {
        find(q, t, false, false)
    }

    #[test]
    fn protein_neighbourhood_recovers_upstream_nonidentical_words() {
        // Exonerate 2.4.0, ungapped protein/protein, score=1, exhaustive=no.
        // No query six-mer occurs literally in the target.
        let q = b"LIVLIVLIVLIV";
        let t = b"LVVLIVLVVLIV";
        assert!(
            q.windows(6)
                .all(|word| !t.windows(6).any(|other| word == other))
        );
        let hsps = find(q, t, true, false);
        let mut observed: Vec<_> = hsps
            .iter()
            .map(|h| (h.query_start, h.target_start, h.length, h.score))
            .collect();
        observed.sort_unstable();
        assert_eq!(observed, vec![(0, 0, 12, 46), (0, 3, 9, 35), (3, 0, 9, 34)]);
        assert_eq!(hsps.iter().find(|h| h.length == 12).unwrap().centre, 5);
    }

    #[test]
    fn positive_alignment_is_not_sufficient_for_an_upstream_seed() {
        assert!(find(&[b'I'; 40], &[b'V'; 40], true, false).is_empty());
        assert!(dna_hsps(&[b'A'; 14], &[b'A'; 14]).is_empty());
        assert_eq!(dna_hsps(&[b'A'; 15], &[b'A'; 15])[0].score, 75);
        assert_eq!(
            dna_hsps(b"acgtacgtacgtacgt", b"ACGTACGTACGTACGT")[0].score,
            80
        );
    }

    #[test]
    fn unsupported_protein_neighbourhoods_request_a_fallback() {
        assert!(has_protein_hsp(b"AAAAAXAAAAA", b"GCTGCTGCTGCTGCTGCTGCT").is_none());
    }
}

// SDP's reverse pass ignores span loops. It records the cells from which a
// seeded alignment can begin; only these cells may thaw intron spans later.
#[derive(Clone, Copy)]
struct Extension {
    score: Score,
    maximum: Score,
    path: usize,
    entry: usize,
}
impl Default for Extension {
    fn default() -> Self {
        Self {
            score: NEG_INF,
            maximum: NEG_INF,
            path: usize::MAX,
            entry: usize::MAX,
        }
    }
}
type Cell = [Extension; 8];
#[derive(Clone, Copy)]
struct PathNode {
    previous: usize,
    step: RawStep,
}
#[derive(Clone, Copy)]
struct Span {
    entry: usize,
    value: Extension,
}

// Every EST transition advances target by at most two bases. Three rolling
// rows and bitsets replace an allocation and tree lookup for every DP cell.
struct Frontier<const ROWS: usize = 3> {
    cells: [Vec<Cell>; ROWS],
    active: [Vec<u64>; ROWS],
}
impl<const ROWS: usize> Frontier<ROWS> {
    fn new(query_len: usize) -> Self {
        Self {
            cells: std::array::from_fn(|_| vec![Cell::default(); query_len + 1]),
            active: std::array::from_fn(|_| vec![0; (query_len + 1).div_ceil(64)]),
        }
    }
    fn at(&mut self, j: usize, i: usize) -> &mut Cell {
        self.active[j % ROWS][i / 64] |= 1 << (i % 64);
        &mut self.cells[j % ROWS][i]
    }
    fn pop(&mut self, j: usize, cursor: &mut usize, reverse: bool) -> Option<(usize, Cell)> {
        let words = &mut self.active[j % ROWS];
        loop {
            if *cursor >= words.len() {
                return None;
            }
            let bits = words[*cursor];
            if bits != 0 {
                let bit = if reverse {
                    63 - bits.leading_zeros() as usize
                } else {
                    bits.trailing_zeros() as usize
                };
                words[*cursor] &= !(1 << bit);
                let i = *cursor * 64 + bit;
                return Some((i, std::mem::take(&mut self.cells[j % ROWS][i])));
            }
            if reverse {
                *cursor = cursor.checked_sub(1).unwrap_or(usize::MAX);
            } else {
                *cursor += 1;
            }
        }
    }
}

fn update_reverse<const ROWS: usize>(
    cells: &mut Frontier<ROWS>,
    position: (usize, usize),
    state: usize,
    source: Extension,
    delta: Score,
) {
    if source.score == NEG_INF {
        return;
    }
    let score = source.score + delta;
    if source.maximum - score > 50 {
        return;
    }
    let cell = cells.at(position.0, position.1);
    if score > cell[state].score {
        cell[state] = Extension {
            score,
            maximum: source.maximum.max(score),
            path: usize::MAX,
            entry: usize::MAX,
        };
    }
}

struct EstSignals {
    values: Vec<[Option<Option<Score>>; 4]>,
    force_gtag: bool,
}
impl EstSignals {
    fn get(&mut self, bases: &[u8], j: usize, lane: usize, close: bool) -> Option<Score> {
        let slot = 2 * lane + usize::from(close);
        if let Some(value) = self.values[j][slot] {
            return value;
        }
        let (opening, closing) = est_splice_types(lane == 1);
        let value = splice_score(
            bases,
            j,
            if close { closing } else { opening },
            self.force_gtag,
        );
        self.values[j][slot] = Some(value);
        value
    }
}

fn est_boundary(
    query: &[u8],
    target: &[u8],
    hsps: &[Hsp],
    intron: IntronScoring,
    signals: &mut EstSignals,
) -> Vec<Vec<usize>> {
    let mut cells: Frontier = Frontier::new(query.len());
    let mut seeds = vec![Vec::new(); target.len() + 1];
    for h in hsps {
        seeds[h.target_start + h.centre].push((h.query_start + h.centre, h.score / 2));
    }
    let mut boundary = vec![Vec::new(); target.len() + 1];
    for j in (0..=target.len()).rev() {
        for &(i, score) in &seeds[j] {
            let cell = cells.at(j, i);
            for state in [0, 4] {
                if score > cell[state].score {
                    cell[state] = Extension {
                        score,
                        maximum: score,
                        path: usize::MAX,
                        entry: usize::MAX,
                    };
                }
            }
        }
        let mut cursor = cells.active[0].len() - 1;
        while let Some((i, mut cell)) = cells.pop(j, &mut cursor, true) {
            for lane in 0..2 {
                let m = lane * 4;
                let ins = m + 1;
                let del = m + 2;
                let span = m + 3;
                // Reverse the gap-return epsilon edges before advancing.
                for state in [ins, del] {
                    if cell[m].score > cell[state].score {
                        cell[state] = cell[m];
                    }
                }
                if (cell[m].score >= 0 || cell[span].score > 0) && boundary[j].last() != Some(&i) {
                    boundary[j].push(i);
                }
                if i > 0 && j > 0 {
                    update_reverse(
                        &mut cells,
                        (j - 1, i - 1),
                        m,
                        cell[m],
                        dna_substitution_score(query[i - 1], target[j - 1]),
                    );
                }
                if i > 0 {
                    update_reverse(&mut cells, (j, i - 1), m, cell[ins], -12);
                    update_reverse(&mut cells, (j, i - 1), ins, cell[ins], -4);
                }
                if j > 0 {
                    update_reverse(&mut cells, (j - 1, i), m, cell[del], -12);
                    update_reverse(&mut cells, (j - 1, i), del, cell[del], -4);
                }
                if j >= 2 {
                    // Closing transitions carry a shadow in C4, so the reverse
                    // scheduler assigns them zero and does not check the length.
                    update_reverse(&mut cells, (j - 2, i), span, cell[m], 0);
                    if let Some(open) = (cell[span].score != NEG_INF)
                        .then(|| signals.get(target, j - 2, lane, false))
                        .flatten()
                    {
                        update_reverse(
                            &mut cells,
                            (j - 2, i),
                            m,
                            cell[span],
                            intron.open_penalty + open,
                        );
                    }
                }
            }
        }
        boundary[j].reverse();
    }
    boundary
}

fn forward_value(source: Extension, delta: Score) -> Option<Extension> {
    if source.score == NEG_INF {
        return None;
    }
    let score = source.score + delta;
    (score >= 0 && source.maximum - score <= 50).then_some(Extension {
        score,
        maximum: source.maximum.max(score),
        path: source.path,
        entry: source.entry,
    })
}

fn record_path(paths: &mut Vec<PathNode>, mut value: Extension, step: RawStep) -> Extension {
    paths.push(PathNode {
        previous: value.path,
        step,
    });
    value.path = paths.len() - 1;
    value
}

fn forward_transition<const ROWS: usize>(
    cells: &mut Frontier<ROWS>,
    paths: &mut Vec<PathNode>,
    position: (usize, usize),
    state: usize,
    source: Extension,
    step: RawStep,
) {
    if let Some(value) = forward_value(source, step.score) {
        let cell = cells.at(position.0, position.1);
        if value.score > cell[state].score {
            cell[state] = record_path(paths, value, step);
        }
    }
}

pub(super) fn align_est(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    strand: Strand,
    intron: IntronScoring,
) -> Option<Alignment> {
    // Analysis_report_func normalises a reverse-query comparison by reversing
    // both sequences. HSPset_revcomp updates starts but preserves cobs; finding
    // seeds afresh on the normalised sequences changes the sparse boundary.
    let hsps = if strand == Strand::Reverse {
        let reverse_query = reverse_complement(&query.bases);
        let mut hsps = find(&reverse_query, &target.bases, false, false);
        for h in &mut hsps {
            h.query_start = query.bases.len() - h.query_start - h.length;
            h.target_start = bases.len() - h.target_start - h.length;
        }
        hsps
    } else {
        find(&query.bases, bases, false, false)
    };
    if hsps.is_empty() {
        return None;
    }
    let intron = canonical_intron_scoring(intron);
    let mut signals = EstSignals {
        values: vec![[None; 4]; bases.len() + 1],
        force_gtag: intron.force_gtag,
    };
    let boundary = est_boundary(&query.bases, bases, &hsps, intron, &mut signals);
    let mut cells: Frontier = Frontier::new(query.bases.len());
    let mut paths: Vec<PathNode> = Vec::new();
    let mut spans: Vec<[Vec<Span>; 2]> = (0..=query.bases.len())
        .map(|_| Default::default())
        .collect();
    let mut end = (0, 0, 0, Extension::default());
    for j in 0..=bases.len() {
        for &i in &boundary[j] {
            cells.at(j, i);
        }
        let mut cursor = 0;
        while let Some((i, mut cell)) = cells.pop(j, &mut cursor, false) {
            let seeded = boundary[j].binary_search(&i).is_ok();
            for lane in 0..2 {
                let m = lane * 4;
                let ins = m + 1;
                let del = m + 2;
                let span = m + 3;
                let id = (lane * 12) as u16;
                if cell[span].score >= 0 {
                    // Freeze only the cell reached by a donor, not a thawed copy.
                    let node = paths[cell[span].path];
                    if node.step.transition_id == id + 6 {
                        spans[i][lane].push(Span {
                            entry: j,
                            value: cell[span],
                        });
                    }
                }
                if seeded && cell[m].score < 0 {
                    cell[m] = record_path(
                        &mut paths,
                        Extension {
                            score: 0,
                            maximum: 0,
                            path: usize::MAX,
                            entry: usize::MAX,
                        },
                        RawStep {
                            transition_id: id,
                            query_advance: 0,
                            target_advance: 0,
                            score: 0,
                        },
                    );
                }
                for (state, edge) in [(ins, 9), (del, 10)] {
                    if cell[state].score > cell[m].score {
                        cell[m] = record_path(
                            &mut paths,
                            cell[state],
                            RawStep {
                                transition_id: id + edge,
                                query_advance: 0,
                                target_advance: 0,
                                score: 0,
                            },
                        );
                    }
                }
                if endpoint_key(cell[m].score, i, j, 0) > endpoint_key(end.3.score, end.1, end.0, 0)
                {
                    end = (j, i, lane, cell[m]);
                }
                if i < query.bases.len() && j < bases.len() {
                    forward_transition(
                        &mut cells,
                        &mut paths,
                        (j + 1, i + 1),
                        m,
                        cell[m],
                        RawStep {
                            transition_id: id + 1,
                            query_advance: 1,
                            target_advance: 1,
                            score: dna_substitution_score(query.bases[i], bases[j]),
                        },
                    );
                }
                if i < query.bases.len() {
                    for (state, edge, score) in [(m, 2, -12), (ins, 4, -4)] {
                        forward_transition(
                            &mut cells,
                            &mut paths,
                            (j, i + 1),
                            ins,
                            cell[state],
                            RawStep {
                                transition_id: id + edge,
                                query_advance: 1,
                                target_advance: 0,
                                score,
                            },
                        );
                    }
                }
                if j < bases.len() {
                    for (state, edge, score) in [(m, 3, -12), (del, 5, -4)] {
                        forward_transition(
                            &mut cells,
                            &mut paths,
                            (j + 1, i),
                            del,
                            cell[state],
                            RawStep {
                                transition_id: id + edge,
                                query_advance: 0,
                                target_advance: 1,
                                score,
                            },
                        );
                    }
                }
                if j + 2 <= bases.len() {
                    if let Some(score) = (cell[m].score != NEG_INF)
                        .then(|| signals.get(bases, j, lane, false))
                        .flatten()
                    {
                        forward_transition(
                            &mut cells,
                            &mut paths,
                            (j + 2, i),
                            span,
                            cell[m],
                            RawStep {
                                transition_id: id + 6,
                                query_advance: 0,
                                target_advance: 2,
                                score: intron.open_penalty + score,
                            },
                        );
                    }
                    if seeded {
                        if let Some(candidate) = spans[i][lane]
                            .iter()
                            .filter(|s| {
                                j + 4 - s.entry >= intron.min_len as usize
                                    && j + 4 - s.entry <= intron.max_len as usize
                            })
                            .max_by_key(|s| (s.value.score, std::cmp::Reverse(s.entry)))
                            .copied()
                        {
                            if candidate.value.score > cell[span].score {
                                cell[span] = record_path(
                                    &mut paths,
                                    candidate.value,
                                    RawStep {
                                        transition_id: id + 7,
                                        query_advance: 0,
                                        target_advance: (j - candidate.entry) as u32,
                                        score: 0,
                                    },
                                );
                            }
                        }
                    }
                    if let Some(score) = (cell[span].score != NEG_INF)
                        .then(|| signals.get(bases, j, lane, true))
                        .flatten()
                    {
                        if cell[span].score != NEG_INF {
                            // The stored opening node identifies the shadow's donor.
                            let mut node = cell[span].path;
                            let mut length = 2usize;
                            while node != usize::MAX {
                                let step = paths[node].step;
                                length += step.target_advance as usize;
                                if step.transition_id == id + 6 {
                                    break;
                                }
                                node = paths[node].previous;
                            }
                            if length >= intron.min_len as usize
                                && length <= intron.max_len as usize
                            {
                                forward_transition(
                                    &mut cells,
                                    &mut paths,
                                    (j + 2, i),
                                    m,
                                    cell[span],
                                    RawStep {
                                        transition_id: id + 8,
                                        query_advance: 0,
                                        target_advance: 2,
                                        score,
                                    },
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    if end.3.score <= 0 {
        return None;
    }
    let mut raw = Vec::new();
    let mut node = end.3.path;
    while node != usize::MAX {
        raw.push(paths[node].step);
        node = paths[node].previous;
    }
    raw.reverse();
    let qa: u64 = raw.iter().map(|s| u64::from(s.query_advance)).sum();
    let ta: u64 = raw.iter().map(|s| u64::from(s.target_advance)).sum();
    let mut trace = Vec::new();
    for step in &raw {
        let op = match step.transition_id % 12 {
            1 => Op::Match,
            2 | 4 => Op::Insert,
            3 | 5 => Op::Delete,
            6 => {
                if end.2 == 0 {
                    Op::Splice5
                } else {
                    Op::Splice3
                }
            }
            7 => Op::Intron,
            8 => {
                if end.2 == 0 {
                    Op::Splice3
                } else {
                    Op::Splice5
                }
            }
            _ => continue,
        };
        TraceFragment::one(TraceRun {
            transition_id: step.transition_id,
            op,
            query_advance: step.query_advance,
            target_advance: step.target_advance,
            repeats: 1,
        })
        .append_to(&mut trace);
    }
    raw.push(RawStep {
        transition_id: (end.2 * 12 + 11) as u16,
        query_advance: 0,
        target_advance: 0,
        score: 0,
    });
    let (ts, te) = (end.0 as u64 - ta, end.0 as u64);
    Some(Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: end.1 as u64 - qa,
        query_end: end.1 as u64,
        query_strand: Strand::Forward,
        target_start: if strand == Strand::Forward {
            ts
        } else {
            bases.len() as u64 - ts
        },
        target_end: if strand == Strand::Forward {
            te
        } else {
            bases.len() as u64 - te
        },
        target_len: bases.len() as u64,
        target_strand: strand,
        score: end.3.score,
        trace,
        raw_trace: raw,
    })
}

fn protein_hsps(query: &[u8], bases: &[u8]) -> Option<Vec<Hsp>> {
    if query
        .iter()
        .any(|b| !b"ACDEFGHIKLMNPQRSTVWY".contains(&b.to_ascii_uppercase()))
    {
        return None;
    }
    let mut hsps = Vec::new();
    for frame in 0..3 {
        for mut h in find(query, &translate_dna(bases, frame), true, false) {
            h.target_start = frame as usize + 3 * h.target_start;
            hsps.push(h);
        }
    }
    hsps.sort_by_key(|h| (h.target_start + 3 * h.centre, h.query_start + h.centre));
    Some(hsps)
}

fn protein_boundary(
    query: &[u8],
    bases: &[u8],
    translated: &[u8],
    hsps: &[Hsp],
    scoring: Scoring,
    intron: IntronScoring,
    signals: &mut EstSignals,
) -> Vec<Vec<usize>> {
    let mut cells = Frontier::<6>::new(query.len());
    let mut seeds = vec![Vec::new(); bases.len() + 1];
    for h in hsps {
        seeds[h.target_start + 3 * h.centre].push((h.query_start + h.centre, h.score / 2));
    }
    let mut boundary = vec![Vec::new(); bases.len() + 1];
    for j in (0..=bases.len()).rev() {
        for &(i, score) in &seeds[j] {
            let cell = cells.at(j, i);
            if score > cell[0].score {
                cell[0] = Extension {
                    score,
                    maximum: score,
                    path: usize::MAX,
                    entry: usize::MAX,
                };
            }
        }
        let mut cursor = cells.active[0].len() - 1;
        while let Some((i, mut cell)) = cells.pop(j, &mut cursor, true) {
            for state in [1, 2] {
                if cell[0].score > cell[state].score {
                    cell[state] = cell[0];
                }
            }
            if cell[0].score >= 0 || (3..6).any(|state| cell[state].score > 0) {
                boundary[j].push(i);
            }
            if i > 0 && j >= 3 {
                update_reverse(
                    &mut cells,
                    (j - 3, i - 1),
                    0,
                    cell[0],
                    protein_substitution_score(query[i - 1], translated[j - 3]),
                );
            }
            if i > 0 {
                update_reverse(&mut cells, (j, i - 1), 0, cell[1], scoring.codon_gap_open);
                update_reverse(&mut cells, (j, i - 1), 1, cell[1], scoring.codon_gap_extend);
            }
            if j >= 3 {
                update_reverse(&mut cells, (j - 3, i), 0, cell[2], scoring.codon_gap_open);
                update_reverse(&mut cells, (j - 3, i), 2, cell[2], scoring.codon_gap_extend);
            }
            for advance in [1, 2, 4, 5] {
                if j >= advance {
                    update_reverse(&mut cells, (j - advance, i), 0, cell[0], scoring.frameshift);
                }
            }
            for phase in 0..3 {
                if phase == 0 && j >= 2 {
                    update_reverse(&mut cells, (j - 2, i), 3, cell[0], 0);
                }
                if phase > 0 && i > 0 && j >= 5 - phase {
                    update_reverse(&mut cells, (j - (5 - phase), i - 1), 3 + phase, cell[0], 0);
                }
                if j >= 2 + phase && cell[3 + phase].score != NEG_INF {
                    if let Some(score) = signals.get(bases, j - 2, 0, false) {
                        update_reverse(
                            &mut cells,
                            (j - 2 - phase, i),
                            0,
                            cell[3 + phase],
                            intron.open_penalty + score,
                        );
                    }
                }
            }
        }
        boundary[j].reverse();
    }
    boundary
}

fn forward_steps<const ROWS: usize>(
    cells: &mut Frontier<ROWS>,
    paths: &mut Vec<PathNode>,
    position: (usize, usize),
    state: usize,
    source: Extension,
    steps: &[RawStep],
) {
    let mut value = source;
    for step in steps {
        let Some(next) = forward_value(value, step.score) else {
            return;
        };
        value = next;
    }
    let cell = cells.at(position.0, position.1);
    if value.score <= cell[state].score {
        return;
    }
    for step in steps {
        value = record_path(paths, value, *step);
    }
    cell[state] = value;
}

// Outer None requests the existing backend for unsupported alphabets. An empty
// inner result means this supported seeded search found no positive alignment.
#[allow(clippy::too_many_arguments)]
pub(super) fn align_protein(
    query: &Sequence,
    target: &Sequence,
    bases: &[u8],
    strand: Strand,
    scoring: Scoring,
    intron: IntronScoring,
) -> Option<Option<Alignment>> {
    let hsps = protein_hsps(&query.bases, bases)?;
    if hsps.is_empty() {
        return Some(None);
    }
    let intron = canonical_intron_scoring(intron);
    let translated: Vec<_> = bases.windows(3).map(translate_codon).collect();
    let mut signals = EstSignals {
        values: vec![[None; 4]; bases.len() + 1],
        force_gtag: intron.force_gtag,
    };
    let boundary = protein_boundary(
        &query.bases,
        bases,
        &translated,
        &hsps,
        scoring,
        intron,
        &mut signals,
    );
    let mut cells = Frontier::<6>::new(query.bases.len());
    let mut paths: Vec<PathNode> = Vec::new();
    let mut spans: Vec<[Vec<Span>; 3]> = (0..=query.bases.len())
        .map(|_| Default::default())
        .collect();
    let mut end = (0, 0, Extension::default());
    let step = |id, qa, ta, score| RawStep {
        transition_id: id,
        query_advance: qa,
        target_advance: ta,
        score,
    };
    for j in 0..=bases.len() {
        for &i in &boundary[j] {
            cells.at(j, i);
        }
        let mut cursor = 0;
        while let Some((i, mut cell)) = cells.pop(j, &mut cursor, false) {
            let seeded = boundary[j].binary_search(&i).is_ok();
            if seeded && cell[0].score < 0 {
                cell[0] = record_path(
                    &mut paths,
                    Extension {
                        score: 0,
                        maximum: 0,
                        path: usize::MAX,
                        entry: usize::MAX,
                    },
                    step(0, 0, 0, 0),
                );
            }
            for (state, id) in [(1, 6), (2, 7)] {
                if cell[state].score > cell[0].score {
                    cell[0] = record_path(&mut paths, cell[state], step(id, 0, 0, 0));
                }
            }
            if endpoint_key(cell[0].score, i, j, 0) > endpoint_key(end.2.score, end.1, end.0, 0) {
                end = (j, i, cell[0]);
            }
            if i < query.bases.len() && j + 3 <= bases.len() {
                forward_transition(
                    &mut cells,
                    &mut paths,
                    (j + 3, i + 1),
                    0,
                    cell[0],
                    step(
                        1,
                        1,
                        3,
                        protein_substitution_score(query.bases[i], translated[j]),
                    ),
                );
            }
            if i < query.bases.len() {
                for (state, id, score) in [
                    (0, 2, scoring.codon_gap_open),
                    (1, 4, scoring.codon_gap_extend),
                ] {
                    forward_transition(
                        &mut cells,
                        &mut paths,
                        (j, i + 1),
                        1,
                        cell[state],
                        step(id, 1, 0, score),
                    );
                }
            }
            if j + 3 <= bases.len() {
                for (state, id, score) in [
                    (0, 3, scoring.codon_gap_open),
                    (2, 5, scoring.codon_gap_extend),
                ] {
                    forward_transition(
                        &mut cells,
                        &mut paths,
                        (j + 3, i),
                        2,
                        cell[state],
                        step(id, 0, 3, score),
                    );
                }
            }
            for advance in [1usize, 2, 4, 5] {
                if j + advance <= bases.len() {
                    let short = advance % 3;
                    forward_steps(
                        &mut cells,
                        &mut paths,
                        (j + advance, i),
                        0,
                        cell[0],
                        &[
                            step(
                                if short == 1 { 8 } else { 9 },
                                0,
                                short as u32,
                                scoring.frameshift,
                            ),
                            step(
                                if advance > 3 { 11 } else { 10 },
                                0,
                                if advance > 3 { 3 } else { 0 },
                                0,
                            ),
                        ],
                    );
                }
            }
            for phase in 0..3 {
                let state = 3 + phase;
                if cell[state].score >= 0 {
                    spans[i][phase].push(Span {
                        entry: j,
                        value: cell[state],
                    });
                }
                if j + phase + 2 <= bases.len() && cell[0].score != NEG_INF {
                    if let Some(score) = signals.get(bases, j + phase, 0, false) {
                        let mut source = cell[0];
                        source.entry = j + phase + 2;
                        let opening = step(110 + phase as u16, 0, 2, intron.open_penalty + score);
                        if phase == 0 {
                            forward_transition(
                                &mut cells,
                                &mut paths,
                                (j + 2, i),
                                state,
                                source,
                                opening,
                            );
                        } else {
                            forward_steps(
                                &mut cells,
                                &mut paths,
                                (j + phase + 2, i),
                                state,
                                source,
                                &[step(100 + phase as u16, 0, phase as u32, 0), opening],
                            );
                        }
                    }
                }
                if seeded {
                    if let Some(candidate) = spans[i][phase]
                        .iter()
                        .filter(|s| {
                            j + 4 - s.entry >= intron.min_len as usize
                                && j + 4 - s.entry <= intron.max_len as usize
                        })
                        .max_by_key(|s| (s.value.score, std::cmp::Reverse(s.entry)))
                        .copied()
                    {
                        if candidate.value.score > cell[state].score {
                            cell[state] = record_path(
                                &mut paths,
                                candidate.value,
                                step(120 + phase as u16, 0, (j - candidate.entry) as u32, 0),
                            );
                        }
                    }
                }
                if j + 2 <= bases.len() && cell[state].score != NEG_INF {
                    let length = j + 4 - cell[state].entry;
                    if length >= intron.min_len as usize && length <= intron.max_len as usize {
                        if let Some(score) = signals.get(bases, j, 0, true) {
                            forward_transition(
                                &mut cells,
                                &mut paths,
                                (j + 2, i),
                                if phase == 0 { 0 } else { 5 + phase },
                                cell[state],
                                step(130 + phase as u16, 0, 2, score),
                            );
                        }
                    }
                }
                if phase > 0
                    && cell[5 + phase].score != NEG_INF
                    && i < query.bases.len()
                    && j + 3 - phase <= bases.len()
                {
                    let source = cell[5 + phase];
                    let donor = source.entry - 2;
                    let mut codon = [0u8; 3];
                    codon[..phase].copy_from_slice(&bases[donor - phase..donor]);
                    codon[phase..].copy_from_slice(&bases[j..j + 3 - phase]);
                    let score = protein_substitution_score(query.bases[i], translate_codon(&codon));
                    forward_transition(
                        &mut cells,
                        &mut paths,
                        (j + 3 - phase, i + 1),
                        0,
                        source,
                        step(140 + phase as u16, 1, (3 - phase) as u32, score),
                    );
                }
            }
        }
    }
    if end.2.score <= 0 {
        return Some(None);
    }
    let mut raw = Vec::new();
    let mut node = end.2.path;
    while node != usize::MAX {
        raw.push(paths[node].step);
        node = paths[node].previous;
    }
    raw.reverse();
    let qa: u64 = raw.iter().map(|s| u64::from(s.query_advance)).sum();
    let ta: u64 = raw.iter().map(|s| u64::from(s.target_advance)).sum();
    let mut trace = Vec::new();
    for s in &raw {
        let op = match s.transition_id {
            1 => Op::Match,
            2 | 4 => Op::Insert,
            3 | 5 => Op::Delete,
            8 | 9 | 11 => Op::Frameshift,
            101 | 102 | 141 | 142 => Op::SplitCodon,
            110..=112 => Op::Splice5,
            120..=122 => Op::Intron,
            130..=132 => Op::Splice3,
            _ => continue,
        };
        TraceFragment::one(TraceRun {
            transition_id: s.transition_id,
            op,
            query_advance: s.query_advance,
            target_advance: s.target_advance,
            repeats: 1,
        })
        .append_to(&mut trace);
    }
    let mut raw_trace: Vec<_> = raw
        .into_iter()
        .flat_map(|s| {
            if (120..=122).contains(&s.transition_id) {
                std::iter::repeat_n(
                    RawStep {
                        target_advance: 1,
                        ..s
                    },
                    s.target_advance as usize,
                )
            } else {
                std::iter::repeat_n(s, 1)
            }
        })
        .collect();
    raw_trace.push(step(12, 0, 0, 0));
    let (ts, te) = (end.0 as u64 - ta, end.0 as u64);
    Some(Some(Alignment {
        query_id: query.id.clone(),
        target_id: target.id.clone(),
        query_start: end.1 as u64 - qa,
        query_end: end.1 as u64,
        query_strand: Strand::Unknown,
        target_start: if strand == Strand::Forward {
            ts
        } else {
            bases.len() as u64 - ts
        },
        target_end: if strand == Strand::Forward {
            te
        } else {
            bases.len() as u64 - te
        },
        target_len: bases.len() as u64,
        target_strand: strand,
        score: end.2.score,
        trace,
        raw_trace,
    }))
}
