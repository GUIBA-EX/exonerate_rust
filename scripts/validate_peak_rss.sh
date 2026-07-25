#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
sequence_length=500
if (( $# != 0 )); then
    echo "usage: $0" >&2
    exit 2
fi
if [[ ! -x /usr/bin/time ]]; then
    echo "/usr/bin/time is required for peak-RSS validation" >&2
    exit 2
fi

work_dir=$(mktemp -d /tmp/exonerate-peak-rss.XXXXXX)
trap 'rm -r "$work_dir"' EXIT

awk -v seq_len="$sequence_length" '
    BEGIN {
        print ">q"
        state = 17
        for (i = 0; i < seq_len; i++) {
            state = (state * 25173 + 13849) % 65536
            printf substr("ACGT", (int(state / 256) % 4) + 1, 1)
        }
        print ""
    }
' > "$work_dir/query.fa"
sed 's/^>q$/>t/' "$work_dir/query.fa" > "$work_dir/target.fa"

cd "$repo_root"
cargo build --release --locked -p exonerate

common=(
    --verbose 0
    --model ner
    --exhaustive yes
    --subopt no
    --forwardonly
    --score 0
    --showalignment no
    --showsugar yes
    --showcigar yes
    --showvulgar yes
    --ryo '{%Pn|%Ps|%Pqa|%Pta\n}'
)
binary=target/release/exonerate-rs

# At the fixed 500x500 fixture, the generic full-table planner is safely
# below 4096 MiB.  Keeping the fixture fixed prevents this reference run from
# silently selecting the checkpoint backend.
/usr/bin/time -f '%M' -o "$work_dir/full.rss" \
    "$binary" "${common[@]}" --dpmemory 4096 \
    "$work_dir/query.fa" "$work_dir/target.fa" > "$work_dir/full.out"
/usr/bin/time -f '%M' -o "$work_dir/checkpoint.rss" \
    "$binary" "${common[@]}" --dpmemory 0 \
    "$work_dir/query.fa" "$work_dir/target.fa" > "$work_dir/checkpoint.out"

cmp "$work_dir/full.out" "$work_dir/checkpoint.out"
full_rss=$(<"$work_dir/full.rss")
checkpoint_rss=$(<"$work_dir/checkpoint.rss")
if (( checkpoint_rss >= full_rss )); then
    echo "checkpoint RSS (${checkpoint_rss} KiB) was not below full-table RSS (${full_rss} KiB)" >&2
    exit 1
fi

echo "NER length:       $sequence_length bases per sequence"
echo "full-table RSS:   $full_rss KiB"
echo "checkpoint RSS:   $checkpoint_rss KiB"
echo "outputs:          byte-identical"
