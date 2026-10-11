#!/usr/bin/env python3
"""Seeded differential audit; requires built Rust and upstream executables.

This is deliberately separate from the passing regression suite: exit status 1
means observed incompatibilities, all of which are retained in the JSON report.
No downloads or third-party Python packages are required.
"""
import argparse
import collections
import json
import pathlib
import random
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODELS = [
    "ungapped", "ungapped:trans", "affine:global", "affine:bestfit",
    "affine:local", "affine:overlap", "ner", "est2genome", "protein2dna",
    "protein2dna:bestfit", "protein2genome", "protein2genome:bestfit",
    "coding2coding", "coding2genome", "cdna2genome", "genome2genome",
]
CODONS = dict(zip("ACDEFGHIKLMNPQRSTVWY", "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))


def execute(binary, args, timeout):
    try:
        result = subprocess.run([str(binary), *args], capture_output=True, text=True, timeout=timeout)
        return {"returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}
    except subprocess.TimeoutExpired:
        return {"returncode": None, "stdout": "", "stderr": "audit timeout"}


def same(left, right):
    return (left["returncode"], left["stdout"]) == (right["returncode"], right["stdout"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rust", type=pathlib.Path, default=ROOT / "target/debug/exonerate-rs")
    parser.add_argument("--upstream", type=pathlib.Path, default=ROOT / "upstream/src/program/exonerate")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--seed", type=int, default=711)
    parser.add_argument("--cases", type=int, default=16)
    parser.add_argument("--suboptimal", action="store_true", help="also audit suboptimal enumeration (some upstream inputs may time out)")
    parser.add_argument("--timeout", type=float, default=15)
    parser.add_argument("--ambiguous", action="store_true", help="include N and lengths 1–64 for translation boundary checks")
    parser.add_argument("--heuristic", action="store_true", help="audit default seeded EST/protein spliced extension")
    parser.add_argument("--raw-steps", action="store_true", help="compare atomic RYO transition names and scores as well")
    options = parser.parse_args()
    if options.cases < 1 or options.timeout <= 0:
        parser.error("cases and timeout must be positive")
    for binary in (options.rust, options.upstream):
        if not binary.is_file():
            parser.error(f"build the executable first: {binary}")
    rng = random.Random(options.seed)
    dna_alphabet = "ACGTN" if options.ambiguous else "ACGT"
    summary = {}
    failures = []
    suites = [(model, False) for model in MODELS]
    suites += [(model, True) for model in MODELS if model == "ungapped" or model.startswith("affine:")]
    if options.heuristic:
        suites = [(model, False) for model in ("est2genome", "protein2genome")]
    with tempfile.TemporaryDirectory(prefix="exonerate-audit-") as directory:
        query_path, target_path = pathlib.Path(directory) / "q.fa", pathlib.Path(directory) / "t.fa"
        for model, protein_pair in suites:
            name = model + (":protein2protein" if protein_pair else "")
            counts = collections.Counter(cases=0, upstream_differences=0, checkpoint_differences=0, execution_errors=0)
            for case in range(options.cases):
                dna = "".join(rng.choices(dna_alphabet, k=rng.randrange(1, 65) if options.ambiguous else rng.randrange(9, 31)))
                query = "".join(rng.choices(list(CODONS), k=rng.randrange(3, 10))) if model.startswith("protein") or protein_pair else dna
                target = query if protein_pair else "".join(CODONS[x] for x in query) if model.startswith("protein") else dna
                if case % 4 == 0:
                    target = "".join(rng.choices(list(CODONS) if protein_pair else dna_alphabet, k=rng.randrange(1, 65) if options.ambiguous else rng.randrange(9, 35)))
                elif case % 4 == 1:
                    target = "GG" + target[:len(target)//2] + "AC" + target[len(target)//2:] + "TT"
                elif case % 4 == 2 and not protein_pair:
                    target = target.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                if options.heuristic:
                    query = "".join(rng.choices(list(CODONS) if model.startswith("protein") else "ACGT", k=80))
                    target = "".join(CODONS[x] for x in query) if model.startswith("protein") else query
                    cut = len(target) // 2 + case % 3
                    intron = "GT" + "".join(rng.choices("ACGT", k=60 + case)) + "AG"
                    if case % 6 < 3:
                        target = target[:cut] + intron + target[cut:]
                    elif case % 6 == 3:
                        target = target[:cut] + "A" + target[cut:]
                    elif case % 6 == 4:
                        target = target[:cut] + "CC" + target[cut:]
                    else:
                        target = target[:cut] + target[cut + 3:]
                    target = "ACGTT" + target + "TGCAT"
                    if case % 2:
                        target = target.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                query_path.write_text(">q\n" + query + "\n")
                target_path.write_text(">t\n" + target + "\n")
                for subopt in (("no", "yes") if options.suboptimal else ("no",)):
                    for revcomp in (("no",) if protein_pair else ("no", "yes")):
                        args = ["--model", model, "--verbose", "0", "--exhaustive", "no" if options.heuristic else "yes", "--subopt", subopt, "--revcomp", revcomp, "--score", "1", "--showalignment", "no", "--showsugar", "yes", "--showcigar", "yes", "--showvulgar", "yes"]
                        if protein_pair or model.startswith("protein"):
                            args += ["--querytype", "protein", "--targettype", "protein" if protein_pair else "dna"]
                        if options.heuristic:
                            args += ["--maxintron", "50" if case % 7 == 6 else "200000", "--forcegtag", "yes" if case % 5 == 4 else "no"]
                        if options.raw_steps:
                            args += ["--ryo", "{%Pn|%Ps|%Pqa|%Pta\\n}"]
                        inputs = [str(query_path), str(target_path)]
                        rust = execute(options.rust.resolve(), [*args, *inputs], options.timeout)
                        low = execute(options.rust.resolve(), [*args, "--dpmemory", "0", *inputs], options.timeout)
                        upstream = execute(options.upstream.resolve(), [*args, *inputs], options.timeout)
                        counts["cases"] += 1
                        counts["upstream_differences"] += not same(rust, upstream)
                        counts["checkpoint_differences"] += not same(rust, low)
                        counts["execution_errors"] += any(r["returncode"] != 0 for r in (rust, low, upstream))
                        if not same(rust, upstream) or not same(rust, low) or any(r["returncode"] != 0 for r in (rust, low, upstream)):
                            failures.append({"suite": name, "case": case, "query": query, "target": target, "args": args, "rust": rust, "checkpoint": low, "upstream": upstream})
            summary[name] = dict(counts)
            print(name, dict(counts), flush=True)
    report = {"seed": options.seed, "cases_per_suite": options.cases, "ambiguous": options.ambiguous, "heuristic": options.heuristic, "raw_steps": options.raw_steps, "summary": summary, "failures": failures}
    options.output.write_text(json.dumps(report, indent=2) + "\n")
    print("Report:", options.output)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
