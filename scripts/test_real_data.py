#!/usr/bin/env python3
"""Bounded comparisons on bundled biological FASTA and pinned RefSeqGene data.

Build release binaries first. Inputs and provenance: docs/audit/real-data/inputs.
No downloads or third-party Python packages. Linux /usr/bin/time required.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import resource
import signal
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sequence_info(path):
    data = path.read_bytes()
    return {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(data).hexdigest(),
            "length": sum(len(line.strip()) for line in data.splitlines() if not line.startswith(b">"))}


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def execute(binary, args, directory, label, timeout):
    stdout, stderr, metrics = [directory / (label + suffix) for suffix in (".out", ".err", ".time")]
    command = ["/usr/bin/time", "-f", "%e %M", "-o", str(metrics), str(binary), *args]
    start = time.perf_counter()
    timed_out = False
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True, preexec_fn=limits)
    try:
        out, err = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGKILL)
        out, err = process.communicate()
    elapsed = time.perf_counter() - start
    stdout.write_bytes(out)
    stderr.write_bytes(err)
    result = {"returncode": process.returncode, "timeout": timed_out,
              "wall_seconds": round(elapsed, 6),
              "peak_rss_kib": None, "command": command,
              "stdout": out.decode(), "stderr": err.decode()}
    if metrics.exists():
        lines = metrics.read_text().splitlines()
        if lines and len(lines[-1].split()) == 2:
            elapsed, rss = lines[-1].split()
            result.update(time_wall_seconds=float(elapsed), peak_rss_kib=int(rss))
    return result



def batch_test(output, timeout):
    directory = output.parent / (output.stem + "-batch")
    directory.mkdir(parents=True, exist_ok=True)
    genes = ["calm.human", "p53.human", "htrt.human", "tube.drome"]
    header = "task_id sample_id locus_id candidate_id expected_family reference_family model query_fasta query_id target_fasta target_id".split()
    rows = ["\t".join(header)]
    for candidate in genes:
        target = ROOT / f"upstream/test/data/cdna/{candidate}.dna.fasta"
        for reference in genes:
            query = ROOT / f"upstream/test/data/protein/{reference}.protein.fasta"
            rows.append("\t".join([candidate + "." + reference, "real-data", candidate, candidate,
                                   candidate, reference, "protein2dna", str(query),
                                   query.read_text().splitlines()[0][1:].split()[0], str(target),
                                   target.read_text().splitlines()[0][1:].split()[0]]))
    manifest = directory / "tasks.tsv"
    manifest.write_text("\n".join(rows) + "\n")
    results = {}
    documents = {}
    for threads in [1, 4]:
        paths = {name: directory / f"threads-{threads}.{suffix}" for name, suffix in
                 [("alignments", "tsv"), ("orthology", "orthology.tsv"), ("evidence", "gff3")]}
        args = ["--tasks", str(manifest), "--threads", str(threads),
                "--audit", "protein-candidate", "--subopt", "no",
                "--score", "100", "--revcomp", "yes",
                "--result-tsv", str(paths["alignments"]),
                "--orthology-report", str(paths["orthology"]),
                "--evidence-gff3", str(paths["evidence"])]
        results[str(threads)] = execute(ROOT / "target/release/exonerate-rs", args, directory, f"threads-{threads}", timeout)
        documents[str(threads)] = {name: path.read_text() if path.exists() else None for name, path in paths.items()}
    succeeded = all(r["returncode"] == 0 for r in results.values())
    report = {"tasks": 16, "threads": [1, 4], "all_succeeded": succeeded,
              "byte_identical_reports": succeeded and documents["1"] == documents["4"],
              "runs": results, "reports": documents}
    if succeeded:
        table = documents["1"]["alignments"].splitlines()
        report["tsv_column_counts_valid"] = all(len(row.split("\t")) == len(table[0].split("\t")) for row in table[1:])
        report["orthology_statuses"] = [row.split("\t")[4] for row in documents["1"]["orthology"].splitlines()[1:]]
        ids, parents = [], []
        for line in documents["1"]["evidence"].splitlines():
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            assert int(fields[3]) <= int(fields[4])
            for attr in fields[8].split(";"):
                if attr.startswith("ID="):
                    ids.append(attr[3:])
                if attr.startswith("Parent="):
                    parents.append(attr[7:])
        report["gff_ids_unique"] = len(set(ids)) == len(ids)
        report["gff_parents_valid"] = all(parent in set(ids) for parent in parents)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print({k: v for k, v in report.items() if k not in ("runs", "reports")})
    return int(not (succeeded and report["byte_identical_reports"] and report["tsv_column_counts_valid"]
                    and report["gff_ids_unique"] and report["gff_parents_valid"]))



def summarize(report):
    for case in report["cases"]:
        runs = case["runs"]
        for run in runs.values():
            run["status"] = ("timeout" if run["timeout"] else
                             "memory_limit_failure" if "memory allocation of" in run["stderr"] else
                             "ok" if run["returncode"] == 0 else "execution_failure")
        case["all_succeeded"] = all(run["status"] == "ok" and run.get("repetitions_consistent", True) for run in runs.values())
        for field, labels in [("upstream_equal", ("upstream", "rust_64mb")),
                              ("memory_modes_equal", ("rust_64mb", "rust_0mb"))]:
            left, right = (runs[label] for label in labels)
            case[field] = left["stdout"] == right["stdout"] if left["status"] == right["status"] == "ok" else None
        case["speedup"] = {
            label: runs["upstream"]["wall_seconds"] / runs[label]["wall_seconds"]
            if runs[label]["wall_seconds"] > 0 and case["all_succeeded"] else None
            for label in ("rust_64mb", "rust_0mb")
        }
        case["performance_accepted"] = (case["all_succeeded"] and case["upstream_equal"] is True
                                        and case["memory_modes_equal"] is True
                                        and all(value is not None and value > 1 for value in case["speedup"].values()))
    cases = report["cases"]
    runs = [run for case in cases for run in case["runs"].values()]
    report["summary"] = {
        "cases": len(cases), "runs": len(runs),
        "executions": sum(len(r.get("samples", [r])) for r in runs),
        "cases_all_succeeded": sum(c["all_succeeded"] for c in cases),
        "execution_failures": sum(not c["all_succeeded"] for c in cases),
        "upstream_comparisons_completed": sum(c["upstream_equal"] is not None for c in cases),
        "upstream_differences": sum(c["upstream_equal"] is False for c in cases),
        "memory_comparisons_completed": sum(c["memory_modes_equal"] is not None for c in cases),
        "memory_differences": sum(c["memory_modes_equal"] is False for c in cases),
        "timed_out_runs": sum(r["status"] == "timeout" for r in runs),
        "memory_limit_failure_runs": sum(r["status"] == "memory_limit_failure" for r in runs),
        "performance_accepted_cases": sum(c["performance_accepted"] for c in cases),
        "all_cases_faster_and_equal": bool(cases) and all(c["performance_accepted"] for c in cases),
    }



def write_summary(report, path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["case", "query_length", "target_length", "backend", "status",
                         "seconds", "peak_rss_kib", "upstream_equal", "memory_modes_equal"])
        for case in report["cases"]:
            for label, run in case["runs"].items():
                writer.writerow([case["id"], case["query"]["length"], case["target"]["length"],
                                 label, run["status"], run["wall_seconds"], run["peak_rss_kib"],
                                 case["upstream_equal"], case["memory_modes_equal"]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, default=ROOT / "docs/audit/real-data/results.json")
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--case", help="only case IDs containing this substring")
    parser.add_argument("--require-faster", action="store_true", help="fail unless every case matches upstream and both Rust modes have lower median time")
    parser.add_argument("--repeat", type=int, default=1, help="repeat each configuration; report median elapsed time")
    parser.add_argument("--batch-only", action="store_true", help="test 16 real reference/candidate tasks at 1 and 4 threads")
    options = parser.parse_args()
    if options.repeat <= 0:
        parser.error("repeat must be positive")
    if options.timeout <= 0:
        parser.error("timeout must be positive")
    options.output = options.output.resolve()
    if options.batch_only:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        return batch_test(options.output, options.timeout)
    upstream = ROOT / "upstream/src/program/exonerate"
    rust = ROOT / "target/release/exonerate-rs"
    cases = []

    def add(name, model, query, target, protein, revcomp="yes", exhaustive="yes"):
        cases.append((name, model, query, target, protein, revcomp, exhaustive))

    for gene in ["calm.human", "p53.human", "htrt.human", "tube.drome"]:
        protein = ROOT / f"upstream/test/data/protein/{gene}.protein.fasta"
        cdna = ROOT / f"upstream/test/data/cdna/{gene}.dna.fasta"
        for model in ["protein2dna", "protein2dna:bestfit", "protein2genome"]:
            add(gene + "." + model, model, protein, cdna, True)
        add(gene + ".dna-self", "ungapped", cdna, cdna, False, "no")
    for gene, accession in [("calm.human", "NG_013338.1"), ("p53.human", "NG_017013.2")]:
        genomic = ROOT / f"docs/audit/real-data/inputs/{accession}.fasta"
        protein = ROOT / f"upstream/test/data/protein/{gene}.protein.fasta"
        cdna = ROOT / f"upstream/test/data/cdna/{gene}.dna.fasta"
        for mode in ["yes", "no"]:
            add(gene + ".genomic.protein2genome.exhaustive-" + mode, "protein2genome", protein, genomic, True, "yes", mode)
            add(gene + ".genomic.est2genome.exhaustive-" + mode, "est2genome", cdna, genomic, False, "yes", mode)
    report = {"reference": "repository upstream Exonerate 2.4.0", "timeout_seconds": options.timeout,
              "address_space_limit_gib": 8, "rust_binary_sha256": hashlib.sha256(rust.read_bytes()).hexdigest(),
              "upstream_binary_sha256": hashlib.sha256(upstream.read_bytes()).hexdigest(),
              "notes": ["Gene-sized sequences, not whole-genome throughput testing.",
                        "subopt=no, score=100; best hit per orientation.",
                        "dpmemory is not an RSS cap; protein2dna retains full traceback parents.",
                        f"{options.repeat} timings per configuration; wall_seconds is the median, including process launch and /usr/bin/time overhead.",
                        "Heuristic search may differ from exhaustive search."], "cases": []}
    options.output.parent.mkdir(parents=True, exist_ok=True)
    for name, model, query, target, protein, revcomp, exhaustive in cases:
        if options.case and options.case not in name:
            continue
        directory = options.output.parent / (options.output.stem + "-runs") / name
        directory.mkdir(parents=True, exist_ok=True)
        args = ["--model", model, "--querytype", "protein" if protein else "dna", "--targettype", "dna",
                "--verbose", "0", "--exhaustive", exhaustive, "--subopt", "no", "--score", "100",
                "--revcomp", revcomp, "--showalignment", "no", "--showsugar", "yes", "--showcigar", "yes", "--showvulgar", "yes"]
        row = {"id": name, "model": model, "exhaustive": exhaustive,
               "query": sequence_info(query), "target": sequence_info(target), "runs": {}}
        for label, binary, extra in [("upstream", upstream, ["--dpmemory", "64"]),
                                     ("rust_64mb", rust, ["--dpmemory", "64"]),
                                     ("rust_0mb", rust, ["--dpmemory", "0"])]:
            samples = [execute(binary, [*args, *extra, str(query), str(target)], directory,
                               f"{label}-{iteration}", options.timeout)
                       for iteration in range(options.repeat)]
            # Preserve failures, even if a later repetition succeeds.
            run = next((sample.copy() for sample in samples
                        if sample["returncode"] != 0 or sample["timeout"]), samples[0].copy())
            run["wall_seconds"] = statistics.median(sample["wall_seconds"] for sample in samples)
            run["peak_rss_kib"] = max((sample["peak_rss_kib"] for sample in samples
                                      if sample["peak_rss_kib"] is not None), default=None)
            run["repetitions_consistent"] = all(sample["stdout"] == samples[0]["stdout"] for sample in samples)
            run["samples"] = samples
            row["runs"][label] = run
            print(name, label, run["returncode"], run["wall_seconds"], run["peak_rss_kib"], flush=True)
        report["cases"].append(row)
        summarize(report)
        options.output.write_text(json.dumps(report, indent=2) + "\n")
        write_summary(report, options.output.with_suffix(".tsv"))
    if not report["cases"]:
        parser.error("no cases matched --case")
    print(report.get("summary", {}), flush=True)
    return int(any(report.get("summary", {}).get(k, 0) for k in ["execution_failures", "upstream_differences", "memory_differences"])
               or (options.require_faster and not report.get("summary", {}).get("all_cases_faster_and_equal", False)))


if __name__ == "__main__":
    raise SystemExit(main())
