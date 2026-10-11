use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_exonerate-rs"))
        .args(["--verbose", "0"])
        .args(args)
        .output()
        .expect("run exonerate-rs");
    assert!(
        output.status.success(),
        "exonerate-rs failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("CLI output is UTF-8")
}

fn run_failure(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_exonerate-rs"))
        .args(["--verbose", "0"])
        .args(args)
        .output()
        .expect("run exonerate-rs");
    assert!(
        !output.status.success(),
        "exonerate-rs unexpectedly succeeded"
    );
    String::from_utf8(output.stderr).expect("CLI stderr is UTF-8")
}

fn run_with_defaults(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_exonerate-rs"))
        .args(args)
        .output()
        .expect("run exonerate-rs");
    assert!(
        output.status.success(),
        "exonerate-rs failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("CLI output is UTF-8")
}

fn test_directory(name: &str) -> PathBuf {
    let directory = std::env::temp_dir().join(format!(
        "exonerate-rs-{name}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock after Unix epoch")
            .as_nanos()
    ));
    fs::create_dir(&directory).expect("create test directory");
    directory
}

#[test]
fn help_and_version_aliases_are_available_without_fasta_inputs() {
    for argument in ["-h", "--shorthelp", "--help"] {
        let output = run(&[argument]);
        assert!(output.starts_with("Usage: exonerate-rs"));
    }
    for argument in ["-v", "--version"] {
        let output = run(&[argument]);
        assert_eq!(
            output,
            format!("exonerate-rs {}\n", env!("CARGO_PKG_VERSION"))
        );
    }
}

#[test]
fn explicit_verbose_header_and_footer_match_upstream_layout() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let output = run(&[
        "--verbose",
        "1",
        "--model",
        "ungapped",
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showvulgar",
        "no",
        "--forwardonly",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    let mut lines = output.lines();
    let command = lines.next().expect("command-line header");
    assert!(command.starts_with("Command line: ["));
    assert!(command.contains("--verbose 1 --model ungapped"));
    assert!(
        lines
            .next()
            .is_some_and(|line| line.starts_with("Hostname: ["))
    );
    assert_eq!(lines.next(), Some("-- completed exonerate analysis"));
    assert_eq!(lines.next(), None);
}

#[test]
fn default_and_signed_verbose_levels_match_upstream_behavior() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let common = [
        "--model",
        "ungapped",
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showvulgar",
        "no",
        "--forwardonly",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ];
    let output = run_with_defaults(&common);
    let expected_hostname = Command::new("hostname")
        .output()
        .expect("hostname command")
        .stdout;
    let expected_hostname = String::from_utf8(expected_hostname)
        .expect("hostname is UTF-8")
        .trim()
        .to_owned();
    assert!(output.starts_with("Command line: ["));
    assert!(output.contains(&format!("Hostname: [{expected_hostname}]\n")));
    assert!(output.ends_with("-- completed exonerate analysis\n"));

    let mut quiet = vec!["--verbose", "-1"];
    quiet.extend(common);
    assert_eq!(run_with_defaults(&quiet), "");
}

#[test]
fn exhaustive_ungapped_reports_and_score_boundary_match_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let common = [
        "--model",
        "ungapped",
        "--exhaustive",
        "yes",
        "--revcomp",
        "no",
        "--subopt",
        "no",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
    ];
    let mut passing = common.to_vec();
    passing.extend([
        "--score",
        "40",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        run(&passing),
        concat!(
            "sugar: q 0 8 + t 2 10 + 40\n",
            "cigar: q 0 8 + t 2 10 + 40  M 8\n",
            "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
        )
    );

    let mut filtered = common.to_vec();
    filtered.extend([
        "--score",
        "41",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&filtered), "");
}

#[test]
fn ungapped_bestn_retains_tied_dna_and_protein_hsps() {
    for (query_name, target_name, query_type, score, expected) in [
        (
            "tie-dna-query.fa",
            "tie-dna-target.fa",
            "dna",
            "40",
            concat!(
                "tie-query tie-target 40 0 8 0 8 1\n",
                "tie-query tie-target 40 0 8 12 20 2\n",
            ),
        ),
        (
            "tie-protein-query.fa",
            "tie-protein-target.fa",
            "protein",
            "44",
            concat!(
                "tie-protein-query tie-protein-target 44 0 4 0 4 1\n",
                "tie-protein-query tie-protein-target 44 0 4 8 12 2\n",
            ),
        ),
    ] {
        let query = fixture(query_name);
        let target = fixture(target_name);
        assert_eq!(
            run(&[
                "--model",
                "ungapped",
                "--querytype",
                query_type,
                "--targettype",
                query_type,
                "--exhaustive",
                "yes",
                "--subopt",
                "yes",
                "--revcomp",
                "no",
                "--score",
                score,
                "--bestn",
                "1",
                "--showalignment",
                "no",
                "--showsugar",
                "no",
                "--showcigar",
                "no",
                "--showvulgar",
                "no",
                "--ryo",
                "%qi %ti %s %qab %qae %tab %tae %r\\n",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected ungapped {query_type} tie order"
        );
    }
}

#[test]
fn protein_to_dna_ungapped_bestn_retains_tied_hsps() {
    let query = fixture("tie-protein-query.fa");
    let target = fixture("tie-protein-dna-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "ungapped",
            "--querytype",
            "protein",
            "--targettype",
            "dna",
            "--exhaustive",
            "yes",
            "--subopt",
            "yes",
            "--revcomp",
            "no",
            "--score",
            "44",
            "--bestn",
            "1",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %qab %qae %tab %tae %r\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "tie-protein-query tie-protein-dna-target 44 0 4 0 12 1\n",
            "tie-protein-query tie-protein-dna-target 44 0 4 15 27 2\n",
        )
    );
}

#[test]
fn bestn_multi_record_order_matches_the_upstream_oracle() {
    let query = fixture("multi-order-query.fa");
    let target = fixture("multi-order-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--exhaustive",
            "yes",
            "--revcomp",
            "no",
            "--score",
            "10",
            "--bestn",
            "1",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %r\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "a-second a-second-target 20 1\n",
            "z-first z-first-target 40 1\n",
        )
    );
}

#[test]
fn protein2genome_bestn_multi_record_order_matches_the_upstream_oracle() {
    let query = fixture("multi-protein-query.fa");
    let target = fixture("multi-protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--exhaustive",
            "yes",
            "--revcomp",
            "no",
            "--minintron",
            "30",
            "--score",
            "20",
            "--bestn",
            "1",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %r\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!("a-protein a-genome 44 1\n", "z-protein z-genome 125 1\n",)
    );
}

#[test]
fn composed_multi_record_order_and_score_boundaries_match_the_upstream_oracle() {
    let query = fixture("multi-coding-query.fa");
    let target = fixture("multi-coding-genome.fa");
    for (model, at_boundary, above_boundary, filtered) in [
        (
            "coding2genome",
            concat!(
                "z-coding z-coding-genome 20 0\n",
                "a-coding a-coding-genome 44 0\n",
            ),
            "21",
            "a-coding a-coding-genome 44 0\n",
        ),
        (
            "cdna2genome",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "z-coding a-coding-genome 28 0\n",
                "a-coding z-coding-genome 28 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
            "29",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
        ),
        (
            "genome2genome",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "z-coding a-coding-genome 28 0\n",
                "a-coding z-coding-genome 28 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
            "29",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
        ),
    ] {
        let common = [
            "--model",
            model,
            "--exhaustive",
            "yes",
            "--revcomp",
            "no",
            "--subopt",
            "no",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %r\\n",
        ];
        let mut boundary_arguments = common.to_vec();
        boundary_arguments.extend([
            "--score",
            "20",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(
            run(&boundary_arguments),
            at_boundary,
            "unexpected {model} multi-record order"
        );

        let mut filtered_arguments = common.to_vec();
        filtered_arguments.extend([
            "--score",
            above_boundary,
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(
            run(&filtered_arguments),
            filtered,
            "unexpected {model} score boundary"
        );
    }
}

#[test]
fn est2genome_multi_record_order_and_score_boundary_match_the_upstream_oracle() {
    let query = fixture("multi-coding-query.fa");
    let target = fixture("multi-coding-genome.fa");
    let common = [
        "--model",
        "est2genome",
        "--exhaustive",
        "yes",
        "--revcomp",
        "no",
        "--subopt",
        "no",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "no",
        "--ryo",
        "%qi %ti %s %r\\n",
    ];
    let mut boundary_arguments = common.to_vec();
    boundary_arguments.extend([
        "--score",
        "20",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        run(&boundary_arguments),
        concat!(
            "z-coding z-coding-genome 60 0\n",
            "z-coding a-coding-genome 28 0\n",
            "a-coding z-coding-genome 28 0\n",
            "a-coding a-coding-genome 60 0\n",
        )
    );

    let mut filtered_arguments = common.to_vec();
    filtered_arguments.extend([
        "--score",
        "29",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        run(&filtered_arguments),
        concat!(
            "z-coding z-coding-genome 60 0\n",
            "a-coding a-coding-genome 60 0\n",
        )
    );
}

#[test]
fn ner_and_ungapped_trans_multi_record_order_match_the_upstream_oracle() {
    let query = fixture("multi-coding-query.fa");
    let target = fixture("multi-coding-genome.fa");
    for (model, boundary, at_boundary, above_boundary, filtered) in [
        (
            "ner",
            "20",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "z-coding a-coding-genome 28 0\n",
                "a-coding z-coding-genome 28 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
            "29",
            concat!(
                "z-coding z-coding-genome 60 0\n",
                "a-coding a-coding-genome 60 0\n",
            ),
        ),
        (
            "ungapped:trans",
            "20",
            concat!(
                "z-coding z-coding-genome 20 0\n",
                "a-coding a-coding-genome 44 0\n",
            ),
            "21",
            "a-coding a-coding-genome 44 0\n",
        ),
    ] {
        let common = [
            "--model",
            model,
            "--exhaustive",
            "yes",
            "--revcomp",
            "no",
            "--subopt",
            "no",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %r\\n",
        ];
        let mut boundary_arguments = common.to_vec();
        boundary_arguments.extend([
            "--score",
            boundary,
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(
            run(&boundary_arguments),
            at_boundary,
            "unexpected {model} multi-record order"
        );

        let mut filtered_arguments = common.to_vec();
        filtered_arguments.extend([
            "--score",
            above_boundary,
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(
            run(&filtered_arguments),
            filtered,
            "unexpected {model} score boundary"
        );
    }
}

#[test]
fn fasta_byte_chunks_match_the_upstream_record_boundaries() {
    let query = fixture("multi-query.fa");
    let target = fixture("multi-target.fa");
    let common = [
        "--model",
        "affine:local",
        "--exhaustive",
        "yes",
        "--revcomp",
        "no",
        "--score",
        "10",
        "--bestn",
        "1",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
    ];
    for (chunk_option, expected) in [
        ("--querychunkid", "vulgar: q2 0 4 + t2 4 8 + 20 M 4 4\n"),
        (
            "--targetchunkid",
            concat!(
                "vulgar: q1 1 3 + t2 3 5 + 10 M 2 2\n",
                "vulgar: q1 5 7 + t2 3 5 + 10 M 2 2\n",
                "vulgar: q2 0 4 + t2 4 8 + 20 M 4 4\n",
            ),
        ),
    ] {
        let mut arguments = common.to_vec();
        arguments.extend([
            chunk_option,
            "2",
            if chunk_option == "--querychunkid" {
                "--querychunktotal"
            } else {
                "--targetchunktotal"
            },
            "2",
        ]);
        arguments.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
        assert_eq!(
            run(&arguments),
            expected,
            "unexpected {chunk_option} output"
        );
    }

    let invalid = run_failure(&[
        "--querychunkid",
        "3",
        "--querychunktotal",
        "2",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(invalid.contains("chunk id should be between 1 and 2"));
}

#[test]
fn exhaustive_bestn_retains_score_ties_like_the_upstream_oracle() {
    let query = fixture("multi-query.fa");
    let target = fixture("multi-target.fa");
    let output = run(&[
        "--model",
        "affine:local",
        "--exhaustive",
        "--subopt",
        "yes",
        "--score",
        "30",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--bestn",
        "1",
        "--ryo",
        "%qi %ti %s %r\\n",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: q1 0 8 + t1 2 10 + 40\n",
            "cigar: q1 0 8 + t1 2 10 + 40  M 8\n",
            "vulgar: q1 0 8 + t1 2 10 + 40 M 8 8\n",
            "q1 t1 40 1\n",
            "sugar: q1 8 0 - t1 2 10 + 40\n",
            "cigar: q1 8 0 - t1 2 10 + 40  M 8\n",
            "vulgar: q1 8 0 - t1 2 10 + 40 M 8 8\n",
            "q1 t1 40 2\n",
        )
    );
}

#[test]
fn percent_filter_matches_upstream_for_translated_and_composed_models() {
    let coding_query = fixture("coding-query.fa");
    let coding_target = fixture("coding-genome.fa");
    let protein_query = fixture("protein-query.fa");
    let protein_target = fixture("protein-genome.fa");
    let expected = [
        (
            "protein2genome",
            protein_query.as_path(),
            protein_target.as_path(),
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
        ),
        (
            "coding2genome",
            coding_query.as_path(),
            coding_target.as_path(),
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n",
        ),
        (
            "cdna2genome",
            coding_query.as_path(),
            coding_target.as_path(),
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        ),
        (
            "genome2genome",
            coding_query.as_path(),
            coding_target.as_path(),
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        ),
    ];
    for (model, query, target, output) in expected {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--exhaustive",
                "yes",
                "--revcomp",
                "no",
                "--minintron",
                "30",
                "--maxintron",
                "10000",
                "--percent",
                "50",
                "--showalignment",
                "no",
                "--showsugar",
                "no",
                "--showcigar",
                "no",
                "--showvulgar",
                "yes",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            output,
            "unexpected {model} percent-filtered output"
        );
    }
}

#[test]
fn ryo_definitions_preserve_full_fasta_headers() {
    let query = fixture("multi-query.fa");
    let target = fixture("multi-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--exhaustive",
            "--forwardonly",
            "--score",
            "30",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi|%qd|%ti|%td\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        "q1|q1 first query|t1|t1 first target\n"
    );
}

#[test]
fn ryo_equivalenced_and_percent_fields_match_upstream_oracles() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:global",
            "--exhaustive",
            "--forwardonly",
            "--gapopen",
            "-1",
            "--gapextend",
            "-1",
            "--score",
            "0",
            "--subopt",
            "no",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%et %ei %es %em %pc %pi %pI %ps %pS\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        "8 8 8 0 100.00 100.00 66.67 100.00 100.00\n"
    );

    let protein = fixture("protein-query.fa");
    let genome = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%et %ei %es %em %pc %pi %pI %ps %pS\\n",
            protein.to_str().unwrap(),
            genome.to_str().unwrap(),
        ]),
        "28 28 28 0 96.55 100.00 100.00 100.00 100.00\n"
    );
}

#[test]
fn ryo_coding_region_fields_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "coding2genome",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qcb|%qce|%qcl|%qcs\\n%tcb|%tce|%tcl|%tcs\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "0|114|117|AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTG\n",
            "CCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT\n",
            "\n",
            "0|148|117|AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAGAGGACAGCACTG\n",
            "CCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT\n",
            "\n",
        )
    );

    let protein = fixture("protein-query.fa");
    let genome = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%tcs",
            protein.to_str().unwrap(),
            genome.to_str().unwrap(),
        ]),
        concat!(
            "ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGGAGGCCTTCTCCCTCTTTGACAAGGATGGAG\n",
            "ATGGCACTATTACCACC\n",
            "AATAGTGCCATCTCCATCCTTGTCAAAGAGGGAGAAGGC\n",
        )
    );
}

#[test]
fn reverse_target_ryo_ranges_and_sequences_match_the_upstream_oracle() {
    let protein = fixture("protein-query.fa");
    let genome = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%tab,%tae,%tS|%tas|%tcs",
            protein.to_str().unwrap(),
            genome.to_str().unwrap(),
        ]),
        concat!(
            "0,133,+|ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n",
            "NNNNNNNNNNNNAGGGAGGCCTTCTCCCTCTTTGACAAGGATGGAGATGGCACTATTACCACC\n",
            "|ATGGCTGACCAGCTGACTGAGCAGATTGCAGAGTTCAAGGAGGCCTTCTCCCTCTTTGACAAGGATGGAG\n",
            "ATGGCACTATTACCACC\n",
            "127,88,-|AATAGTGCCATCTCCATCCTTGTCAAAGAGGGAGAAGGC\n",
            "|AATAGTGCCATCTCCATCCTTGTCAAAGAGGGAGAAGGC\n",
        )
    );
}

#[test]
fn invalid_cli_options_and_ryo_tokens_fail_explicitly() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let base = [query.to_str().unwrap(), target.to_str().unwrap()];

    let percent = run_failure(&["--percent", "101", base[0], base[1]]);
    assert!(percent.contains("percent threshold must be between 0 and 100"));

    let word_len = run_failure(&["--wordlen", "0", base[0], base[1]]);
    assert!(word_len.contains("word length must be between 1 and 31"));

    let invalid_type = run_failure(&["--querytype", "nope", base[0], base[1]]);
    assert!(invalid_type.contains("query type must be dna or protein"));

    let model_type_conflict = run_failure(&[
        "--model",
        "coding2genome",
        "--querytype",
        "protein",
        base[0],
        base[1],
    ]);
    assert!(model_type_conflict.contains("model coding2genome requires DNA query and target"));

    let repeated_model = run_failure(&[
        "--model",
        "affine:local",
        "--model",
        "affine:global",
        base[0],
        base[1],
    ]);
    assert!(repeated_model.contains("model was already specified"));

    let ryo = run_failure(&[
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "no",
        "--score",
        "0",
        "--ryo",
        "%x",
        base[0],
        base[1],
    ]);
    assert!(ryo.contains("unknown ryo token %x"));

    let coding = run_failure(&[
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "no",
        "--score",
        "0",
        "--ryo",
        "%qcs",
        base[0],
        base[1],
    ]);
    assert!(coding.contains("%qc requires a coding region"));
}

#[test]
fn pretty_affine_alignment_matches_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--exhaustive",
            "--forwardonly",
            "--score",
            "0",
            "--subopt",
            "no",
            "--showalignment",
            "yes",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "C4 Alignment:\n",
            "------------\n",
            "         Query: q\n",
            "        Target: t\n",
            "     Raw score: 40\n",
            "   Query range: 0 -> 8\n",
            "  Target range: 2 -> 10\n",
            "\n",
            "ACGTACGT\n",
            "||||||||\n",
            "ACGTACGT\n",
        )
    );
}

#[test]
fn pretty_composed_alignment_identifies_the_upstream_model_name() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let output = run(&[
        "--model",
        "coding2genome",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        "--showalignment",
        "yes",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "no",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(output.starts_with(concat!(
        "C4 Alignment:\n",
        "------------\n",
        "         Query: coding\n",
        "        Target: genome\n",
        "         Model: coding2genome\n",
        "     Raw score: 194\n",
    )));
}

#[test]
fn exhaustive_affine_subopt_order_matches_the_upstream_oracle() {
    let query = fixture("multi-query.fa");
    let target = fixture("multi-target.fa");
    let output = run(&[
        "--model",
        "affine:local",
        "--exhaustive",
        "--subopt",
        "yes",
        "--score",
        "15",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "vulgar: q1 0 8 + t1 2 10 + 40 M 8 8\n",
            "vulgar: q1 3 8 + t1 1 6 + 25 M 5 5\n",
            "vulgar: q1 0 4 + t1 6 10 + 20 M 4 4\n",
            "vulgar: q1 8 0 - t1 2 10 + 40 M 8 8\n",
            "vulgar: q1 5 0 - t1 1 6 + 25 M 5 5\n",
            "vulgar: q1 8 4 - t1 6 10 + 20 M 4 4\n",
            "vulgar: q2 0 4 + t2 4 8 + 20 M 4 4\n",
            "vulgar: q2 1 4 + t2 4 7 + 15 M 3 3\n",
            "vulgar: q2 0 3 + t2 5 8 + 15 M 3 3\n",
            "vulgar: q2 4 0 - t2 0 4 + 20 M 4 4\n",
            "vulgar: q2 3 0 - t2 0 3 + 15 M 3 3\n",
            "vulgar: q2 4 1 - t2 1 4 + 15 M 3 3\n",
        )
    );
}

#[test]
fn subopt_is_enabled_by_default_like_the_upstream_cli() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--exhaustive",
            "yes",
            "--revcomp",
            "no",
            "--score",
            "10",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "yes",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
            "vulgar: q 3 8 + t 1 6 + 25 M 5 5\n",
            "vulgar: q 0 4 + t 6 10 + 20 M 4 4\n",
        )
    );
}

#[test]
fn exhaustive_percent_self_filter_matches_the_upstream_oracle() {
    let query = fixture("multi-query.fa");
    let target = fixture("multi-target.fa");
    let output = run(&[
        "--model",
        "affine:local",
        "--exhaustive",
        "--subopt",
        "yes",
        "--score",
        "15",
        "--percent",
        "100",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "vulgar: q1 0 8 + t1 2 10 + 40 M 8 8\n",
            "vulgar: q1 8 0 - t1 2 10 + 40 M 8 8\n",
            "vulgar: q2 0 4 + t2 4 8 + 20 M 4 4\n",
            "vulgar: q2 4 0 - t2 0 4 + 20 M 4 4\n",
        )
    );
}

#[test]
fn protein2genome_splice_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    let output = run(&[
        "--model",
        "protein2genome",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: protein 0 29 . genome 0 133 + 125\n",
            "cigar: protein 0 29 . genome 0 133 + 125  M 36 D 48 M 49\n",
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
        )
    );
}

#[test]
fn protein2genome_reverse_target_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    let output = run(&[
        "--model",
        "protein2genome",
        "--exhaustive",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: protein 0 29 . genome 0 133 + 125\n",
            "cigar: protein 0 29 . genome 0 133 + 125  M 36 D 48 M 49\n",
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
            "sugar: protein 12 25 . genome 127 88 - 28\n",
            "cigar: protein 12 25 . genome 127 88 - 28  M 39\n",
            "vulgar: protein 12 25 . genome 127 88 - 28 M 13 39\n",
        )
    );
}

#[test]
fn forwardcoordinates_no_reports_reverse_oriented_coordinates() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--forwardcoordinates",
            "no",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: protein 0 29 . genome 0 133 + 125\n",
            "cigar: protein 0 29 . genome 0 133 + 125  M 36 D 48 M 49\n",
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
            "sugar: protein 12 25 . genome 6 45 - 28\n",
            "cigar: protein 12 25 . genome 6 45 - 28  M 39\n",
            "vulgar: protein 12 25 . genome 6 45 - 28 M 13 39\n",
        )
    );
}

#[test]
fn protein2genome_bestfit_splice_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    let output = run(&[
        "--model",
        "protein2genome:bestfit",
        "--exhaustive",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: protein 0 29 . genome 0 133 + 125\n",
            "cigar: protein 0 29 . genome 0 133 + 125  M 36 D 48 M 49\n",
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
        )
    );
}

#[test]
fn exhaustive_protein2genome_zero_dpmemory_matches_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    let output = run(&[
        "--model",
        "protein2genome",
        "--exhaustive",
        "--dpmemory",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: protein 0 29 . genome 0 133 + 125\n",
            "cigar: protein 0 29 . genome 0 133 + 125  M 36 D 48 M 49\n",
            "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n",
        )
    );
}

#[test]
fn coding2genome_split_codon_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let output = run(&[
        "--model",
        "coding2genome",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 194\n",
            "cigar: coding 0 117 + genome 0 151 + 194  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n",
        )
    );
}

#[test]
fn coding_suboptimal_order_matches_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--exhaustive",
        "--forwardonly",
        "--subopt",
        "yes",
        "--score",
        "20",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];

    let mut coding2coding = vec!["--model", "coding2coding"];
    coding2coding.extend(common);
    coding2coding.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    assert_eq!(
        run(&coding2coding),
        concat!(
            "vulgar: coding 30 117 + genome 64 151 + 155 C 87 87\n",
            "vulgar: coding 0 27 + genome 0 27 + 51 C 27 27\n",
            "vulgar: coding 27 54 + genome 112 139 + 20 C 27 27\n",
        )
    );

    let mut coding2genome = vec!["--model", "coding2genome"];
    coding2genome.extend(common);
    coding2genome.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    assert_eq!(
        run(&coding2genome),
        concat!(
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n",
            "vulgar: coding 27 54 + genome 112 139 + 20 C 27 27\n",
        )
    );
}

#[test]
fn coding2genome_all_strand_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "coding2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 194\n",
            "cigar: coding 0 117 + genome 0 151 + 194  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n",
            "sugar: coding 9 48 + genome 125 86 - 22\n",
            "cigar: coding 9 48 + genome 125 86 - 22  M 39\n",
            "vulgar: coding 9 48 + genome 125 86 - 22 C 39 39\n",
            "sugar: coding 62 44 - genome 66 84 + 21\n",
            "cigar: coding 62 44 - genome 66 84 + 21  M 18\n",
            "vulgar: coding 62 44 - genome 66 84 + 21 C 18 18\n",
            "sugar: coding 115 1 - genome 149 1 - 185\n",
            "cigar: coding 115 1 - genome 149 1 - 185  M 89 D 34 M 25\n",
            "vulgar: coding 115 1 - genome 149 1 - 185 C 87 87 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 24 24\n",
        )
    );
}

#[test]
fn coding2genome_all_strand_zero_dpmemory_matches_full_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "coding2genome",
        "--exhaustive",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut checkpointed = common.to_vec();
    checkpointed.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&checkpointed), run(&full));
}

#[test]
fn exhaustive_coding2genome_zero_dpmemory_matches_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let output = run(&[
        "--model",
        "coding2genome",
        "--exhaustive",
        "--dpmemory",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 194\n",
            "cigar: coding 0 117 + genome 0 151 + 194  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n",
        )
    );
}

#[test]
fn exhaustive_cdna2genome_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "cdna2genome",
            "--exhaustive",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        )
    );
}

#[test]
fn exhaustive_cdna2genome_zero_dpmemory_matches_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "cdna2genome",
            "--exhaustive",
            "--dpmemory",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        )
    );
}

#[test]
fn cdna2genome_all_strand_zero_dpmemory_matches_full_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "cdna2genome",
        "--exhaustive",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut checkpointed = common.to_vec();
    checkpointed.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&checkpointed), run(&full));
}

#[test]
fn cdna2genome_all_strand_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "cdna2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
            "sugar: coding 35 91 + genome 132 75 - 57\n",
            "cigar: coding 35 91 + genome 132 75 - 57  M 5 D 1 M 7 I 1 M 5 D 1 M 19 M 6 M 1 I 1 M 5 D 1 M 6\n",
            "vulgar: coding 35 91 + genome 132 75 - 57 M 5 5 G 0 1 M 7 7 G 1 0 M 5 5 G 0 1 M 19 19 C 6 6 M 1 1 G 1 0 M 5 5 G 0 1 M 6 6\n",
            "sugar: coding 98 22 - genome 69 147 + 58\n",
            "cigar: coding 98 22 - genome 69 147 + 58  M 5 I 1 M 7 D 1 M 5 I 1 M 26 D 1 M 5 I 1 M 6 D 3 M 2 M 12 M 5\n",
            "vulgar: coding 98 22 - genome 69 147 + 58 M 5 5 G 1 0 M 7 7 G 0 1 M 5 5 G 1 0 M 26 26 G 0 1 M 5 5 G 1 0 M 6 6 G 0 3 M 2 2 C 12 12 M 5 5\n",
            "sugar: coding 117 0 - genome 151 0 - 545\n",
            "cigar: coding 117 0 - genome 151 0 - 545  M 88 D 34 M 29\n",
            "vulgar: coding 117 0 - genome 151 0 - 545 M 88 88 5 0 2 I 0 30 3 0 2 M 29 29\n",
        )
    );
}

#[test]
fn exhaustive_genome2genome_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "genome2genome",
            "--exhaustive",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        )
    );
}

#[test]
fn genome2genome_all_strand_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "genome2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
            "sugar: coding 35 91 + genome 132 75 - 57\n",
            "cigar: coding 35 91 + genome 132 75 - 57  M 5 D 1 M 7 I 1 M 5 D 1 M 19 M 6 M 1 I 1 M 5 D 1 M 6\n",
            "vulgar: coding 35 91 + genome 132 75 - 57 M 5 5 G 0 1 M 7 7 G 1 0 M 5 5 G 0 1 M 19 19 C 6 6 M 1 1 G 1 0 M 5 5 G 0 1 M 6 6\n",
            "sugar: coding 98 22 - genome 69 147 + 58\n",
            "cigar: coding 98 22 - genome 69 147 + 58  M 5 I 1 M 7 D 1 M 5 I 1 M 26 D 1 M 5 I 1 M 6 D 3 M 2 M 12 M 5\n",
            "vulgar: coding 98 22 - genome 69 147 + 58 M 5 5 G 1 0 M 7 7 G 0 1 M 5 5 G 1 0 M 26 26 G 0 1 M 5 5 G 1 0 M 6 6 G 0 3 M 2 2 C 12 12 M 5 5\n",
            "sugar: coding 117 0 - genome 151 0 - 545\n",
            "cigar: coding 117 0 - genome 151 0 - 545  M 88 D 34 M 29\n",
            "vulgar: coding 117 0 - genome 151 0 - 545 M 88 88 5 0 2 I 0 30 3 0 2 M 29 29\n",
        )
    );
}

#[test]
fn heuristic_genome2genome_preserves_all_strand_oracle_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "genome2genome",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut heuristic = common.to_vec();
    heuristic.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut exhaustive = common.to_vec();
    exhaustive.splice(2..2, ["--exhaustive"]);
    exhaustive.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    assert_eq!(run(&heuristic), run(&exhaustive));
}

#[test]
fn heuristic_genome2genome_zero_dpmemory_matches_full_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "genome2genome",
        "--subopt",
        "no",
        "--score",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut rolling = common.to_vec();
    rolling.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&rolling), run(&full));
}

#[test]
fn genome2genome_all_strand_suboptimal_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "genome2genome",
            "--exhaustive",
            "--subopt",
            "yes",
            "--score",
            "500",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
            "sugar: coding 117 0 - genome 151 0 - 545\n",
            "cigar: coding 117 0 - genome 151 0 - 545  M 88 D 34 M 29\n",
            "vulgar: coding 117 0 - genome 151 0 - 545 M 88 88 5 0 2 I 0 30 3 0 2 M 29 29\n",
        )
    );
}

#[test]
fn exhaustive_genome2genome_zero_dpmemory_preserves_oracle_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "genome2genome",
            "--exhaustive",
            "--dpmemory",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        )
    );
}

#[test]
fn exhaustive_genome2genome_small_dpmemory_matches_full_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "genome2genome",
        "--exhaustive",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut low_memory = common.to_vec();
    low_memory.extend([
        "--dpmemory",
        "1",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&low_memory), run(&full));
}

#[test]
fn genome2genome_all_strand_zero_dpmemory_matches_full_output() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "genome2genome",
        "--exhaustive",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--minintron",
        "30",
        "--maxintron",
        "10000",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut rolling = common.to_vec();
    rolling.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&rolling), run(&full));
}

#[test]
fn exhaustive_ner_reports_match_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "ner",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: q 0 8 + t 2 10 + 40\n",
            "cigar: q 0 8 + t 2 10 + 40  M 8\n",
            "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
        )
    );
}

#[test]
fn ner_zero_dpmemory_matches_full_generic_output() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let common = [
        "--model",
        "ner",
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
    ];
    let mut full = common.to_vec();
    full.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let mut checkpointed = common.to_vec();
    checkpointed.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(run(&checkpointed), run(&full));
}

#[test]
fn exhaustive_ungapped_trans_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "ungapped:trans",
            "--exhaustive",
            "--subopt",
            "no",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 30 117 + genome 64 151 + 155\n",
            "cigar: coding 30 117 + genome 64 151 + 155  M 87\n",
            "vulgar: coding 30 117 + genome 64 151 + 155 C 87 87\n",
        )
    );
}

#[test]
fn ungapped_translated_zero_dpmemory_uses_linear_generic_backend() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let output = run(&[
        "--model",
        "ungapped:trans",
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--dpmemory",
        "0",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(output.contains("sugar: coding 30 117 + genome 64 151 + 155"));
}

#[test]
fn exhaustive_protein2dna_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2dna",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: protein 12 29 . genome 82 133 + 85\n",
            "cigar: protein 12 29 . genome 82 133 + 85  M 51\n",
            "vulgar: protein 12 29 . genome 82 133 + 85 M 17 51\n",
        )
    );
}

#[test]
fn exhaustive_protein2dna_reverse_target_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2dna",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: protein 12 29 . genome 82 133 + 85\n",
            "cigar: protein 12 29 . genome 82 133 + 85  M 51\n",
            "vulgar: protein 12 29 . genome 82 133 + 85 M 17 51\n",
            "sugar: protein 12 25 . genome 127 88 - 28\n",
            "cigar: protein 12 25 . genome 127 88 - 28  M 39\n",
            "vulgar: protein 12 25 . genome 127 88 - 28 M 13 39\n",
        )
    );
}

#[test]
fn exhaustive_protein2dna_bestfit_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2dna:bestfit",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: protein 0 29 . genome 46 133 + 76\n",
            "cigar: protein 0 29 . genome 46 133 + 76  M 87\n",
            "vulgar: protein 0 29 . genome 46 133 + 76 M 29 87\n",
        )
    );
}

#[test]
fn exhaustive_coding2coding_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "coding2coding",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 30 117 + genome 64 151 + 155\n",
            "cigar: coding 30 117 + genome 64 151 + 155  M 87\n",
            "vulgar: coding 30 117 + genome 64 151 + 155 C 87 87\n",
        )
    );
}

#[test]
fn exhaustive_affine_scopes_match_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    for (model, expected) in [
        (
            "affine:global",
            concat!(
                "sugar: q 0 8 + t 0 12 + 8\n",
                "cigar: q 0 8 + t 0 12 + 8 D 2 M 8 D 2\n",
                "vulgar: q 0 8 + t 0 12 + 8 G 0 2 M 8 8 G 0 2\n",
            ),
        ),
        (
            "affine:bestfit",
            concat!(
                "sugar: q 0 8 + t 2 10 + 40\n",
                "cigar: q 0 8 + t 2 10 + 40  M 8\n",
                "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
            ),
        ),
        (
            "affine:local",
            concat!(
                "sugar: q 0 8 + t 2 10 + 40\n",
                "cigar: q 0 8 + t 2 10 + 40  M 8\n",
                "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
            ),
        ),
        (
            "affine:overlap",
            concat!(
                "sugar: q 0 8 + t 2 10 + 40\n",
                "cigar: q 0 8 + t 2 10 + 40  M 8\n",
                "vulgar: q 0 8 + t 2 10 + 40 M 8 8\n",
            ),
        ),
    ] {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--exhaustive",
                "--subopt",
                "no",
                "--score",
                "0",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "yes",
                "--showvulgar",
                "yes",
                "--forwardonly",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected {model} output",
        );
    }
}

#[test]
fn affine_global_and_bestfit_subopt_match_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    for (model, expected) in [
        ("affine:global", "sugar: q 0 8 + t 0 12 + 8\n"),
        (
            "affine:bestfit",
            concat!(
                "sugar: q 0 8 + t 2 10 + 40\n",
                "sugar: q 0 8 + t 0 6 + 5\n",
                "sugar: q 0 8 + t 6 12 + 5\n",
            ),
        ),
    ] {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--exhaustive",
                "yes",
                "--subopt",
                "yes",
                "--score",
                "0",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "no",
                "--showvulgar",
                "no",
                "--revcomp",
                "no",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected {model} suboptimal output",
        );
    }
}

#[test]
fn affine_global_and_bestfit_bestn_match_the_upstream_oracle() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    for (model, expected) in [
        ("affine:global", "sugar: q 0 8 + t 0 12 + 8\n"),
        ("affine:bestfit", "sugar: q 0 8 + t 2 10 + 40\n"),
    ] {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--exhaustive",
                "yes",
                "--subopt",
                "yes",
                "--bestn",
                "1",
                "--score",
                "0",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "no",
                "--showvulgar",
                "no",
                "--revcomp",
                "no",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected {model} bestn output",
        );
    }
}

#[test]
fn protein_affine_reports_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-target.fa");
    let output = run(&[
        "--model",
        "affine:local",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(
        output,
        concat!(
            "sugar: protein 0 29 . protein-target 2 31 . 146\n",
            "cigar: protein 0 29 . protein-target 2 31 . 146  M 29\n",
            "vulgar: protein 0 29 . protein-target 2 31 . 146 M 29 29\n",
        )
    );
}

#[test]
fn pretty_protein_affine_identifies_the_upstream_model_name() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-target.fa");
    let output = run(&[
        "--model",
        "affine:local",
        "--showalignment",
        "yes",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "no",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(output.contains("         Model: affine:local:protein2protein\n"));
}

#[test]
fn exhaustive_protein_affine_subopt_matches_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--querytype",
            "protein",
            "--targettype",
            "protein",
            "--exhaustive",
            "yes",
            "--subopt",
            "yes",
            "--score",
            "20",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "yes",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        "vulgar: protein 0 29 . protein-target 2 31 . 146 M 29 29\n"
    );
}

#[test]
fn protein_affine_nonlocal_scopes_subopt_match_the_upstream_oracle() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-target.fa");
    for (model, expected) in [
        (
            "affine:global",
            "sugar: protein 0 29 . protein-target 0 33 . 114\n",
        ),
        (
            "affine:bestfit",
            "sugar: protein 0 29 . protein-target 2 31 . 146\n",
        ),
        (
            "affine:overlap",
            "sugar: protein 0 29 . protein-target 2 31 . 146\n",
        ),
    ] {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--querytype",
                "protein",
                "--targettype",
                "protein",
                "--exhaustive",
                "yes",
                "--subopt",
                "yes",
                "--score",
                "20",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "no",
                "--showvulgar",
                "no",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected {model} protein suboptimal output",
        );
    }
}

#[test]
fn protein_affine_bestn_retains_tied_suboptimal_loci() {
    let query = fixture("tie-protein-query.fa");
    let target = fixture("tie-protein-target.fa");
    assert_eq!(
        run(&[
            "--model",
            "affine:local",
            "--querytype",
            "protein",
            "--targettype",
            "protein",
            "--exhaustive",
            "yes",
            "--subopt",
            "yes",
            "--score",
            "44",
            "--bestn",
            "1",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "%qi %ti %s %qab %qae %tab %tae %r\\n",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "tie-protein-query tie-protein-target 44 0 4 0 4 1\n",
            "tie-protein-query tie-protein-target 44 0 4 8 12 2\n",
        )
    );
}

#[test]
fn zero_dp_memory_preserves_affine_cli_reports() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let common = [
        "--model",
        "affine:global",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        "--forwardonly",
    ];
    let mut full_args = common.to_vec();
    full_args.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
    let full = run(&full_args);

    let mut checkpointed_args = common.to_vec();
    checkpointed_args.extend([
        "--dpmemory",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    let checkpointed = run(&checkpointed_args);
    assert_eq!(checkpointed, full);
}

#[test]
fn zero_dp_memory_preserves_est2genome_cli_reports() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let common = [
        "--model",
        "est2genome",
        "--exhaustive",
        "--forwardonly",
        "--minintron",
        "2",
        "--maxintron",
        "10000",
        "--showalignment",
        "no",
        "--showsugar",
        "yes",
        "--showcigar",
        "yes",
        "--showvulgar",
        "yes",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ];
    let full = run(&common);
    let mut checkpointed = common.to_vec();
    checkpointed.splice(3..3, ["--dpmemory", "0"]);
    assert_eq!(run(&checkpointed), full);
}

#[test]
fn zero_dp_memory_preserves_raw_transition_reports() {
    for (model, query_name, target_name, extra) in [
        (
            "ner",
            "checkpoint-zero-query.fa",
            "checkpoint-zero-target.fa",
            Vec::<&str>::new(),
        ),
        (
            "genome2genome",
            "checkpoint-splice-query.fa",
            "checkpoint-splice-target.fa",
            vec![
                "--minintron",
                "6",
                "--maxintron",
                "20",
                "--intronpenalty",
                "-1",
                "--forcegtag",
                "yes",
            ],
        ),
    ] {
        let query = fixture(query_name);
        let target = fixture(target_name);
        let mut common = vec![
            "--model",
            model,
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--forwardonly",
            "--score",
            "0",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "{%Pn|%Ps|%Pqa|%Pta\\n}",
        ];
        common.extend(extra);
        let mut full = common.clone();
        full.extend([
            "--dpmemory",
            "4096",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        let mut checkpointed = common;
        checkpointed.extend([
            "--dpmemory",
            "0",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(
            run(&checkpointed),
            run(&full),
            "raw transition mismatch for {model}"
        );
    }
}

#[test]
fn exhaustive_est2genome_reports_match_the_upstream_oracle() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "est2genome",
            "--exhaustive",
            "--subopt",
            "no",
            "--score",
            "0",
            "--forwardonly",
            "--minintron",
            "30",
            "--maxintron",
            "10000",
            "--showalignment",
            "no",
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]),
        concat!(
            "sugar: coding 0 117 + genome 0 151 + 564\n",
            "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
            "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
        )
    );
}

#[test]
fn exhaustive_accepts_upstream_boolean_syntax() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let common = [
        "--model",
        "affine:local",
        "--showalignment",
        "no",
        "--showsugar",
        "no",
        "--showcigar",
        "no",
        "--showvulgar",
        "yes",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ];
    let mut short_form = common.to_vec();
    short_form.splice(2..2, ["-E"]);
    let mut upstream_form = common.to_vec();
    upstream_form.splice(2..2, ["--exhaustive", "yes"]);
    assert_eq!(run(&upstream_form), run(&short_form));

    let mut disabled_form = common.to_vec();
    disabled_form.splice(2..2, ["--exhaustive", "no"]);
    assert_eq!(run(&disabled_form), run(&common));
}

#[test]
fn short_minimum_intron_is_safe_and_matches_complex_model_oracles() {
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let expected = concat!(
        "sugar: coding 0 117 + genome 0 151 + 564\n",
        "cigar: coding 0 117 + genome 0 151 + 564  M 29 D 34 M 88\n",
        "vulgar: coding 0 117 + genome 0 151 + 564 M 29 29 5 0 2 I 0 30 3 0 2 M 88 88\n",
    );
    for model in ["est2genome", "cdna2genome", "genome2genome"] {
        assert_eq!(
            run(&[
                "--model",
                model,
                "--exhaustive",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "yes",
                "--showvulgar",
                "yes",
                "--forwardonly",
                "--minintron",
                "2",
                "--maxintron",
                "10000",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]),
            expected,
            "unexpected {model} output"
        );
    }
}

#[test]
fn short_minimum_intron_is_safe_for_phase_intron_models() {
    let protein_query = fixture("protein-query.fa");
    let protein_target = fixture("protein-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--exhaustive",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "2",
            "--maxintron",
            "10000",
            protein_query.to_str().unwrap(),
            protein_target.to_str().unwrap(),
        ]),
        "vulgar: protein 0 29 . genome 0 133 + 125 M 12 36 S 0 2 5 0 2 I 0 42 3 0 2 S 1 1 M 16 48\n"
    );

    let coding_query = fixture("coding-query.fa");
    let coding_target = fixture("coding-genome.fa");
    assert_eq!(
        run(&[
            "--model",
            "coding2genome",
            "--exhaustive",
            "--score",
            "190",
            "--showalignment",
            "no",
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "yes",
            "--forwardonly",
            "--minintron",
            "2",
            "--maxintron",
            "10000",
            coding_query.to_str().unwrap(),
            coding_target.to_str().unwrap(),
        ]),
        "vulgar: coding 0 117 + genome 0 151 + 194 C 27 27 S 2 2 5 0 2 I 0 30 3 0 2 S 1 1 C 87 87\n"
    );
}

#[test]
fn machine_readable_reports_cover_single_batch_and_audit_workflows() {
    let directory = test_directory("machine-reports");
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    let single_tsv = directory.join("single.tsv");
    let single_gff3 = directory.join("single.gff3");
    let output = run(&[
        "--model",
        "affine:local",
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--score",
        "0",
        "--forwardonly",
        "--showalignment",
        "no",
        "--showvulgar",
        "no",
        "--result-tsv",
        single_tsv.to_str().unwrap(),
        "--evidence-gff3",
        single_gff3.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert_eq!(output, "");
    let single_rows = fs::read_to_string(&single_tsv).expect("read single TSV");
    assert!(single_rows.starts_with("task_id\tstatus\tmodel\trank\t"));
    assert!(single_rows.contains("\taligned\taffine:local:dna2dna\t"));
    assert!(single_rows.contains("query_aligned"));
    let single_evidence = fs::read_to_string(&single_gff3).expect("read single GFF3");
    assert!(single_evidence.starts_with("##gff-version 3\n"));
    assert!(single_evidence.contains("\tmatch\t"));
    assert!(!single_evidence.contains("Task=;"));

    let manifest = directory.join("tasks.tsv");
    fs::write(
        &manifest,
        format!(
            "task_id\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\nfirst\taffine:local\t{}\tq\t{}\tt\nsecond\taffine:local\t{}\tq\t{}\tt\n",
            query.display(), target.display(), query.display(), target.display()
        ),
    )
    .expect("write task manifest");
    let batch_tsv = directory.join("batch.tsv");
    let batch_gff3 = directory.join("batch.gff3");
    assert_eq!(
        run(&[
            "--tasks",
            manifest.to_str().unwrap(),
            "--result-tsv",
            batch_tsv.to_str().unwrap(),
            "--evidence-gff3",
            batch_gff3.to_str().unwrap(),
            "--threads",
            "2",
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--score",
            "0",
            "--forwardonly",
        ]),
        ""
    );
    let batch_rows = fs::read_to_string(&batch_tsv).expect("read batch TSV");
    let rows = batch_rows.lines().collect::<Vec<_>>();
    assert_eq!(rows.len(), 3);
    assert!(rows[1].starts_with("first\taligned\t"));
    assert!(rows[2].starts_with("second\taligned\t"));
    let batch_evidence = fs::read_to_string(&batch_gff3).expect("read batch GFF3");
    assert!(batch_evidence.contains("Task=first"));
    assert!(batch_evidence.contains("Task=second"));
    assert!(batch_evidence.contains("ID=alignment.5:first.1:t.1"));
    assert!(batch_evidence.contains("ID=alignment.6:second.1:t.1"));

    let protein = fixture("protein-query.fa");
    let genome = fixture("protein-genome.fa");
    let audit_tsv = directory.join("audit.tsv");
    assert_eq!(
        run(&[
            "--model",
            "protein2genome",
            "--audit",
            "protein-candidate",
            "--result-tsv",
            audit_tsv.to_str().unwrap(),
            "--forwardonly",
            protein.to_str().unwrap(),
            genome.to_str().unwrap(),
        ]),
        ""
    );
    assert!(
        fs::read_to_string(&audit_tsv)
            .expect("read audit TSV")
            .contains("\taligned\tprotein2genome:local\t")
    );

    let competitor = directory.join("competitor.fa");
    fs::write(&competitor, ">Fam_014\nWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n")
        .expect("write competing protein");
    let orthology_manifest = directory.join("orthology-tasks.tsv");
    fs::write(
        &orthology_manifest,
        format!(
            "task_id\tsample_id\tlocus_id\tcandidate_id\texpected_family\treference_family\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\nexpected\tS1\tLocus_001\tgenome\tFam_001\tFam_001\tprotein2genome\t{}\tprotein\t{}\tgenome\ncompetitor\tS1\tLocus_001\tgenome\tFam_001\tFam_014\tprotein2genome\t{}\tFam_014\t{}\tgenome\n",
            protein.display(),
            genome.display(),
            competitor.display(),
            genome.display(),
        ),
    )
    .expect("write orthology task manifest");
    let orthology_tsv = directory.join("orthology.tsv");
    assert_eq!(
        run(&[
            "--tasks",
            orthology_manifest.to_str().unwrap(),
            "--audit",
            "protein-candidate",
            "--result-tsv",
            directory.join("orthology-alignments.tsv").to_str().unwrap(),
            "--orthology-report",
            orthology_tsv.to_str().unwrap(),
            "--orthology-min-coverage",
            "0.8",
            "--orthology-min-delta",
            "20",
            "--forwardonly",
        ]),
        ""
    );
    let orthology_rows = fs::read_to_string(&orthology_tsv).expect("read orthology report");
    let orthology_rows = orthology_rows.lines().collect::<Vec<_>>();
    assert_eq!(
        orthology_rows[0],
        "sample_id\tlocus_id\tcandidate_id\texpected_family\tstatus\texpected_score\texpected_coverage\tbest_competitor_family\tbest_competitor_score\tscore_delta\tframeshift_bases\tintron_count\texpected_task_id"
    );
    let values = orthology_rows[1].split('\t').collect::<Vec<_>>();
    assert_eq!(values.len(), 13);
    assert_eq!(
        values[0..5],
        ["S1", "Locus_001", "genome", "Fam_001", "accepted"]
    );
    assert_eq!(values[12], "expected");
    fs::remove_dir_all(directory).expect("remove test directory");
}

#[test]
fn orthology_uses_reference_coverage_and_the_winning_task() {
    let directory = test_directory("orthology-reference-coverage");
    let weak = directory.join("weak.fa");
    fs::write(&weak, ">weak\nM\n").unwrap();
    let genome = directory.join("flanked.fa");
    let original = fs::read_to_string(fixture("protein-genome.fa")).unwrap();
    let bases = original.lines().skip(1).collect::<String>();
    fs::write(
        &genome,
        format!(">genome\n{}{}{}\n", "N".repeat(100), bases, "N".repeat(100)),
    )
    .unwrap();
    let manifest = directory.join("tasks.tsv");
    fs::write(&manifest, format!(
        "task_id\tsample_id\tlocus_id\tcandidate_id\texpected_family\treference_family\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\nweak\tS\tL\tC\tF\tF\tprotein2genome\t{}\tweak\t{}\tgenome\nstrong\tS\tL\tC\tF\tF\tprotein2genome\t{}\tprotein\t{}\tgenome\n",
        weak.display(), genome.display(), fixture("protein-query.fa").display(), genome.display()
    )).unwrap();
    let evidence = directory.join("evidence.tsv");
    let report = directory.join("orthology.tsv");
    run(&[
        "--tasks",
        manifest.to_str().unwrap(),
        "--audit",
        "protein-candidate",
        "--result-tsv",
        evidence.to_str().unwrap(),
        "--orthology-report",
        report.to_str().unwrap(),
        "--forwardonly",
    ]);
    let report_text = fs::read_to_string(report).unwrap();
    let row = report_text
        .lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(row[4], "accepted");
    assert_eq!(row[6], "1.000000");
    assert_eq!(row[12], "strong");
    let evidence_text = fs::read_to_string(evidence).unwrap();
    let best = evidence_text
        .lines()
        .find(|line| line.starts_with("strong\t"))
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(row[5], best[12]);
    assert!(best[18].parse::<f64>().unwrap() < 0.8);
    assert_eq!(best[20], "0", "an intron is not an ordinary gap");
    assert_eq!(best[22], "1");
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn batch_failures_preserve_tsv_columns_and_order_across_thread_counts() {
    let directory = test_directory("batch-failures");
    let manifest = directory.join("tasks.tsv");
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    fs::write(&manifest, format!("task_id\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\nok\taffine:local\t{}\tq\t{}\tt\nbad\taffine:local\t{}\tmissing\t{}\tt\n", query.display(), target.display(), query.display(), target.display())).unwrap();
    let mut previous = None;
    for threads in ["1", "2"] {
        let result = directory.join(format!("{threads}.tsv"));
        run(&[
            "--tasks",
            manifest.to_str().unwrap(),
            "--result-tsv",
            result.to_str().unwrap(),
            "--threads",
            threads,
            "--subopt",
            "no",
            "--forwardonly",
            "--score",
            "0",
        ]);
        let text = fs::read_to_string(result).unwrap();
        let rows = text
            .lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|row| row.len() == 24));
        assert_eq!(&rows[1][..2], &["ok", "aligned"]);
        assert_eq!(&rows[2][..2], &["bad", "failed"]);
        assert!(rows[2][23].contains("was not found"));
        if let Some(previous) = &previous {
            assert_eq!(&text, previous);
        }
        previous = Some(text);
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn duplicate_task_ids_are_rejected_before_outputs_are_created() {
    let directory = test_directory("duplicate-task-ids");
    let manifest = directory.join("tasks.tsv");
    fs::write(&manifest, "task_id\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\nsame\tungapped\tx\tq\ty\tt\nsame\tungapped\tx\tq\ty\tt\n").unwrap();
    let output = directory.join("out.tsv");
    assert!(
        run_failure(&[
            "--tasks",
            manifest.to_str().unwrap(),
            "--result-tsv",
            output.to_str().unwrap()
        ])
        .contains("duplicate task ID")
    );
    assert!(!output.exists());
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn explicit_sequence_types_are_checked_independently_of_option_order() {
    let query = fixture("protein-query.fa");
    let target = fixture("protein-genome.fa");
    for options in [
        ["--querytype", "dna", "--model", "protein2dna"],
        ["--model", "protein2dna", "--querytype", "dna"],
        ["--targettype", "protein", "--model", "protein2dna"],
        ["--model", "protein2dna", "--targettype", "protein"],
    ] {
        let mut args = options.to_vec();
        args.extend([query.to_str().unwrap(), target.to_str().unwrap()]);
        assert!(run_failure(&args).contains("requires a protein query and DNA target"));
    }
}

#[test]
fn heuristic_tied_regions_are_deterministic_and_padding_cannot_overflow() {
    let directory = test_directory("heuristic-ties");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    fs::write(&query, ">q\nAAAAAAAAAAAA\n").unwrap();
    fs::write(
        &target,
        format!(
            ">t\n{}{}{}\n",
            "A".repeat(12),
            "G".repeat(100),
            "A".repeat(12)
        ),
    )
    .unwrap();
    let mut previous = None;
    for _ in 0..32 {
        let output = run(&[
            "--model",
            "affine:local",
            "--wordlen",
            "12",
            "--seedpadding",
            "0",
            "--subopt",
            "no",
            "--forwardonly",
            "--score",
            "0",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        if let Some(previous) = &previous {
            assert_eq!(&output, previous);
        }
        previous = Some(output);
    }
    run(&[
        "--model",
        "affine:local",
        "--seedpadding",
        &usize::MAX.to_string(),
        "--subopt",
        "no",
        "--forwardonly",
        "--score",
        "0",
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    fs::remove_dir_all(directory).unwrap();
}

#[cfg(target_os = "linux")]
#[test]
fn buffered_output_errors_are_reported_for_hits_and_no_hits() {
    let query = fixture("dna-query.fa");
    let target = fixture("dna-target.fa");
    for option in ["--result-tsv", "--evidence-gff3"] {
        for score in ["0", "2147483647"] {
            let error = run_failure(&[
                "--model",
                "affine:local",
                "--subopt",
                "no",
                "--score",
                score,
                "--forwardonly",
                option,
                "/dev/full",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]);
            assert!(error.contains("No space left on device"));
        }
    }
}

#[test]
fn protein2dna_local_alignment_can_start_in_each_target_reading_frame() {
    let directory = test_directory("protein-reading-frame-starts");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    fs::write(&query, ">q\nHYR\n").unwrap();
    fs::write(&target, ">t\nACGATAATG\n").unwrap();
    for subopt in ["no", "yes"] {
        let output = run(&[
            "--model",
            "protein2dna",
            "--exhaustive",
            "yes",
            "--subopt",
            subopt,
            "--score",
            "5",
            "--forwardonly",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(output, "vulgar: q 2 3 . t 1 4 + 5 M 1 3\n");
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn gff3_reverse_queries_stay_in_bounds_and_use_unique_ids() {
    let directory = test_directory("gff3-reverse-query");
    let query = fixture("coding-query.fa");
    let target = fixture("coding-genome.fa");
    let mut previous = None;
    for forward_coordinates in ["yes", "no"] {
        let evidence = directory.join(format!("{forward_coordinates}.gff3"));
        let output = run(&[
            "--model",
            "coding2genome",
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--score",
            "0",
            "--showvulgar",
            "no",
            "--showtargetgff",
            "yes",
            "--showquerygff",
            "yes",
            "--forwardcoordinates",
            forward_coordinates,
            "--evidence-gff3",
            evidence.to_str().unwrap(),
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        let evidence = fs::read_to_string(evidence).unwrap();
        for document in [&output, &evidence] {
            let mut ids = std::collections::HashSet::new();
            let mut reverse_query_seen = false;
            for line in document.lines().filter(|line| !line.starts_with('#')) {
                let columns = line.split('\t').collect::<Vec<_>>();
                assert_eq!(columns.len(), 9);
                let reference_length = if columns[0] == "coding" { 117 } else { 151 };
                let start = columns[3].parse::<usize>().unwrap();
                let end = columns[4].parse::<usize>().unwrap();
                assert!(
                    start >= 1 && start <= end && end <= reference_length,
                    "{line}"
                );
                let attributes = columns[8].split(';').collect::<Vec<_>>();
                let id = attributes
                    .iter()
                    .find_map(|attribute| attribute.strip_prefix("ID="))
                    .unwrap();
                assert!(ids.insert(id), "duplicate GFF3 ID {id}");
                let target = attributes
                    .iter()
                    .find_map(|attribute| attribute.strip_prefix("Target="))
                    .unwrap()
                    .split_whitespace()
                    .collect::<Vec<_>>();
                let target_length = if target[0] == "coding" { 117 } else { 151 };
                let start = target[1].parse::<usize>().unwrap();
                let end = target[2].parse::<usize>().unwrap();
                assert!(start >= 1 && start <= end && end <= target_length, "{line}");
                reverse_query_seen |= target[0] == "coding" && target.get(3) == Some(&"-");
            }
            assert!(reverse_query_seen);
        }
        if let Some((old_output, old_evidence)) = &previous {
            assert_eq!(&output, old_output);
            assert_eq!(&evidence, old_evidence);
        }
        previous = Some((output, evidence));
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn gff3_identifiers_encode_utf8_and_reserved_characters_once() {
    let directory = test_directory("gff3-utf8");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    let evidence = directory.join("out.gff3");
    fs::write(&query, ">查询%\nATGATGATG\n").unwrap();
    fs::write(&target, ">目;标\nATGATGATG\n").unwrap();
    let output = run(&[
        "--model",
        "ungapped",
        "--subopt",
        "no",
        "--score",
        "0",
        "--forwardonly",
        "--showvulgar",
        "no",
        "--showgff",
        "yes",
        "--evidence-gff3",
        evidence.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    let evidence = fs::read_to_string(evidence).unwrap();
    for document in [&output, &evidence] {
        assert!(document.is_ascii());
        assert!(document.contains("%E7%9B%AE%3B%E6%A0%87\texonerate-rs"));
        assert!(document.contains("Target=%E6%9F%A5%E8%AF%A2%25 1 9"));
        assert!(!document.contains("%25E7"));
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn report_paths_cannot_overwrite_inputs_or_each_other() {
    let directory = test_directory("output-path-collisions");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    let output = directory.join("out");
    fs::write(&query, ">q\nACGT\n").unwrap();
    fs::write(&target, ">t\nACGT\n").unwrap();
    let error = run_failure(&[
        "--result-tsv",
        query.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(error.contains("would overwrite an input"));
    assert_eq!(fs::read_to_string(&query).unwrap(), ">q\nACGT\n");
    let error = run_failure(&[
        "--result-tsv",
        output.to_str().unwrap(),
        "--evidence-gff3",
        output.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(error.contains("more than one report"));
    assert!(!output.exists());
    let manifest = directory.join("tasks.tsv");
    fs::write(&manifest, format!("task_id\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\none\tungapped\t{}\tq\t{}\tt\n", query.display(), target.display())).unwrap();
    let before = fs::read_to_string(&manifest).unwrap();
    let error = run_failure(&[
        "--tasks",
        manifest.to_str().unwrap(),
        "--result-tsv",
        manifest.to_str().unwrap(),
    ]);
    assert!(error.contains("would overwrite an input"));
    assert_eq!(fs::read_to_string(&manifest).unwrap(), before);
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn empty_local_alignments_do_not_emit_invalid_gff3_intervals() {
    let directory = test_directory("gff3-empty");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    let output = directory.join("out.gff3");
    fs::write(&query, ">q\nA\n").unwrap();
    fs::write(&target, ">t\nC\n").unwrap();
    let tsv = directory.join("out.tsv");
    let stdout = run(&[
        "--subopt",
        "no",
        "--score",
        "0",
        "--forwardonly",
        "--showvulgar",
        "no",
        "--showgff",
        "yes",
        "--showquerygff",
        "yes",
        "--result-tsv",
        tsv.to_str().unwrap(),
        "--evidence-gff3",
        output.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(stdout.is_empty());
    assert_eq!(
        fs::read_to_string(tsv)
            .unwrap()
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(1),
        Some("no_hit")
    );
    assert_eq!(fs::read_to_string(output).unwrap(), "##gff-version 3\n");
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn ner_span_hints_do_not_limit_exhaustive_paths_or_atomic_reports() {
    let directory = test_directory("ner-span-hints");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    fs::write(&query, ">q\nAGATACTTGGCTTTATCT\n").unwrap();
    fs::write(&target, ">t\nAGATAAAGCCAAGTATCT\n").unwrap();
    for (minimum, maximum) in [("10", "50000"), ("100", "100"), ("1", "1")] {
        for budget in ["32", "0"] {
            let output = run(&[
                "--model",
                "ner",
                "--exhaustive",
                "yes",
                "--subopt",
                "no",
                "--revcomp",
                "no",
                "--score",
                "1",
                "--showalignment",
                "no",
                "--minner",
                minimum,
                "--maxner",
                maximum,
                "--dpmemory",
                budget,
                "--ryo",
                "{%Pn|%Ps|%Pqa|%Pta\\n}",
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]);
            assert_eq!(output, include_str!("fixtures/ner_short_span_ryo.expected"));
        }
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn upstream_audit_counterexamples_match_in_both_memory_modes() {
    let directory = test_directory("upstream-audit");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    for (index, line) in include_str!("fixtures/upstream_audit.tsv")
        .lines()
        .enumerate()
    {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 7);
        fs::write(&query, format!(">q\n{}\n", fields[4])).unwrap();
        fs::write(&target, format!(">t\n{}\n", fields[5])).unwrap();
        for budget in ["32", "0"] {
            let output = run(&[
                "--model",
                fields[0],
                "--revcomp",
                fields[1],
                "--querytype",
                fields[2],
                "--targettype",
                fields[3],
                "--exhaustive",
                "yes",
                "--subopt",
                "no",
                "--score",
                "1",
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "yes",
                "--showvulgar",
                "yes",
                "--dpmemory",
                budget,
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]);
            assert_eq!(
                output,
                fields[6].replace("\\n", "\n"),
                "upstream fixture line {}, model {}, budget {}",
                index + 1,
                fields[0],
                budget
            );
        }
    }
    fs::remove_dir_all(directory).unwrap();
}

#[cfg(unix)]
#[test]
fn hardlinked_outputs_cannot_overwrite_inputs_or_other_reports() {
    let directory = test_directory("hardlink-protection");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    let output = directory.join("out.tsv");
    let evidence = directory.join("out.gff");
    fs::write(&query, ">q\nACGT\n").unwrap();
    fs::write(&target, ">t\nACGT\n").unwrap();
    fs::hard_link(&query, &output).unwrap();
    let error = run_failure(&[
        "--result-tsv",
        output.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(error.contains("would overwrite an input"));
    assert_eq!(fs::read_to_string(&query).unwrap(), ">q\nACGT\n");
    fs::remove_file(&output).unwrap();
    fs::write(&output, "preserve me").unwrap();
    fs::hard_link(&output, &evidence).unwrap();
    let error = run_failure(&[
        "--result-tsv",
        output.to_str().unwrap(),
        "--evidence-gff3",
        evidence.to_str().unwrap(),
        query.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(error.contains("more than one report"));
    assert_eq!(fs::read_to_string(&output).unwrap(), "preserve me");
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn batch_gff_ids_have_unambiguous_component_boundaries() {
    let directory = test_directory("gff-id-boundaries");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    let tasks = directory.join("tasks.tsv");
    let result = directory.join("out.tsv");
    let evidence = directory.join("out.gff");
    fs::write(&query, ">q\nACGT\n").unwrap();
    fs::write(&target, ">c\nACGT\n>b.c\nACGT\n").unwrap();
    fs::write(&tasks, format!("task_id\tmodel\tquery_fasta\tquery_id\ttarget_fasta\ttarget_id\na.b\tungapped\t{}\tq\t{}\tc\na\tungapped\t{}\tq\t{}\tb.c\n", query.display(), target.display(), query.display(), target.display())).unwrap();
    run(&[
        "--tasks",
        tasks.to_str().unwrap(),
        "--result-tsv",
        result.to_str().unwrap(),
        "--evidence-gff3",
        evidence.to_str().unwrap(),
        "--exhaustive",
        "yes",
        "--subopt",
        "no",
        "--revcomp",
        "no",
        "--score",
        "1",
    ]);
    let document = fs::read_to_string(evidence).unwrap();
    let mut ids = std::collections::HashSet::new();
    let mut parents = Vec::new();
    for line in document.lines().filter(|line| !line.starts_with('#')) {
        for attribute in line.split('\t').nth(8).unwrap().split(';') {
            if let Some(id) = attribute.strip_prefix("ID=") {
                assert!(ids.insert(id));
            }
            if let Some(parent) = attribute.strip_prefix("Parent=") {
                parents.push(parent);
            }
        }
    }
    assert_eq!(ids.len(), 4);
    assert!(parents.iter().all(|parent| ids.contains(parent)));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn pretty_reverse_coordinates_preserve_aligned_sequence() {
    let directory = test_directory("pretty-coordinates");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    fs::write(&query, ">q\nAAAGTC\n").unwrap();
    fs::write(&target, ">t\nGACT\n").unwrap();
    for (forward, range, sugar) in [("yes", "6 -> 2", "6 2"), ("no", "0 -> 4", "0 4")] {
        let output = run(&[
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--score",
            "10",
            "--showalignment",
            "yes",
            "--showsugar",
            "yes",
            "--forwardcoordinates",
            forward,
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert!(output.contains(&format!("Query range: {range}")));
        assert!(output.contains(&format!("sugar: q {sugar} - t 0 4 + 20")));
        assert!(output.contains("GACT\n||||\nGACT"));
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn ambiguous_codon_preserves_upstream_protein_hit() {
    let directory = test_directory("ambiguous-codon");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    fs::write(&query, ">q\nA\n").unwrap();
    fs::write(&target, ">t\nGCN\n").unwrap();
    for memory in ["32", "0"] {
        let output = run(&[
            "--model",
            "protein2dna",
            "--querytype",
            "protein",
            "--targettype",
            "dna",
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--revcomp",
            "no",
            "--score",
            "1",
            "--dpmemory",
            memory,
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ]);
        assert_eq!(output, "vulgar: q 0 1 . t 0 3 + 4 M 1 3\n");
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn all_iupac_codons_match_original_translation_table() {
    // Captured from upstream Translate_create(FALSE), standard genetic code.
    let alphabet = b"-ACGTRYSWKMBDHVNXU";
    let expected = include_str!("fixtures/iupac_translation.expected")
        .trim()
        .as_bytes();
    let mut index = 0;
    for &a in alphabet {
        for &b in alphabet {
            for &c in alphabet {
                let codon = [a, b, c];
                assert_eq!(
                    exonerate_core::translate_dna(&codon, 0),
                    [expected[index]],
                    "codon {:?}",
                    codon
                );
                assert_eq!(
                    exonerate_core::translate_dna(&codon.map(|base| base.to_ascii_lowercase()), 0),
                    [expected[index]]
                );
                index += 1;
            }
        }
    }
    assert_eq!(index, expected.len());
}

#[test]
fn est_reverse_splice_direction_matches_upstream_with_both_memory_backends() {
    let query = fixture("est-reverse-query.fa");
    let target = fixture("est-reverse-target.fa");
    for budget in ["0", "64"] {
        let common = [
            "--model",
            "est2genome",
            "--querytype",
            "dna",
            "--targettype",
            "dna",
            "--exhaustive",
            "yes",
            "--subopt",
            "no",
            "--score",
            "0",
            "--revcomp",
            "no",
            "--dpmemory",
            budget,
            "--showalignment",
            "no",
            query.to_str().unwrap(),
            target.to_str().unwrap(),
        ];
        let mut reports = common.to_vec();
        reports.extend([
            "--showsugar",
            "yes",
            "--showcigar",
            "yes",
            "--showvulgar",
            "yes",
        ]);
        assert_eq!(run(&reports), include_str!("fixtures/est-reverse.expected"));
        let mut raw = common.to_vec();
        raw.extend([
            "--showsugar",
            "no",
            "--showcigar",
            "no",
            "--showvulgar",
            "no",
            "--ryo",
            "{%Pn|%Ps|%Pqa|%Pta\\n}",
        ]);
        assert_eq!(run(&raw), include_str!("fixtures/est-reverse-ryo.expected"));
    }
}

#[test]
fn seeded_spliced_search_does_not_invent_hits_without_hsps() {
    let directory = test_directory("seedless-spliced");
    let query = directory.join("query.fa");
    let target = directory.join("target.fa");
    let cases = [
        ("est2genome", "dna", "CAGCCCTTCAACAGGGAACTCCAAACGAAGCGGGCTTGAGCCAACCTATGCATTACAAACACCTCTCGGTTACATATCCTTAAGATGGTTGCCGTATCCT".to_owned(),
         "AAGCCCTTAAACAGGGCACTCCAACCGAAGCGAGCTTGAGACAACCTAAGCATTACCAACACCTATCGGTTAAATATCCTAAAGATGGATGCCGTAACCT".to_owned()),
        ("protein2genome", "protein", "I".repeat(40), "GTG".repeat(40)),
    ];
    for (model, querytype, q, t) in cases {
        fs::write(&query, format!(">q\n{q}\n")).unwrap();
        fs::write(&target, format!(">t\n{t}\n")).unwrap();
        for memory in ["0", "64"] {
            // Original 2.4.0 emits no report for either input, even though
            // exhaustive alignment scores exceed the reporting threshold.
            assert_eq!(
                run(&[
                    "--model",
                    model,
                    "--querytype",
                    querytype,
                    "--targettype",
                    "dna",
                    "--exhaustive",
                    "no",
                    "--subopt",
                    "no",
                    "--revcomp",
                    "no",
                    "--score",
                    "100",
                    "--dpmemory",
                    memory,
                    "--showalignment",
                    "no",
                    "--showsugar",
                    "yes",
                    "--showvulgar",
                    "yes",
                    query.to_str().unwrap(),
                    target.to_str().unwrap(),
                ]),
                ""
            );
        }
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn seeded_splice_phases_frameshifts_and_raw_steps_match_upstream() {
    let directory = test_directory("seeded-extension-oracles");
    let query = directory.join("q.fa");
    let target = directory.join("t.fa");
    for (index, line) in include_str!("fixtures/seeded_extension.tsv")
        .lines()
        .enumerate()
    {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 5);
        fs::write(&query, format!(">q\n{}\n", fields[1])).unwrap();
        fs::write(&target, format!(">t\n{}\n", fields[2])).unwrap();
        for budget in ["0", "64"] {
            let output = run(&[
                "--model",
                fields[0],
                "--querytype",
                if fields[0] == "est2genome" {
                    "dna"
                } else {
                    "protein"
                },
                "--targettype",
                "dna",
                "--exhaustive",
                "no",
                "--subopt",
                "no",
                "--revcomp",
                "yes",
                "--score",
                "1",
                "--maxintron",
                fields[3],
                "--showalignment",
                "no",
                "--showsugar",
                "yes",
                "--showcigar",
                "yes",
                "--showvulgar",
                "yes",
                "--ryo",
                "{%Pn|%Ps|%Pqa|%Pta\\n}",
                "--dpmemory",
                budget,
                query.to_str().unwrap(),
                target.to_str().unwrap(),
            ]);
            assert_eq!(
                output,
                fields[4].replace("\\n", "\n"),
                "oracle {index}, budget {budget}"
            );
        }
    }
    fs::remove_dir_all(directory).unwrap();
}
