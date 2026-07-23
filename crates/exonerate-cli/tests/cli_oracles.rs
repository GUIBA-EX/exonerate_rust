use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(args: &[&str]) -> String {
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

fn run_failure(args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_exonerate-rs"))
        .args(args)
        .output()
        .expect("run exonerate-rs");
    assert!(
        !output.status.success(),
        "exonerate-rs unexpectedly succeeded"
    );
    String::from_utf8(output.stderr).expect("CLI stderr is UTF-8")
}

#[test]
fn help_and_version_aliases_are_available_without_fasta_inputs() {
    for argument in ["-h", "--shorthelp", "--help"] {
        let output = run(&[argument]);
        assert!(output.starts_with("Usage: exonerate-rs"));
    }
    for argument in ["-v", "--version"] {
        let output = run(&[argument]);
        assert_eq!(output, "exonerate-rs 0.1.0\n");
    }
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
