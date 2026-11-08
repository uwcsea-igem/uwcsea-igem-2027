# V2 projection stops on missing copied-protein candidates

Updated 2026-11-08. **The v2 feature batch stopped at its seventh task because two identical copied proteins had no external-query search candidates.** Six tasks passed, one failed and eleven were never started. No biological model was fitted or scored, and the six passing tasks do not form an admissible partial accuracy evaluation.

The [failed host1-phage0 host report](../data/processed/secondary-features-v2-host1-phage0-host.json) preserves eight observed zero command exits and the failed mapping control. Exact executed commands and normalized inputs match the pre-execution plan; every listed artifact hash was verified in the [read-only audit](../data/processed/secondary-feature-v2-mapping-failure-audit.json). The batch controller exits with code 2 and its ledger records `HOLD_feature_execution_failure`. This is a technical failure before prediction, not measured host-prediction failure.

## Observed discrepancy

The copied ECOR12 proteome contains 4,962 proteins. Two fail family equality; both source/copy amino acid strings were read from the actual normalized FASTA files and verified identical.

| Training protein | Length | Training self-hit E-value | Copied-query hit | Retained family feature |
| --- | ---: | ---: | --- | --- |
| htg0000p0004747 | 46 aa | 1.572E+03 | None | None |
| htg0000p0004956 | 81 aa | 2.125E-44 | None | None |

Neither family survives the fixed requirement for occurrence in more than one training genome. The source and copied-control feature vectors therefore agree. The stronger predeclared family-and-vector control still fails; vector equality is not substituted retrospectively for that control. The failure and all six accepted receipts remain unchanged.

## Saved intermediate evidence

The [prefilter audit](../data/processed/secondary-feature-v2-prefilter-failure-audit.json) binds the saved numeric-key lookups, prefilter index files and all four data chunks for both search contexts. It records the index offsets, chunk offsets and exact entry bytes for each affected query. Both external-query entries consist of a single NUL byte and contain no candidates. Each corresponding training-self entry contains only its own numeric target key, prefilter score 65535 and diagonal 0. No new search was run to obtain this evidence.

This localizes both missing copied hits to the candidate stage. It does not establish the precise seed/masking reason for either absence. The pinned [MMseqs Alignment.cpp](https://github.com/soedinglab/MMseqs2/blob/8cc5ce367b5638c4306c2d7cfc652dd099a4643f/src/alignment/Alignment.cpp#L548) also explicitly admits identity hits outside ordinary hit criteria. The observed training self hit with E-value 1572 exceeds the logged 0.001 threshold; such a self hit is not evidence that an ordinary external query would pass. Alignment mode 3 corrected the prior identity-estimation mismatch but does not guarantee symmetric candidate availability.

## Consequence

Source/native fitting remains on technical HOLD. A further repair would require a separately recorded, label-free diagnostic and a versioned projection procedure with symmetric treatment of training and external queries, followed by full copied-proteome validation. No control is waived, failed task is retried, or phenotype is removed here. The missing candidate cause must not be inferred from the absence of final alignments alone.

Original Ando new-method coverage remains 0/10. WIsH's earlier measured negative result remains in the decision evidence, and the four exact exchange-construct cases continue to abstain. This work establishes neither binding nor tail-fiber assembly compatibility.

## Focal external-projection diagnostic

A [separately frozen six-command protocol](../software/tipcraft/benchmarks/receptor-resource/external-projection-probe-policy.json), SHA256 `d40a7fc3b2d9138b721262e603334cf7ed01c6f3303f9527fbe44ed6f7994884`, tests exactly the two discrepant amino acid sequences using independent query DB prefixes for both training-named and copied proteins. Actual source/copy sequences agree. The exact existing target database and v2 search-option vector are retained; no threshold, phenotype, vocabulary or model setting is changed. The [diagnostic helper](../scripts/probe-mmseqs-external-projection.py) verifies all bound inputs before and after the run and preserves observed commands and logs.

The [completed receipt](../data/processed/mmseqs-external-projection-probe.json) records six observed zero command exits. Both external contexts return zero hit rows for both proteins. This reproduces equal absence when neither side uses the training-self search path. It supports using symmetric projection as a candidate technical repair, but establishes neither recovered family assignments nor whole-proteome consistency.

The v2 family-and-vector control remains failed. A prospective full-proteome diagnostic must separately count equal assigned families, both-unassigned proteins and assignment discrepancies; absence cannot be renamed a successful family identification. Before production admission, any revised control contract and training projection must be versioned and validated on the complete copied proteome. No training or native scoring follows this two-protein probe.
