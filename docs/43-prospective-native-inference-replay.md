# Prospective native inference replay

Recorded 2026-11-08, before any biological model fit or score. The [replay protocol](../software/tipcraft/benchmarks/receptor-resource/native-inference-replay-policy.json) fixes two checks on all six native pairs and all three saved model scopes: a fresh-process numeric replay in the original row order, and a separate replay using opaque row identifiers, a new input filename and a fixed row permutation. Both require exact finite raw-margin equality and unchanged authenticated model artifacts. Differences and failed process exits will be preserved. No tolerance or model setting will be selected from replay outcomes.

These checks operate at the inference boundary. They reuse the frozen sequence-derived features and saved models; they do not rerun annotation, protein searches, vocabulary construction or training. The scorer input carries no target phenotype vector. Its runtime is not a security boundary, and model provenance remains readable. Passing these checks would establish reproducible inference under the specified input transformation, not an anonymous end-to-end experiment, label inaccessibility, biological accuracy or blinding.

A complete actual-input/command manifest will precede replay execution. The replay has not run, no biological scores exist at this checkpoint, and new-method Ando coverage remains 0/10. The four exact exchange-construct cases continue to abstain. This supplementary procedure changes neither the frozen primary evaluation nor the source/native fitting parameters.

## Executable worker validation

The [inference-only worker](../scripts/replay-native-inference-worker.py) accepts an externally hashed packet with exactly six opaque rows and exactly the three fixed model scopes. It rejects additional packet fields, changed bundle manifests, wrong scope bindings, invalid binary feature matrices, and missing or nonfinite scores. It preserves row order and checks packet and manifest hashes again after inference. The caller remains responsible for authenticating the complete source/native run, supplying the actual runtime/command manifest, and restoring the separately saved identity mapping before comparing scores.

Four focused packet tests pass. The locked learning-runtime native integration test also passes with real saved CatBoost bundles trained on synthetic data: two separate worker processes load those bundles and reproduce all 18 saved raw margins exactly in original order and under the prospectively fixed permutation. These are technical execution checks, not a biological replay or a new host-prediction result. Existing feature/fitting/evaluation consumers and their frozen settings are unchanged.

Worker SHA256 at this checkpoint: `bb72d0ace6ab29cd1c897968e2553f02117eeb8ddcd1a688a81d7b36188be22d`. The actual biological replay manifest must bind this implementation or explicitly record a revised implementation before execution.
