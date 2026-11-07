# Isolated fitting runtime

This environment pins Python 3.12.5, CatBoost 1.2.7, scikit-learn 1.5.2, NumPy 1.26.4, pandas 2.2.2 and SciPy 1.14.1. The lockfile records the complete resolved environment and distribution hashes. The lightweight TipCraft package keeps its existing dependency-free installation.

From the repository root:

```bash
uv sync --project software/tipcraft/experiments/learning --locked --python 3.12.5
uv run --project software/tipcraft/experiments/learning --no-sync python -m unittest discover -s software/tipcraft/experiments/learning/tests -v
uv run --project software/tipcraft/experiments/learning --no-sync python scripts/check-learning-fit.py --workdir data/raw/learning-engine-check/reproduction-1 --output /tmp/learning-engine-check.json
```

Use new output paths. The runner generates its own fixed 32-row synthetic binary fixture and accepts no biological input. It actually fits RFE/CatBoost, checks frozen-bundle loading, metadata-renaming and host-only controls, and rejects changed artifacts. Model files remain in ignored raw storage. No accuracy or biological generalization estimate is calculated.

The fitting component uses fixed iterations without an evaluation set or hyperparameter search. All host columns must remain consistent per host ID, and all phage columns per phage ID. Scoring requires the externally frozen manifest hash and accepts no label vector or sample metadata. Scores are uncalibrated raw margins.

A biological runner still requires authenticated input membership, relatedness exclusions, fold-specific feature vocabularies, explicit splits and a complete frozen training specification. These are not supplied or certified by this environment or by passing its synthetic checks. See the [training-flow review](../../../../docs/18-training-protocol-and-component.md).
