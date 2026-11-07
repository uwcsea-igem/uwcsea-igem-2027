# A host encoder that cannot distinguish the target references

Updated 2026-11-07. **NO-GO for the two-reference, pattern-merged global 5-mer host encoder with union projection as the sole host input to a paired ranking model.** The three original target-host proteomes project to exactly the same nonzero vector. This is an executed representation-capacity result, not a host-infection accuracy estimate or a rejection of all GenoPHI methods.

## Frozen label-free experiment

The [policy](../software/tipcraft/benchmarks/receptor-resource/host-feature-distinguishability-policy.json) was written before fitting the first host vocabulary or projecting these target proteomes. It binds the feature kernel, source inventories and existing target annotations. Only two author-browser host proteomes fit the vocabulary; the three original reference proteomes are projection-only. No infection, receptor or donor labels enter the script. Method developers had previously seen target outcomes, so this is not a blinded study.

The profile uses within-protein 5-mers, a minimum of two distinct training-protein occurrences, training-only occurrence-pattern merging and the declared union projection. Deposited U residues are retained and X-containing windows skipped. At most one terminal stop-marker asterisk is removed as an explicit FASTA-format normalization; internal asterisks are rejected. Source sequences and annotation provenance remain unchanged.

The [executed report](../data/processed/receptor-host-feature-distinguishability.json) records **772,195 retained k-mers merged into three features**, with 768,395, 2,194 and 1,606 member k-mers, respectively:

| Declared annotated reference | Role | F1 | F2 | F3 |
| --- | --- | ---: | ---: | ---: |
| Author BL21 | Vocabulary fitting only | 1 | 0 | 1 |
| Author Keio | Vocabulary fitting only | 1 | 1 | 0 |
| Original BL21 | Projection only | 1 | 1 | 1 |
| Original BW25113 | Projection only | 1 | 1 | 1 |
| Original MG1655 | Projection only | 1 | 1 | 1 |

No target is a zero vector, and no target windows were skipped for X. The complete training vocabulary remains in ignored raw storage and is hash-bound by the versioned report. The run took about 5.6 seconds with roughly 403 MiB peak process RSS on this machine; no interaction model was fitted.

## Why a stronger classifier cannot recover the missing distinction

With two fitting samples, a present k-mer has one of three nonzero binary occurrence patterns: `10`, `01` or `11`. Merging every identical pattern therefore yields at most `2² − 1 = 3` features. This bound alone does not prove that three new hosts must collide; the actual projection above establishes that they do in this run.

Union projection activates a merged feature when any member k-mer occurs in the new proteome. Each of the three original targets activates all three groups. For a fixed phage representation `p`, a paired predictor receiving only `(p, host_features)` receives identical inputs for all three hosts. It therefore cannot rank them differently, regardless of downstream classifier complexity. Metadata or separately selected host models would be additional inputs outside this declared sequence-feature configuration.

The initial source and target annotations have different provenance: author-browser deposited proteins versus the existing Prodigal 2.6.3 calls on the original public references. The fixed uniform-caller control below tests whether that difference explains this collision. The result remains conditional on these reference inputs and this encoder; it does not establish biological or stock equivalence, universal encoder failure or equivalence to the authors' novel-genome assignment workflow.

## Uniform-caller control

After observing the initial collision, we froze a [single control policy](../software/tipcraft/benchmarks/receptor-resource/host-feature-uniform-control-policy.json) before new calls or feature projection. We applied the **same Prodigal 2.6.3 binary and single-mode/table-11 parameters** as the original three target annotations to the two author reference genomes. The three original target calls were reused unchanged. The 5-mer, rare-filter, merging and union-projection rules were unchanged; no parameter search, residue substitution or interaction model was introduced.

The [executed control](../data/processed/receptor-host-feature-uniform-control.json) gives:

| Host-reference input representation | Two fitting proteomes: protein count | Retained k-mers | Merged features | Each original target vector |
| --- | ---: | ---: | ---: | --- |
| Author-deposited annotations | 8,341 | 772,195 | 3 | `[1, 1, 1]` |
| Uniform Prodigal 2.6.3 annotations | 8,587 | 787,904 | 3 | `[1, 1, 1]` |

Caller uniformity **did not remove the observed collision**. The initial and control reports, deposited annotations and reference genomes remain separate and unchanged. This control strengthens the NO-GO for this encoder on the current panel; it does not show that all possible annotations or encoders must collide. Its timing is explicitly post hoc relative to the initial result, with its own parameters fixed before control execution.

## Consequence for the December gate

Do not spend a biological model-fit run on this host configuration and then treat ties or host priors as strain-specific evidence. Preserve this negative result. A future candidate needs a separately specified representation that retains a distinction among the actual target inputs, followed by source admission, grouped validation and matched-subset baselines. Changing an encoder to preserve input information is not evidence that its infection predictions will improve.

The receptor-only annotation backend remains a different endpoint. The original WIsH failure, six native observations, four construct-source abstentions and ten-case denominator remain unchanged. New independent-method coverage remains **0/10**; there is no new infection, binding or assembly score.

## Reproduction

```bash
uv run --project software/tipcraft python scripts/check-receptor-host-features.py --workdir data/raw/receptor-resource/host-feature-capacity-reproduction-1 --output /tmp/receptor-host-feature-distinguishability.json
```

The target calls can be restored from the public annotation receipt and the three original reference genomes alone:

```bash
uv run --project software/tipcraft python scripts/restore-reference-host-proteomes.py
```

The [fresh restoration check](../data/processed/reference-host-proteome-restoration-check.json) actually reconstructed all three target protein FASTAs and reproduced their recorded hashes exactly. No Picard host collection or old phage corpus was needed. The copied historical annotation receipt supplies provenance; it is not a claim that its full collection was rerun. Verified protein caches can be used without invoking the caller; fresh calls require the pinned executable.

Reproduce the uniform-caller control separately:

```bash
uv run --project software/tipcraft python scripts/check-receptor-host-uniform-control.py --workdir data/raw/receptor-resource/host-feature-uniform-reproduction-1 --output /tmp/receptor-host-feature-uniform-control.json
```

Use fresh output paths and the authenticated caches. Tests reproduce the collision with distinct synthetic proteomes without phenotypes and separately verify that a zero vector is not mistaken for the observed nonzero collision.
