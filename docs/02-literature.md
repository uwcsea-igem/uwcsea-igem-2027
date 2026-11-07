# Literature verification

Reviewed 2026-11-07 using primary literature. ✅ means the stated, limited claim is supported by a publication, not that our team reproduced it. ⚠️ means unresolved or overextended. ❌ rejects the stated inference, not the entire project.

| Claim | Source | Status |
| --- | --- | --- |
| gp17 residues 371–553 are an established, freely replaceable recognition-only module | [Garcia-Doval & van Raaij, 2012, Structure of the receptor-binding carboxy-terminal domain, abstract and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386108/); [earlier crystallization paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3274394/) | ⚠️ The interval is a real crystallized fragment, including structural elements as well as the distal tip. ❌ Its existence does not demonstrate a universal replaceable socket. |
| A T7/T3 tail-fiber exchange has been published and can serve as a literature reference | [Ando et al., 2015, Swapping tail fibers section, Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC4785837/); [independent BASEL study, 2021](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001424) | ✅ Published engineered exchanges and a later related chimera exist. ⚠️ These establish precedent, not local reproducibility or equivalence to replacing only the proposed fragment. |
| ΔlamB / ΔfhuA / ΔbtuB alone provide a justified native T7/T3 receptor panel | [González-García et al., 2015, abstract](https://pubmed.ncbi.nlm.nih.gov/25697363/); [Etz et al., 2001, abstract](https://pmc.ncbi.nlm.nih.gov/articles/PMC95534/); [FhuA study, 2023](https://pubmed.ncbi.nlm.nih.gov/37310294/) | ❌ As a sufficient native-receptor panel, unsupported: T7 recognizes rough LPS; LamB, FhuA and BtuB are associated with other phage receptors. ⚠️ Protein-targeting intent must be confirmed with Kyungmin. The exact three-gene list comes from the user's question, not a verified strain order in the chat. |
| BindCraft binder success rates estimate success after grafting onto gp17 | [Pacesa et al., 2025, abstract and Conclusions](https://www.nature.com/articles/s41586-025-09429-6); [related T3/T7 protein-receptor study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12804162/) | ❌ Binder success is not fiber-assembly success. BindCraft reports 10–100% across tested targets, averaging 46.3%; no transferable gp17 grafting rate was established in the sources reviewed. Related display work is not direct evidence for arbitrary BindCraft grafts. ⚠️ Assembly success rate remains unknown. |

## Interpretation for the decision gate

The published exchange precedents justify retrospective evaluation, not a statement that a universal socket works. The November 4 proposal starts with a known swap, while the November 5 Level 1 is natural-phage mapping; those are different activities. The project must choose which definition it means. The audit does not identify a fatal flaw in natural-phage characterization.

Keep an explicit distinction between binding, oligomer formation, incorporation into a fiber and functional infection. A plausible model cannot establish all four. No comprehensive systematic review or independent experimental reproduction has been completed here.

## Questions for Kyungmin — draft, not sent

1. Does the knockout panel test native LPS recognition or future recognition of named protein receptors? Please share the exact strain IDs and source list.
2. Does “T3 swap” mean the published full/C-terminal exchange or a new boundary? Which primary figure supports the intended scope?
3. What observation will count as direct evidence for assembly, separately from binding and infection?
4. Which claims, if any, have already been checked against original figures and supplementary data?

## Independent strain-scoring follow-up

The [next-method feasibility review](09-next-method-feasibility.md) adds a bounded primary-source comparison of GenoPHI, Escherichia per-phage classifiers, PhageHostLearn, PhageMind, CHERRY and iPHoP, with explicit receptor, endpoint and training-data limits. GenoPHI is a conditional lead; no second method has produced independent target scores. The original four claim checks above remain in force.
