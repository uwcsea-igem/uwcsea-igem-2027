"""Track direct assembly evidence without inventing a structural score.

Binding, high monomer confidence and infectivity are not direct assembly assays.
A structural prediction backend and calibrated thresholds remain unimplemented.
"""
from dataclasses import dataclass
from typing import Literal

@dataclass(frozen=True)
class AssemblyEvidence:
    assembled: bool | None = None
    source: str = ""
    assay: str = ""

    def __post_init__(self):
        if self.assembled is not None and type(self.assembled) is not bool:
            raise ValueError("assembled must be bool or None")
        if self.assembled is not None and (not self.source.strip() or not self.assay.strip()):
            raise ValueError("Direct evidence requires source and assay")

    @property
    def status(self) -> Literal["supported", "unsupported", "unknown"]:
        if self.assembled is None:
            return "unknown"
        return "supported" if self.assembled else "unsupported"

def summarize(evidence: list[AssemblyEvidence]) -> dict:
    counts = {s: sum(e.status == s for e in evidence)
              for s in ("supported", "unsupported", "unknown")}
    return {**counts, "total": len(evidence),
            "direct_evidence_coverage": ((counts["supported"] + counts["unsupported"])
                                         / len(evidence) if evidence else None)}
