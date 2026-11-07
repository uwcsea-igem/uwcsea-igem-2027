"""Structure reference validation. Retrieval and biological suitability are separate."""
import re

def pdb_reference(code: str) -> dict:
    if not re.fullmatch(r"[1-9][A-Za-z0-9]{3}", code):
        raise ValueError("Expected a legacy four-character PDB accession")
    return {"pdb_id": code.upper(), "retrieved": False, "suitability": "unreviewed"}
