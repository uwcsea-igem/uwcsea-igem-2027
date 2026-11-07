"""Run the external WIsH scorer without loading phenotype labels.

WIsH (Galiez et al., 2017) is GPL-3.0; no upstream source is copied here.
This wrapper only launches its executable and parses numeric output.
"""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Expected FASTA")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("This benchmark requires one complete sequence per FASTA")
    bases = "".join(lines[1:]).upper()
    if not bases or set(bases) - set("ACGT"):
        raise ValueError("This run requires unambiguous ACGT reference sequences")
    return bases


def parse_matrix(path: Path) -> dict[str, dict[str, float]]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("Empty WIsH matrix")
    phages = lines[0].split()
    if len(set(phages)) != len(phages):
        raise ValueError("Duplicate phage columns")
    result = {phage: {} for phage in phages}
    hosts = set()
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != len(phages) + 1 or fields[0] in hosts:
            raise ValueError("Malformed matrix or duplicate host")
        hosts.add(fields[0])
        for phage, value in zip(phages, fields[1:]):
            score = float(value)
            if not math.isfinite(score):
                raise ValueError("Nonfinite WIsH likelihood")
            result[phage][fields[0]] = score
    return result


def run_scores(manifest_path: Path, genome_dir: Path, executable: Path,
               output_dir: Path, anonymize: bool = False) -> dict:
    # Explicit inputs only: there is no cases/labels argument or phenotype database.
    manifest = json.loads(manifest_path.read_text())
    protected = {manifest_path.resolve(), executable.resolve()}
    protected.update((genome_dir / group / (name + ".fna")).resolve()
                     for group in ("hosts", "phages") for name in manifest[group])
    outputs = {(output_dir / name).resolve() for name in
               ("scores.json", "llikelihood.matrix", "prediction.list", "run.log")}
    if protected & outputs:
        raise ValueError("Output cannot overwrite scoring inputs")
    for group in ("hosts", "phages"):
        for name in manifest[group]:
            bases = sequence(genome_dir / group / (name + ".fna"))
            expected = manifest["genome_sources"][name]
            if (len(bases) != expected["length"] or
                hashlib.sha256(bases.encode()).hexdigest() != expected["sequence_sha256"]):
                raise ValueError(f"Reference sequence mismatch: {name}")
    executable = executable.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tipcraft-wish-") as temp:
        work = Path(temp)
        inverse = {}
        for group in ("hosts", "phages"):
            (work / group).mkdir()
            inverse[group] = {}
            for i, name in enumerate(manifest[group]):
                alias = f"{group[0]}{i}" if anonymize else name
                inverse[group][alias] = name
                # Replace header even in the named run: metadata never reaches scoring.
                bases = sequence(genome_dir / group / (name + ".fna"))
                (work / group / (alias + ".fna")).write_text(f">sequence\n{bases}\n")
        for directory in ("models", "results"):
            (work / directory).mkdir()
        commands = [
            [str(executable), "-c", "build", "-g", str(work / "hosts"),
             "-m", str(work / "models"), "-k", str(manifest["order"]),
             "-a", str(manifest["alpha"]), "-t", "1"],
            [str(executable), "-c", "predict", "-g", str(work / "phages"),
             "-m", str(work / "models"), "-r", str(work / "results"), "-b", "3", "-t", "1"],
        ]
        logs = []
        for command in commands:
            done = subprocess.run(command, capture_output=True, text=True, check=True)
            logs.append(done.stdout + done.stderr)
        matrix_file = work / "results/llikelihood.matrix"
        raw = parse_matrix(matrix_file)
        scores = {inverse["phages"][p]: {inverse["hosts"][h]: v for h, v in hosts.items()}
                  for p, hosts in raw.items()}
        (output_dir / "llikelihood.matrix").write_bytes(matrix_file.read_bytes())
        (output_dir / "prediction.list").write_bytes((work / "results/prediction.list").read_bytes())
        # Replace temporary directory names to make logs inspectable and stable.
        (output_dir / "run.log").write_text("\n".join(logs).replace(str(work), "<workdir>"))
    result = {
        "method": "WIsH upstream 1.1", "order": manifest["order"], "alpha": manifest["alpha"],
        "upstream_commit": manifest["upstream_commit"],
        "executable_sha256": sha256(executable), "manifest_sha256": sha256(manifest_path),
        "score_type": "mean conditional log likelihood; higher is better; not probability",
        "anonymized": anonymize,
        "scores": {p: dict(sorted(s.items())) for p, s in sorted(scores.items())},
        "excluded_constructs": manifest["excluded_constructs"],
        "no_labels_loaded": True,
    }
    (output_dir / "scores.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
