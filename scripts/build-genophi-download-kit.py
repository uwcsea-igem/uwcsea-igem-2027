"""Build a public-data-only download kit for a machine with Figshare access."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ("prepare-genophi-genomes.py", "map-genophi-training-sources.py", "audit-genophi-sources.py")


def build(metadata: Path, output: Path) -> str:
    spec = importlib.util.spec_from_file_location("prepare", ROOT / "scripts/prepare-genophi-genomes.py")
    prepare = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare)
    inventory, provenance = prepare.load_inventory(metadata)
    if output.suffix != ".zip" or output.exists():
        raise ValueError("Use a new .zip output path")
    with tempfile.TemporaryDirectory() as temporary:
        stage = Path(temporary)
        for name in SCRIPTS:
            target = stage / "scripts" / name
            target.parent.mkdir(exist_ok=True)
            target.write_bytes((ROOT / "scripts" / name).read_bytes())
        (stage / "LICENSE").write_bytes((ROOT / "LICENSE").read_bytes())
        (stage / "metadata").mkdir()
        for name in ("figshare-v1.json", "coli-tree.json"):
            (stage / "metadata" / name).write_bytes((metadata / name).read_bytes())
        (stage / "README.md").write_text("""# TipCraft public-genome download kit

This kit downloads already published genomes for a retrospective input audit.
It contains no private chat, new recognition sequences, model training or orders.

Requirements: Python 3.11+; curl on PATH; internet access; about 6 GB free disk
for the 2.07 GB dataset and transfer/archive working space. No GPU is needed.
The source inventory is 403 host assemblies and 97 phage FASTA files.

Unzip into an empty folder, open a terminal there, then run:

```sh
python3 verify-kit.py
python3 scripts/prepare-genophi-genomes.py --metadata-cache metadata --genomes genomes --output genome-integrity.json --workers 4
```

On Windows use `python` or `py -3` if that is your Python 3.11+ command.
The downloader verifies each file's deposited checksum and DNA FASTA content
before publishing it. Rerun the same command to resume verified downloads.
Exit 0 and `complete: true` mean all 500 source files passed the input audit.
Exit 1 means some files need attention; inspect `error_files` and per-file errors.
Repeated source headers are listed for later normalization; all contig sequences
and source bytes are retained without merging or redesign.

After completion, copy the `genomes` folder and `genome-integrity.json` back.
The project expects genome contents under `data/raw/genophi-genomes/`. Existing
verified files can be kept; merge without overwriting them and run the project's
downloader with `--offline` to verify the combined cache. Keep the returned
report for provenance. Neither hosts' names nor successful downloads establish
experimental-stock identity or model performance.

Sources: Figshare 10.6084/m9.figshare.25941691.v1 (CC BY 4.0) and the pinned
mdmparis/coli_phage_interactions_2023 repository. Source data retain their terms;
the original TipCraft helper scripts use the included MIT license.
""")
        (stage / "verify-kit.py").write_text("""import hashlib, json
from pathlib import Path
root = Path(__file__).resolve().parent
manifest = json.loads((root / 'kit-manifest.json').read_text())
for name, expected in manifest['file_sha256'].items():
    path = root / name
    if not path.resolve().is_relative_to(root):
        raise ValueError('Invalid kit path')
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('Kit file differs: ' + name)
print('Kit files verified; expected public genome files:', manifest['expected_genome_files'])
""")
        manifest = {"schema_version": 1, "expected_genome_files": len(inventory),
                    "expected_genome_bytes": sum(x["expected_bytes"] for x in inventory),
                    "source_provenance": provenance,
                    "file_sha256": {str(p.relative_to(stage)): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(stage.rglob("*")) if p.is_file()}}
        (stage / "kit-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    # Fixed ZIP timestamps make the kit reproducible on each build.
                    info = zipfile.ZipInfo(str(path.relative_to(stage)), date_time=(2026, 11, 7, 0, 0, 0))
                    info.compress_type = zipfile.ZIP_DEFLATED
                    info.external_attr = 0o100644 << 16
                    archive.writestr(info, path.read_bytes())
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(".zip.sha256").write_text(f"{digest}  {output.name}\n")
    return digest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(build(args.metadata_cache, args.output))
