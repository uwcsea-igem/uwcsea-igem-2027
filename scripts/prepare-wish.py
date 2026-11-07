#!/usr/bin/env python3
"""Fetch pinned public references and compile WIsH locally; no phenotype labels."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument("--manifest", type=Path, required=True)
parser.add_argument("--workdir", type=Path, required=True)
args = parser.parse_args()
x = json.loads(args.manifest.read_text())
work = args.workdir.resolve()
work.mkdir(parents=True, exist_ok=True)

def download(url, destination):
    subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location", "--retry", "2",
                    "--max-time", "30", url, "--output", str(destination)], check=True)

for group in ("hosts", "phages"):
    (work / group).mkdir(exist_ok=True)
    for name in x[group]:
        record = x["genome_sources"][name]
        dest = work / group / (name + ".fna")
        if not dest.exists():
            download(record["url"], dest)
        lines = dest.read_text().splitlines()
        if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
            raise ValueError(f"Invalid reference FASTA for {name}")
        seq = "".join(lines[1:]).upper()
        if len(seq) != record["length"] or hashlib.sha256(seq.encode()).hexdigest() != record["sequence_sha256"]:
            raise ValueError(f"Reference hash mismatch: {name}")
        print(f"Verified {name}: {record['accession']} ({len(seq)} bp)", flush=True)
archive = work / "pinned-source.tar.gz"
if not archive.exists():
    download(x["upstream_archive_url"], archive)
with tarfile.open(archive) as tar:
    if tar.pax_headers.get("comment") != x["upstream_commit"]:
        raise ValueError("Upstream archive commit mismatch")
    source = work / tar.getmembers()[0].name.split("/")[0]
    if not source.exists():
        tar.extractall(work, filter="data")
build = work / "pinned-build"
commands = [["cmake", "-S", str(source), "-B", str(build),
             "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", "-DCMAKE_BUILD_TYPE=Release"],
            ["cmake", "--build", str(build), "-j", "2"]]
for cmd in commands:
    subprocess.run(cmd, check=True)
provenance = {"upstream_commit": x["upstream_commit"],
              "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
              "binary_sha256": hashlib.sha256((build / "WIsH").read_bytes()).hexdigest(),
              "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0],
              "build_type": "Release", "openmp_requested": "auto-detect",
              "cmake_commands": [[part.replace(str(work), "<workdir>") for part in cmd] for cmd in commands]}
(work / "build-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(f"Executable: {build / 'WIsH'}")
