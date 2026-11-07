"""Compare original pinned VIRIDIC R execution with frozen adapter outputs.

No alignments are launched. The original GPL function is parsed and evaluated
from its authenticated source file, with actual upstream package APIs.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tarfile


R_DRIVER = r'''
options(digits=17, warn=1)
suppressPackageStartupMessages({
  library(magrittr); library(dplyr); library(seqinr)
  library(stringr); library(IRanges); library(data.table)
})
args <- commandArgs(trailingOnly=TRUE)
expressions <- parse(args[1])
selected <- which(vapply(expressions, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("nident_fun"))
}, logical(1)))
stopifnot(length(selected) == 1)
# Evaluate the original expression unchanged; never source the full workflow.
eval(expressions[[selected]], envir=.GlobalEnv)
saveRDS(nident_fun, paste0(args[3], "/original-function.rds"))
DF <- data.table::fread(args[2], sep="\t", header=FALSE, drop=c(3,4))
colnames(DF) <- c("qseqid","sseqid","qlen","slen","qstart","qend",
                 "sstart","send","qseq","sseq","nident","gaps")
groups <- unique(DF[, .(qseqid,sseqid,qlen,slen)])
indices <- split(seq_len(nrow(DF)), paste(DF$qseqid, DF$sseqid, sep="\t"))
counts <- data.frame(query=groups$qseqid, subject=groups$sseqid,
                     qlen=groups$qlen, slen=groups$slen,
                     nident=NA_real_, original_aligned_query_bases=NA_real_,
                     error=rep("",nrow(groups)))
for (i in seq_len(nrow(groups))) {
  ix <- indices[[paste(groups$qseqid[i], groups$sseqid[i], sep="\t")]]
  input <- tibble::as_tibble(DF[ix, .(qstart,qend,sstart,send,qseq,sseq,nident,gaps)])
  result <- tryCatch(nident_fun(input, groups$qlen[i]), error=function(e) e)
  if (inherits(result,"error")) counts$error[i] <- conditionMessage(result) else {
    counts$nident[i] <- result[1]
    counts$original_aligned_query_bases[i] <- result[2]
  }
}
write.table(counts, paste0(args[3],"/directional.tsv"), sep="\t", quote=TRUE, row.names=FALSE)
lengths <- read.delim(paste0(args[3],"/lengths.tsv"), stringsAsFactors=FALSE)
keys <- paste(counts$query, counts$subject, sep="\t")
get_count <- function(a,b) {
  i <- match(paste(a,b,sep="\t"),keys)
  if (is.na(i)) 0 else counts$nident[i]
}
pair_results <- vector("list",nrow(lengths)*(nrow(lengths)-1)/2)
k <- 0
for (i in seq_len(nrow(lengths)-1)) for (j in seq.int(i+1,nrow(lengths))) {
  a <- lengths$id[i]; b <- lengths$id[j]
  sim <- ((get_count(a,b)+get_count(b,a))*100)/(lengths$length[i]+lengths$length[j])
  k <- k+1
  pair_results[[k]] <- data.frame(query=a,subject=b,similarity_percent_unrounded=sim,
                                  similarity_percent=round(sim,digits=3))
}
write.table(do.call(rbind,pair_results), paste0(args[3],"/pairs.tsv"),
            sep="\t", quote=TRUE, row.names=FALSE)
capture.output(sessionInfo(), file=paste0(args[3],"/session-info.txt"))
packages <- c("magrittr","dplyr","tibble","seqinr","stringr","IRanges","S4Vectors","BiocGenerics","data.table")
package_info <- data.frame(package=packages,
                          version=vapply(packages,function(p) as.character(packageVersion(p)),character(1)))
write.table(package_info,paste0(args[3],"/packages.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_runtime(r_home, r_library, manifest):
    """Bind actual supplied runtime and installed packages to authenticated bytes."""
    for path, key in ((r_home/"bin/exec/R","r_executable_sha256"),
                      (r_home/"lib/libR.dylib","libR_sha256")):
        if sha(path) != manifest[key]:
            raise ValueError("Actual supplied R runtime differs from authenticated manifest: "+str(path))
    checked_files = 0
    packages = manifest["packages"]
    if not packages or len({p["name"] for p in packages}) != len(packages):
        raise ValueError("Runtime package inventory is empty or duplicated")
    for package in packages:
        archive = Path(manifest["runtime_root"])/"packages"/f"{package['name']}_{package['version']}.tgz"
        if not package.get("md5_verified") or sha(archive) != package["sha256"]:
            raise ValueError("Authenticated R package archive changed: "+package["name"])
        description = r_library/package["name"]/"DESCRIPTION"
        versions = [line.removeprefix("Version:").strip() for line in description.read_text().splitlines()
                    if line.startswith("Version:")]
        if versions != [package["version"]]:
            raise ValueError("Installed R package DESCRIPTION version changed: "+package["name"])
        with tarfile.open(archive) as bundle:
            for member in bundle.getmembers():
                relative = Path(member.name)
                if relative.is_absolute() or ".." in relative.parts or relative.parts[0] != package["name"]:
                    raise ValueError("Unexpected authenticated R package archive path")
                installed = r_library/relative
                if member.isfile():
                    with bundle.extractfile(member) as stream:
                        expected = hashlib.file_digest(stream,"sha256").hexdigest()
                    if not installed.is_file() or sha(installed) != expected:
                        raise ValueError("Installed R package file differs from authenticated archive: "+member.name)
                    checked_files += 1
                elif member.issym():
                    if not installed.is_symlink() or os.readlink(installed) != member.linkname:
                        raise ValueError("Installed R package symlink differs from authenticated archive")
                elif not member.isdir():
                    raise ValueError("Unsupported authenticated R package archive entry")
    return {"actual_R_executable_sha256":sha(r_home/"bin/exec/R"),
            "actual_libR_sha256":sha(r_home/"lib/libR.dylib"),
            "installed_packages_verified":len(packages), "installed_regular_files_verified":checked_files,
            "verification":"Actual supplied executable/libR hashes; DESCRIPTION versions; every regular installed package file and symlink checked against the authenticated binary archives"}


def compare(policy, workdir):
    lengths = json.loads(Path(policy["inputs"]["genome_lengths"]["path"]).read_text())
    expected = json.loads(Path(policy["inputs"]["adapter_output"]["path"]).read_text())["pairs"]
    original = read_tsv(workdir / "pairs.tsv")
    counts = read_tsv(workdir / "directional.tsv")
    observed = {}
    for row in counts:
        key = row["query"], row["subject"]
        if key in observed or key[0] not in lengths or key[1] not in lengths:
            raise ValueError("Duplicate or undeclared original R directional key")
        if int(row["qlen"]) != lengths[key[0]] or int(row["slen"]) != lengths[key[1]]:
            raise ValueError("Original R lengths differ from frozen genomes")
        observed[key] = row
    errors = [row for row in counts if row["error"]]
    expected_pairs = {(p["query"],p["subject"]):p for p in expected}
    original_pairs = {(p["query"],p["subject"]):p for p in original}
    ids = sorted(lengths)
    all_keys = {(a,b) for i,a in enumerate(ids) for b in ids[i+1:]}
    if (len(expected) != policy["scope"]["distinct_unordered_pairs"] or len(original) != len(expected)
            or set(expected_pairs) != all_keys or set(original_pairs) != all_keys):
        raise ValueError("Original/adapter unordered pair matrices are incomplete or duplicated")
    discrepancies = []
    max_difference = 0.0
    directions_checked = 0
    for key in sorted(all_keys):
        exp, actual = expected_pairs[key], original_pairs[key]
        for a,b,prefix in ((key[0],key[1],"query"),(key[1],key[0],"subject")):
            row = observed.get((a,b))
            directions_checked += 1
            if row is not None and row["error"]:
                discrepancies.append({"query":a,"subject":b,"type":"original_R_error","error":row["error"]})
                continue
            for output, expected_field in (("nident",prefix+"_nident"),
                                           ("original_aligned_query_bases",prefix+"_original_aligned_bases")):
                value = float(row[output]) if row is not None else 0.0
                if not math.isfinite(value) or not value.is_integer() or value != exp[expected_field]:
                    discrepancies.append({"query":a,"subject":b,"type":output,
                                          "original":value,"adapter":exp[expected_field]})
        try:
            value = float(actual["similarity_percent_unrounded"])
            rounded = float(actual["similarity_percent"])
            expected_value = float(exp["similarity_percent_unrounded"])
            expected_rounded = float(exp["similarity_percent"])
        except (TypeError, ValueError):
            discrepancies.append({"query":key[0],"subject":key[1],"type":"invalid_numeric_similarity",
                                  "original_unrounded":actual["similarity_percent_unrounded"],
                                  "original_rounded":actual["similarity_percent"]})
            continue
        if not all(math.isfinite(x) for x in (value,rounded,expected_value,expected_rounded)):
            discrepancies.append({"query":key[0],"subject":key[1],"type":"nonfinite_numeric_similarity",
                                  "original_unrounded":str(actual["similarity_percent_unrounded"]),
                                  "original_rounded":str(actual["similarity_percent"])})
            continue
        difference = abs(value-expected_value)
        max_difference = max(max_difference,difference)
        if not math.isfinite(value) or difference > policy["agreement_requirements"]["unrounded_similarity_absolute_tolerance_percentage_points"]:
            discrepancies.append({"query":key[0],"subject":key[1],"type":"unrounded_similarity",
                                  "original":value,"adapter":exp["similarity_percent_unrounded"]})
        if rounded != expected_rounded:
            discrepancies.append({"query":key[0],"subject":key[1],"type":"three_decimal_similarity",
                                  "original":float(actual["similarity_percent"]),"adapter":exp["similarity_percent"]})
    self_failures = []
    for identifier,length in sorted(lengths.items()):
        row = observed.get((identifier,identifier))
        if (row is None or row["error"] or float(row["nident"]) != length
                or float(row["original_aligned_query_bases"]) != length):
            self_failures.append(identifier)
    passed = not discrepancies and not errors and not self_failures
    return {"status":"passed" if passed else "discrepancy", "original_R_function_executed":True,
            "unordered_pairs_compared":len(original), "ordered_nonself_directions_compared":directions_checked,
            "original_nonempty_directional_groups":len(counts), "self_controls_checked":len(lengths),
            "self_control_failures":self_failures, "original_R_errors":errors,
            "maximum_unrounded_difference_percentage_points":max_difference,
            "discrepancies":discrepancies, "frozen_dataset_directional_equivalence_verified":passed,
            "universal_upstream_equivalence_verified":False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("policy","r-home","r-library","runtime-manifest","workdir","output"):
        parser.add_argument("--"+name,type=Path,required=True)
    args = parser.parse_args()
    if args.workdir.exists() or args.output.exists():
        parser.error("Use new comparison outputs; preserve existing reports")
    policy = json.loads(args.policy.read_text())
    for item in policy["inputs"].values():
        if sha(item["path"]) != item["sha256"]:
            parser.error("Comparison input differs from frozen policy: "+item["path"])
    if policy["runtime_manifest_sha256"] != sha(args.runtime_manifest) or policy["comparison_helper_sha256"] != sha(Path(__file__)):
        parser.error("Runtime/helper differs from frozen comparison policy")
    runtime = json.loads(args.runtime_manifest.read_text())
    runtime_verification = validate_runtime(args.r_home,args.r_library,runtime)
    length_bindings = policy["equivalent_genome_length_inputs"]
    if len(length_bindings) < 2:
        parser.error("Both normalized length inputs must be explicitly bound")
    for item in length_bindings:
        if sha(item["path"]) != item["sha256"] or item["sha256"] != policy["inputs"]["genome_lengths"]["sha256"]:
            parser.error("Normalized length-input bindings differ")
    args.workdir.mkdir(parents=True)
    driver = args.workdir / "run-original.R"
    driver.write_text(R_DRIVER)
    lengths = json.loads(Path(policy["inputs"]["genome_lengths"]["path"]).read_text())
    (args.workdir/"lengths.tsv").write_text("id\tlength\n"+"".join(f"{k}\t{v}\n" for k,v in sorted(lengths.items())))
    environment = os.environ.copy()
    environment.update(R_HOME=str(args.r_home.resolve()), R_LIBS_USER=str(args.r_library.resolve()),
                       DYLD_LIBRARY_PATH=str((args.r_home/"lib").resolve()))
    command = [str(args.r_home/"bin/exec/R"),"--vanilla","--slave","-f",str(driver),"--args",
               policy["inputs"]["upstream_source"]["path"],policy["inputs"]["blast_tsv"]["path"],str(args.workdir)]
    try:
        result = subprocess.run(command,env=environment,capture_output=True,text=True,timeout=600)
        (args.workdir/"execution.log").write_text(result.stdout+result.stderr)
        comparison = compare(policy,args.workdir) if result.returncode == 0 else {
            "status":"blocked", "original_R_function_executed":False,
            "reason":"R comparison process failed; inspect execution.log", "returncode":result.returncode,
            "frozen_dataset_directional_equivalence_verified":False,"universal_upstream_equivalence_verified":False}
    except subprocess.TimeoutExpired as exc:
        (args.workdir/"execution.log").write_text(str(exc))
        comparison = {"status":"blocked","reason":"R calculation exceeded the bounded 600-second runtime",
                      "frozen_dataset_directional_equivalence_verified":False,"universal_upstream_equivalence_verified":False}
    report = {"schema_version":1, "purpose":policy["purpose"], "policy_sha256":sha(args.policy),
              "runtime_manifest_sha256":sha(args.runtime_manifest), "comparison_helper_sha256":sha(Path(__file__)),
              "runtime_manifest":runtime, "runtime_content_verification":runtime_verification,
              "equivalent_genome_length_inputs":length_bindings, "command":command,
              "environment_overrides":{k:environment[k] for k in ("R_HOME","R_LIBS_USER","DYLD_LIBRARY_PATH")},
              "inputs":policy["inputs"], "comparison":comparison,
              "raw_artifacts":{p.name:sha(p) for p in sorted(args.workdir.iterdir()) if p.is_file()},
              "blast_rerun":False, "parameters_retuned":False, "phenotype_labels_loaded":False,
              "infection_or_assembly_predictions_produced":False, "interpretation":policy["interpretation"]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    return 0 if comparison["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
