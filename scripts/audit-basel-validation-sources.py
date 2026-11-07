"""Map published secondary-panel IDs to accessions without inspecting outcomes."""
import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

UPSTREAM_COMMIT = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"
MATRIX_PATH = "data/experimental_validation/BASEL_ECOR_interaction_matrix.csv"
MATRIX_BLOB = "dcbbca4b4e4b1dd0bc543c5e2348189813e55bb0"
MATRIX_SHA256 = "bc0426a80b392281815ab24066e97cdc445399a62d67bb4f8e56511088947ba5"
ARCHIVE_MD5 = "157de25b9c075a93b94781db04182d7b"
ARCHIVE_SHA256 = "1d98c77fc4068315246a4eb56aa185fc769c06bf70fd20e5090856b25f94a00d"
TRAIN_SHA256 = "c4795d786a1e13ac998e263dbf4aaf95fc0f94c52b842d02b6f7fc6e99ccfef2"
MEMBERS = {
    "phages": ("GenoPHI_processed_data/paper_supplementary_tables/SuppTable8_phage_validation_list.csv",
               "889a17c9188b4ef96292b48bdcebeff92a2ac39e1fc549d9c6391a607881b443"),
    "hosts": ("GenoPHI_processed_data/ECOR_validation_strain_assemblies.csv",
              "fb5ff1a75f4caaeafb0ea1a48361a68572c6249b2f391d86d9b7773d58774927"),
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def unique_mapping(rows, identifier, accession):
    result = {}
    for row in rows:
        name, record = row[identifier].strip(), row[accession].strip()
        if not name or not record or name in result:
            raise ValueError("Empty or ambiguous metadata identifier/accession")
        result[name] = record
    return result


def reconcile(matrix_rows, phage_rows, host_rows, training_rows):
    if not matrix_rows or not training_rows:
        raise ValueError("Panel and training rosters must contain records")
    pairs = [(row["strain"], row["phage"]) for row in matrix_rows]
    if len(set(pairs)) != len(pairs) or any(not a.strip() or not b.strip() for a, b in pairs):
        raise ValueError("Empty identifier or duplicate panel pair")
    hosts, phages = sorted({a for a, _ in pairs}), sorted({b for _, b in pairs})
    pm = unique_mapping(phage_rows, "Phage ID", "Genome Accession")
    hm = unique_mapping(host_rows, "Strain", "NCBI accession")
    if set(phages) - pm.keys() or set(hosts) - hm.keys():
        raise ValueError("Missing exact panel-ID accession mapping")
    flagged = {r["Phage ID"].strip() for r in phage_rows if r["Interaction matrix"].strip()}
    train_hosts = {row["strain"] for row in training_rows}
    train_phages = {row["phage"] for row in training_rows}
    accession_ids = defaultdict(list)
    for name in phages:
        accession_ids[pm[name]].append(name)
    return {
        "panel_pairs": len(pairs), "panel_hosts": len(hosts), "panel_phages": len(phages),
        "panel_roster_authority": "Every ID in the authenticated matrix, independent of outcome values and metadata flags",
        "phage_id_to_accession": {name: pm[name] for name in phages},
        "host_id_to_wgs_master_accession": {name: hm[name] for name in hosts},
        "metadata_phage_rows": len(phage_rows), "metadata_flagged_interaction_phages": len(flagged),
        "flagged_phages_absent_from_matrix": sorted(flagged - set(phages)),
        "matrix_phages_without_metadata_flag": sorted(set(phages) - flagged),
        "metadata_host_rows": len(host_rows), "metadata_hosts_absent_from_matrix": sorted(hm.keys() - set(hosts)),
        "duplicate_panel_phage_accession_groups": {a: ids for a, ids in accession_ids.items() if len(ids) > 1},
        "training_host_id_overlap_count": len(set(hosts) & train_hosts),
        "training_host_id_overlap": sorted(set(hosts) & train_hosts),
        "training_phage_id_overlap_count": len(set(phages) & train_phages),
        "training_phage_id_overlap": sorted(set(phages) & train_phages),
    }


def run(matrix, archive, training):
    mb, ab, tb = matrix.read_bytes(), archive.read_bytes(), training.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(mb)).encode() + b"\0" + mb).hexdigest()
    if blob != MATRIX_BLOB or digest(mb) != MATRIX_SHA256:
        raise ValueError("Panel matrix differs from pinned source")
    if hashlib.md5(ab).hexdigest() != ARCHIVE_MD5 or digest(ab) != ARCHIVE_SHA256:
        raise ValueError("Processed archive differs from reviewed deposit")
    if digest(tb) != TRAIN_SHA256:
        raise ValueError("Training table differs from pinned source")
    member_rows = {}
    with zipfile.ZipFile(io.BytesIO(ab)) as z:
        for key, (name, expected) in MEMBERS.items():
            data = z.read(name)
            if digest(data) != expected:
                raise ValueError("Archive metadata member differs from reviewed source")
            member_rows[key] = list(csv.DictReader(io.StringIO(data.decode())))
    result = reconcile(list(csv.DictReader(io.StringIO(mb.decode()))), member_rows["phages"],
                       member_rows["hosts"], list(csv.DictReader(io.StringIO(tb.decode()))))
    if (result["panel_pairs"], result["panel_hosts"], result["panel_phages"]) != (1240, 25, 52):
        raise ValueError("Unexpected published panel dimensions")
    return {"schema_version": 1, "status": "FUTURE_SECONDARY_PANEL_NOT_SELECTED_FOR_CURRENT_TARGET_EVALUATION",
            "purpose": __doc__, "helper_sha256": digest(Path(__file__).read_bytes()),
            "upstream_commit": UPSTREAM_COMMIT,
            "matrix_url": f"https://raw.githubusercontent.com/Noonanav/GenoPHI/{UPSTREAM_COMMIT}/{MATRIX_PATH}",
            "matrix_git_blob_sha1": blob, "matrix_sha256": digest(mb),
            "archive_record": "https://zenodo.org/records/21696332", "archive_md5": ARCHIVE_MD5,
            "archive_sha256": digest(ab),
            "archive_metadata_members": {key: {"path": name, "sha256": value} for key, (name, value) in MEMBERS.items()},
            "training_table_sha256": digest(tb), **result,
            "assay_source": "https://doi.org/10.1038/s41564-026-02482-5",
            "assay_locator": "Validation on unseen phages and host genes; Methods: Experimental validation of model performance; Figure 4 caption",
            "assay_interpretation": "Published spot-clearance/no-clearance observations. Paper reports 1300 tested pairs and60 unclear outcomes excluded, leaving1240; clearance need not establish productive infection. Do not fill absent pairs as negatives.",
            "validation_domain": "Method's already-published novel-phage/known-host validation panel; our reuse would be retrospective, not a fresh independent method-development holdout",
            "accession_versions_and_sequence_checksums_frozen": False,
            "full_sequence_identity_or_homology_audit_run": False,
            "host_stock_identity_verified": False, "phage_sequence_novelty_verified": False,
            "outcome_values_used_for_roster_selection_or_tuning": False,
            "prediction_output_archive_members_inspected": False,
            "genomes_downloaded": False, "training_or_inference_executed": False,
            "current_t7_t3_evaluation_replaced": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("matrix", "archive", "training-table", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve() in {args.matrix.resolve(), args.archive.resolve(), args.training_table.resolve()}:
        parser.error("Use a new report path separate from sources")
    report = run(args.matrix, args.archive, args.training_table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
