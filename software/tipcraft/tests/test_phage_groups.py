import importlib.util
from decimal import Decimal
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location("phage_groups", Path(__file__).resolve().parents[3] / "scripts/group-genophi-phages.py")
GROUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GROUP)


def node(name, role="external_collection"):
    return {"id": "n_" + name, "source_id": name, "role": role, "retained_bases": 100000000}


def matrix(nodes, values=None, reverse=False):
    values = values or {}
    rows = []
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            value = values.get(frozenset((a["source_id"], b["source_id"])), 0)
            total = int(Decimal(str(value)) * Decimal(a["retained_bases"] + b["retained_bases"]) / 100)
            ai, bi = total // 2, total - total // 2
            actual = 100 * total / (a["retained_bases"] + b["retained_bases"])
            rows.append("\t".join(map(str, [a["id"], b["id"], a["source_id"], b["source_id"],
                                            actual, ai, bi, a["retained_bases"], b["retained_bases"]])))
    if reverse:
        rows.reverse()
    return "\t".join(GROUP.FIELDS) + "\n" + "\n".join(rows) + "\n"


class PhageGroupTests(unittest.TestCase):
    def test_transitive_target_and_near_gp17_exclusions_keep_nonmodel_bridges(self):
        nodes = [node("T7", "target_reference"), node("bridge"), node("modelA"),
                 node("seed"), node("modelB"), node("retained"), node("T3", "target_reference")]
        edges = {frozenset(("T7", "bridge")): 80, frozenset(("bridge", "modelA")): 80,
                 frozenset(("seed", "modelB")): 80}
        result = GROUP.group_phages(matrix(nodes, edges), nodes, ["modelA", "modelB", "retained"], ["seed"])
        self.assertEqual(result["excluded_model_table_phages"], ["modelA", "modelB"])
        self.assertEqual(result["retained_model_table_phages"], ["retained"])
        ga = result["model_phage_groups"]["modelA"]
        self.assertEqual(result["groups"][ga]["source_members"], ["T7", "bridge", "modelA"])
        gb = result["model_phage_groups"]["modelB"]
        self.assertEqual(result["groups"][gb]["near_gp17_seed_members"], ["seed"])
        self.assertEqual(result["group_count"], 4)
        self.assertIn(result["target_phage_groups"]["T3"], result["groups"])

    def test_exact_buffer_boundary_is_included_and_below_boundary_is_separate(self):
        nodes = [node("A"), node("B"), node("C")]
        values = {frozenset(("A", "B")): 69.99, frozenset(("B", "C")): 69.989999}
        result = GROUP.group_phages(matrix(nodes, values), nodes, ["A", "B", "C"], [])
        self.assertEqual(result["source_phage_groups"]["A"], result["source_phage_groups"]["B"])
        self.assertNotEqual(result["source_phage_groups"]["B"], result["source_phage_groups"]["C"])
        self.assertEqual(len(result["boundary_edges"]), 1)

    def test_node_matrix_model_and_seed_reordering_do_not_change_groups_or_ids(self):
        nodes = [node("target", "target_reference"), node("bridge"), node("A"), node("B")]
        values = {frozenset(("target", "bridge")): 75, frozenset(("bridge", "A")): 75}
        first = GROUP.group_phages(matrix(nodes, values), nodes, ["A", "B"], ["B"])
        second = GROUP.group_phages(matrix(nodes[::-1], values, reverse=True), nodes[::-1], ["B", "A"], ["B"])
        self.assertEqual(first, second)

    def test_incomplete_duplicate_unknown_nonfinite_and_inconsistent_matrix_rejected(self):
        nodes = [node("A"), node("B"), node("C")]
        valid = matrix(nodes)
        lines = valid.splitlines()
        bad = ["\n".join(lines[:-1]) + "\n", valid + lines[1] + "\n",
               valid.replace("n_A\tn_B", "unknown\tn_B", 1),
               valid.replace("\t0.0\t", "\tnan\t", 1),
               valid.replace("\t0.0\t", "\t101\t", 1),
               valid.replace("\t0.0\t", "\t5\t", 1),
               valid.replace("\tA\tB\t", "\twrong\tB\t", 1),
               valid.replace("\t100000000\t100000000", "\t10\t100000000", 1)]
        for text in bad:
            with self.subTest(text=text), self.assertRaises(ValueError):
                GROUP.group_phages(text, nodes, ["A"], [])
        for models, seeds in [(["unknown"], []), (["A", "A"], []), (["A"], ["unknown"])]:
            with self.assertRaises(ValueError):
                GROUP.group_phages(valid, nodes, models, seeds)

    def test_policy_input_checksum_rejects_changed_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "input"
            path.write_bytes(b"changed")
            with self.assertRaises(ValueError):
                GROUP.checked(path, "0" * 64)


if __name__ == "__main__":
    unittest.main()
