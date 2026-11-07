import importlib.util
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location("host_groups", Path(__file__).resolve().parents[3] / "scripts/group-genophi-hosts.py")
GROUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GROUP)


def matrix(names, values, row_order=None):
    order = row_order or names
    return "bacteria\t" + "\t".join(names) + "\n" + "".join(
        name + "\t" + "\t".join(str(values[name][other]) for other in names) + "\n" for name in order)


class HostGroupTests(unittest.TestCase):
    def setUp(self):
        self.names = ["A-1", "bridge", "C", "D"]
        self.values = {a: {b: 0 if a == b else 0.1 for b in self.names} for a in self.names}
        for a, b, value in [("A-1", "bridge", 0.00009), ("bridge", "C", 0.00009), ("C", "D", 0.0001)]:
            self.values[a][b] = self.values[b][a] = value

    def test_transitive_components_include_nontraining_bridges_and_strict_threshold(self):
        result = GROUP.group_hosts(matrix(self.names, self.values), ["A1", "C", "D"])
        groups = result["training_host_groups"]
        self.assertEqual(groups["A1"], groups["C"])
        self.assertNotEqual(groups["C"], groups["D"])
        self.assertEqual(result["groups"][groups["A1"]]["source_members"], ["A-1", "C", "bridge"])
        self.assertEqual(result["groups"][groups["A1"]]["training_members"], ["A1", "C"])

    def test_input_order_does_not_change_groups(self):
        first = GROUP.group_hosts(matrix(self.names, self.values), ["C", "A1", "D"])
        second = GROUP.group_hosts(matrix(self.names[::-1], self.values, ["C", "D", "bridge", "A-1"]), ["D", "C", "A1"])
        self.assertEqual(first, second)

    def test_corrupt_matrices_rejected(self):
        original = matrix(self.names, self.values)
        cases = [original.replace("0.0001", "nan", 1), original.replace("0.0001", "-0.1", 1),
                 original.replace("0.0001", "0.2", 1), original.replace("A-1\t0\t", "A-1\t1\t"),
                 "\n".join(original.splitlines()[:-1]), original.replace("A-1\t0\t", "A-1\t\t")]
        for value in cases:
            with self.subTest(value=value), self.assertRaises(ValueError):
                GROUP.group_hosts(value, ["A1"])

    def test_name_collisions_and_unmapped_training_ids_rejected(self):
        text = matrix(self.names, self.values)
        for hosts in [["A1", "A-1"], ["unknown"], [""]]:
            with self.subTest(hosts=hosts), self.assertRaises(ValueError):
                GROUP.group_hosts(text, hosts)
        with self.assertRaises(ValueError):
            GROUP.group_hosts("bacteria\tA-1\tA1\nA-1\t0\t0\nA1\t0\t0\n", ["A1"])

    def test_missing_target_names_are_not_claimed_as_sequence_exclusions(self):
        result = GROUP.group_hosts("bacteria\tBL21\nBL21\t0\n", ["BL21"])
        self.assertEqual(result["named_target_groups"]["BL21"]["status"], "provisional_name_match")
        self.assertIsNone(result["named_target_groups"]["MG1655"]["group_id"])
