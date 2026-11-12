"""Regression tests for HW1 Exercise 1 (2026/27). Run: python -m unittest discover hw1/tests"""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw1_problem1", ROOT / "src" / "run_problem1.py")
M = importlib.util.module_from_spec(SPEC)
import sys
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


class TestHW1_2026_27(unittest.TestCase):
    def test_base_network(self):
        self.assertEqual([e.capacity for e in M.BASE_EDGES], [2, 2, 2, 2, 3, 1])
        self.assertEqual(M.max_flow_value(M.BASE_EDGES), 4)
        cuts = M.cut_table(M.BASE_EDGES, M.NODES)
        self.assertEqual(len(cuts), 8)
        self.assertEqual(cuts.capacity.min(), 4)

    def test_legal_candidate_links(self):
        links = M.candidate_added_links(M.NODES, M.BASE_EDGES)
        self.assertEqual(len(links), 13)
        self.assertNotIn(("o", "d"), links)
        self.assertIn(("o", "c"), links)

    def test_closed_form_and_construction(self):
        for x in range(21):
            k, r = divmod(x, 2)
            existing = [
                M.Edge(e.label, e.tail, e.head,
                       e.capacity + (k if e.label in {"e1", "e2"} else r if e.label=="e6" else 0))
                for e in M.BASE_EDGES
            ]
            self.assertEqual(M.max_flow_value(existing), 4+k)
            with_link = [
                M.Edge(e.label, e.tail, e.head,
                       e.capacity + (k if e.label=="e5" else r if e.label=="e6" else 0))
                for e in M.BASE_EDGES
            ] + [M.Edge("e7", "o", "c", 1+k)]
            self.assertEqual(M.max_flow_value(with_link), 5+k)


if __name__ == "__main__":
    unittest.main()
