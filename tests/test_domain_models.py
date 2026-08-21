import unittest

from qpcr_pipeline.models import DiscoverySet, EvaluationSet, TargetSequenceSet


class DomainSetContractTests(unittest.TestCase):
    def test_discovery_and_evaluation_sets_are_distinct_views_of_target(self):
        target = TargetSequenceSet(sequence_ids=("seq-1", "seq-2", "seq-3"))
        discovery = DiscoverySet(sequence_ids=("seq-1", "seq-3"))
        evaluation = EvaluationSet(sequence_ids=target.sequence_ids)

        self.assertEqual(target.sequence_ids, ("seq-1", "seq-2", "seq-3"))
        self.assertEqual(discovery.sequence_ids, ("seq-1", "seq-3"))
        self.assertEqual(evaluation.sequence_ids, target.sequence_ids)
        self.assertIsNot(type(discovery), type(evaluation))
        self.assertLess(len(discovery.sequence_ids), len(evaluation.sequence_ids))


if __name__ == "__main__":
    unittest.main()
