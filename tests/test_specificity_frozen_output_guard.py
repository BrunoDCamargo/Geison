import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import OffTargetConfig, SpecificityConfig
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.specificity import SpecificityError, evaluate_specificity


class SpecificityFrozenOutputGuardTests(unittest.TestCase):
    def test_rejects_output_inside_frozen_off_target_before_writing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            frozen = Path(tmpdir) / "frozen"
            frozen.mkdir()
            primer_design = PrimerDesignResult(
                status="COMPLETE",
                reference_id="ref",
                candidates=(),
                assays=(),
                candidate_regions_path=None,
                assays_path=None,
                primer3_input_path=None,
                primer3_output_path=None,
                report_path=Path("unused-primer-report.json"),
            )

            with self.assertRaisesRegex(
                SpecificityError,
                "output.*frozen off-target",
            ):
                evaluate_specificity(
                    primer_design,
                    (OffTargetConfig(name="neighbors", frozen_dataset=frozen),),
                    SpecificityConfig(enabled=True),
                    frozen,
                )

            self.assertFalse((frozen / "specificity").exists())


if __name__ == "__main__":
    unittest.main()
