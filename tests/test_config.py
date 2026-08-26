import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import (
    AlignmentConfig,
    ClusteringConfig,
    ConservationConfig,
    InclusivityConfig,
    NcbiInputConfig,
    OligoConstraints,
    PipelineConfig,
    PrimerDesignConfig,
    load_config,
    validate_inclusivity_config,
    validate_primer_design_config,
)


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class PipelineConfigTests(unittest.TestCase):
    def _load_yaml(self, text: str) -> PipelineConfig:
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(text, encoding="utf-8")
            return load_config(config_path)

    def test_loads_ncbi_query_configuration(self):
        config = self._load_yaml(
            "target:\n"
            "  name: query-target\n"
            "input:\n"
            "  ncbi:\n"
            "    query: example[Organism]\n"
            "    batch_size: 25\n"
            "    retries: 2\n"
            "    max_records: 50\n"
        )

        self.assertEqual(config.input_ncbi.query, "example[Organism]")
        self.assertEqual(config.input_ncbi.accessions, ())
        self.assertIsNone(config.input_ncbi.frozen_dataset)
        self.assertEqual(config.input_ncbi.batch_size, 25)
        self.assertEqual(config.input_ncbi.retries, 2)
        self.assertEqual(config.input_ncbi.max_records, 50)
        self.assertEqual(config.selected_input, config.input_ncbi)

    def test_loads_ncbi_accessions_and_frozen_modes(self):
        accessions = self._load_yaml(
            "target:\n  name: accessions\ninput:\n  ncbi:\n"
            "    accessions: [NC_000001.11, AB123456.2]\n"
        )
        frozen = self._load_yaml(
            "target:\n  name: frozen\ninput:\n  ncbi:\n"
            "    frozen_dataset: datasets/frozen\n"
        )

        self.assertEqual(accessions.input_ncbi.accessions, ("NC_000001.11", "AB123456.2"))
        self.assertEqual(accessions.input_ncbi.batch_size, 100)
        self.assertEqual(accessions.input_ncbi.retries, 3)
        self.assertEqual(frozen.input_ncbi.frozen_dataset, Path("datasets/frozen"))

    def test_rejects_invalid_ncbi_configurations(self):
        invalid_cases = (
            (
                "multiple top-level sources",
                "input:\n  fasta: target.fasta\n  ncbi:\n    query: example[Organism]\n",
                "Exactly one local sequence input",
            ),
            (
                "no NCBI mode",
                "input:\n  ncbi: {}\n",
                "exactly one of query, accessions, or frozen_dataset",
            ),
            (
                "multiple NCBI modes",
                "input:\n  ncbi:\n    query: example[Organism]\n    accessions: [NC_000001.11]\n",
                "exactly one of query, accessions, or frozen_dataset",
            ),
            (
                "empty accessions",
                "input:\n  ncbi:\n    accessions: []\n",
                "input.ncbi.accessions",
            ),
            (
                "duplicate accessions",
                "input:\n  ncbi:\n    accessions: [NC_000001.11, NC_000001.11]\n",
                "input.ncbi.accessions",
            ),
            (
                "non-string accession",
                "input:\n  ncbi:\n    accessions: [NC_000001.11, 7]\n",
                "input.ncbi.accessions",
            ),
            (
                "max records without query",
                "input:\n  ncbi:\n    accessions: [NC_000001.11]\n    max_records: 5\n",
                "input.ncbi.max_records.*query",
            ),
            (
                "frozen mode batch size",
                "input:\n  ncbi:\n    frozen_dataset: datasets/frozen\n    batch_size: 2\n",
                "input.ncbi.batch_size.*frozen_dataset",
            ),
            (
                "frozen mode retries",
                "input:\n  ncbi:\n    frozen_dataset: datasets/frozen\n    retries: 2\n",
                "input.ncbi.retries.*frozen_dataset",
            ),
            (
                "frozen mode max records",
                "input:\n  ncbi:\n    frozen_dataset: datasets/frozen\n    max_records: 2\n",
                "input.ncbi.max_records.*frozen_dataset",
            ),
            (
                "batch size too low",
                "input:\n  ncbi:\n    query: example[Organism]\n    batch_size: 0\n",
                "input.ncbi.batch_size",
            ),
            (
                "batch size too high",
                "input:\n  ncbi:\n    query: example[Organism]\n    batch_size: 501\n",
                "input.ncbi.batch_size",
            ),
            (
                "null batch size",
                "input:\n  ncbi:\n    query: example[Organism]\n    batch_size: null\n",
                "input.ncbi.batch_size.*integer",
            ),
            (
                "retries too low",
                "input:\n  ncbi:\n    query: example[Organism]\n    retries: -1\n",
                "input.ncbi.retries",
            ),
            (
                "retries too high",
                "input:\n  ncbi:\n    query: example[Organism]\n    retries: 11\n",
                "input.ncbi.retries",
            ),
            (
                "null retries",
                "input:\n  ncbi:\n    query: example[Organism]\n    retries: null\n",
                "input.ncbi.retries.*integer",
            ),
            (
                "max records not positive",
                "input:\n  ncbi:\n    query: example[Organism]\n    max_records: 0\n",
                "input.ncbi.max_records",
            ),
        )
        for name, input_yaml, message in invalid_cases:
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_yaml("target:\n  name: invalid\n" + input_yaml)

    def test_rejects_unknown_ncbi_configuration_keys(self):
        for key in ("api_key", "credential", "arbitrary"):
            with self.subTest(key=key), self.assertRaisesRegex(
                ValueError, rf"input\.ncbi.*{key}.*unrecognized"
            ):
                self._load_yaml(
                    "target:\n"
                    "  name: invalid\n"
                    "input:\n"
                    "  ncbi:\n"
                    "    query: example[Organism]\n"
                    f"    {key}: secret-value\n"
                )

    def test_selected_input_rejects_ambiguous_direct_pipeline_sources(self):
        cases = (
            PipelineConfig(
                target_name="ambiguous",
                input_fasta=Path("target.fasta"),
                input_genbank=Path("target.gb"),
            ),
            PipelineConfig(
                target_name="ambiguous",
                input_fasta=Path("target.fasta"),
                input_ncbi=NcbiInputConfig(query="example[Organism]"),
            ),
        )

        for config in cases:
            with self.subTest(config=config), self.assertRaisesRegex(
                ValueError, "Exactly one sequence input"
            ):
                _ = config.selected_input

    def test_selected_input_rejects_invalid_direct_ncbi_modes_and_fields(self):
        cases = (
            (NcbiInputConfig(query="virus", accessions=("NC_1",)), "exactly one"),
            (
                NcbiInputConfig(query="virus", frozen_dataset=Path("frozen")),
                "exactly one",
            ),
            (
                NcbiInputConfig(accessions=("NC_1",), frozen_dataset=Path("frozen")),
                "exactly one",
            ),
            (NcbiInputConfig(accessions=("NC_1", "NC_1")), "unique"),
            (NcbiInputConfig(accessions=(" ",)), "non-blank"),
            (NcbiInputConfig(query="virus", batch_size=True), "batch_size"),
            (NcbiInputConfig(query="virus", retries=1.5), "retries"),
            (NcbiInputConfig(accessions=("NC_1",), max_records=1), "max_records"),
            (
                NcbiInputConfig(frozen_dataset=Path("frozen"), batch_size=1),
                "frozen_dataset.*batch_size",
            ),
        )

        for ncbi_config, error in cases:
            with self.subTest(ncbi_config=ncbi_config), self.assertRaisesRegex(
                ValueError, error
            ):
                _ = PipelineConfig(
                    target_name="invalid", input_ncbi=ncbi_config
                ).selected_input

    def test_loads_minimal_yaml_configuration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n",
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIsInstance(config, PipelineConfig)
        self.assertEqual(config.target_name, "synthetic-target")
        self.assertEqual(config.input_fasta, Path("tests/fixtures/target_small.fasta"))

    def test_loads_genbank_input_and_optional_qc_thresholds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  genbank: tests/fixtures/target.gb\n"
                "qc:\n"
                "  min_length: 100\n"
                "  max_ambiguous_fraction: 0.05\n"
                "  expected_length: 150\n"
                "  length_tolerance_fraction: 0.10\n",
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIsNone(config.input_fasta)
        self.assertEqual(config.input_genbank, Path("tests/fixtures/target.gb"))
        self.assertEqual(config.selected_input, (Path("tests/fixtures/target.gb"), "genbank"))
        self.assertEqual(config.qc.min_length, 100)
        self.assertEqual(config.qc.max_ambiguous_fraction, 0.05)
        self.assertEqual(config.qc.expected_length, 150)
        self.assertEqual(config.qc.length_tolerance_fraction, 0.10)

    def test_rejects_fractional_min_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "qc:\n"
                "  min_length: 100.5\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "qc.min_length.*integer"):
                load_config(config_path)

    def test_rejects_fractional_expected_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "qc:\n"
                "  expected_length: 150.5\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "qc.expected_length.*integer"):
                load_config(config_path)

    def test_requires_exactly_one_local_sequence_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "  genbank: tests/fixtures/target.gb\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Exactly one local sequence input"):
                load_config(config_path)

    def test_loads_clustering_configuration_and_defaults_when_omitted(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.clustering, ClusteringConfig())

        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "clustering:\n"
            "  enabled: true\n"
            "  identity: 0.9\n"
            "  threads: 4\n"
            "  memory_mb: 2048\n"
        )
        self.assertEqual(
            config.clustering,
            ClusteringConfig(enabled=True, identity=0.9, threads=4, memory_mb=2048),
        )

    def test_loads_alignment_configuration_and_defaults_when_omitted(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.alignment, AlignmentConfig())

        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n"
            "  enabled: true\n"
            "  threads: 4\n"
            "  reference_id: seq-1\n"
        )
        self.assertEqual(
            config.alignment,
            AlignmentConfig(enabled=True, threads=4, reference_id="seq-1"),
        )

    def test_loads_conservation_configuration_and_defaults_when_omitted(self):
        minimal_yaml = (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
        )
        config = self._load_yaml(
            minimal_yaml
            + "conservation:\n"
            + "  enabled: true\n"
            + "  window_size: 120\n"
            + "  step_size: 12\n"
        )
        self.assertEqual(
            config.conservation,
            ConservationConfig(enabled=True, window_size=120, step_size=12),
        )
        self.assertEqual(self._load_yaml(minimal_yaml).conservation, ConservationConfig())

    def test_loads_primer_design_configuration_and_defaults_when_omitted(self):
        minimal_yaml = (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
        )
        config = self._load_yaml(minimal_yaml)

        self.assertEqual(config.primer_design.max_candidate_regions, 10)
        self.assertEqual(config.primer_design.primer.opt_tm, 60.0)
        self.assertEqual(config.primer_design.probe.opt_tm, 70.0)

        loaded = self._load_yaml(
            minimal_yaml
            + "primer_design:\n"
            "  enabled: true\n"
            "  max_candidate_regions: 3\n"
            "  assays_per_region: 4\n"
            "  candidate_region_length: 220\n"
            "  product_size_max: 180\n"
            "  primer:\n"
            "    opt_tm: 59.5\n"
            "  probe:\n"
            "    min_size: 20\n"
        )

        self.assertTrue(loaded.primer_design.enabled)
        self.assertEqual(loaded.primer_design.max_candidate_regions, 3)
        self.assertEqual(loaded.primer_design.primer.opt_tm, 59.5)
        self.assertEqual(loaded.primer_design.probe.min_size, 20)

    def test_loads_inclusivity_configuration_and_defaults_when_omitted(self):
        base = (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
            "primer_design:\n  enabled: true\n"
        )
        self.assertEqual(self._load_yaml(base).inclusivity, InclusivityConfig())

        loaded = self._load_yaml(
            base
            + "inclusivity:\n"
            "  enabled: true\n"
            "  search_flank: 40\n"
            "  max_hits_per_oligo: 7\n"
            "  max_primer_mismatches: 1\n"
            "  max_probe_mismatches: 0\n"
            "  reject_primer_3_prime_mismatch: false\n"
            "  primer_3_prime_bases: 4\n"
            "  max_primer_degeneracy: 8\n"
            "  max_probe_degeneracy: 2\n"
            "  allow_primer_3_prime_degeneracy: true\n"
            "  max_amplicon_size_delta: 5\n"
        )
        self.assertEqual(
            loaded.inclusivity,
            InclusivityConfig(
                enabled=True,
                search_flank=40,
                max_hits_per_oligo=7,
                max_primer_mismatches=1,
                max_probe_mismatches=0,
                reject_primer_3_prime_mismatch=False,
                primer_3_prime_bases=4,
                max_primer_degeneracy=8,
                max_probe_degeneracy=2,
                allow_primer_3_prime_degeneracy=True,
                max_amplicon_size_delta=5,
            ),
        )

    def test_rejects_invalid_inclusivity_yaml_configuration(self):
        base = (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
            "primer_design:\n  enabled: true\n"
        )
        invalid_cases = (
            ("non-mapping", "inclusivity: true\n", "inclusivity.*mapping"),
            ("unknown", "inclusivity:\n  extra: 1\n", "inclusivity.*extra.*unrecognized"),
            ("integer boolean", "inclusivity:\n  enabled: 1\n", "enabled.*boolean"),
            ("negative zero allowed", "inclusivity:\n  search_flank: -1\n", "search_flank.*non-negative integer"),
            ("zero positive only", "inclusivity:\n  max_hits_per_oligo: 0\n", "max_hits_per_oligo.*positive integer"),
            ("fractional integer", "inclusivity:\n  max_probe_degeneracy: 1.5\n", "max_probe_degeneracy.*positive integer"),
        )
        for name, section, message in invalid_cases:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                self._load_yaml(base + section)

    def test_rejects_invalid_direct_inclusivity_configuration(self):
        invalid = (
            InclusivityConfig(enabled=1),
            InclusivityConfig(search_flank=-1),
            InclusivityConfig(max_hits_per_oligo=0),
            InclusivityConfig(max_probe_degeneracy=1.5),
        )
        for config in invalid:
            with self.subTest(config=config), self.assertRaises(ValueError):
                validate_inclusivity_config(config)

        with self.assertRaisesRegex(ValueError, "inclusivity.*requires enabled primer design"):
            _ = PipelineConfig(
                target_name="target",
                input_fasta=FIXTURE_FASTA,
                inclusivity=InclusivityConfig(enabled=True),
            ).selected_input

    def test_rejects_invalid_primer_design_yaml_configuration(self):
        invalid_cases = (
            ("non-mapping", "primer_design: true\n", "primer_design.*mapping"),
            ("unknown key", "primer_design:\n  extra: true\n", "primer_design.*extra.*unrecognized"),
            ("enabled integer", "primer_design:\n  enabled: 1\n", "enabled.*boolean"),
            ("candidate count bool", "primer_design:\n  max_candidate_regions: true\n", "max_candidate_regions.*positive integer"),
            ("candidate count zero", "primer_design:\n  max_candidate_regions: 0\n", "max_candidate_regions.*positive integer"),
            ("fraction non-finite", "primer_design:\n  min_mean_conservation: .nan\n", "min_mean_conservation.*finite.*0.*1"),
            ("fraction high", "primer_design:\n  min_usable_fraction: 1.01\n", "min_usable_fraction.*0.*1"),
            ("unknown primer field", "primer_design:\n  primer:\n    extra: 1\n", "primer.*extra.*unrecognized"),
            ("primer size bool", "primer_design:\n  primer:\n    min_size: true\n", "primer.*min_size.*positive integer"),
            ("primer unordered sizes", "primer_design:\n  primer:\n    min_size: 21\n    opt_size: 20\n", "primer.*size"),
            ("probe non-finite tm", "primer_design:\n  probe:\n    opt_tm: .inf\n", "probe.*opt_tm.*finite"),
            ("probe unordered tm", "primer_design:\n  probe:\n    min_tm: 71\n    opt_tm: 70\n", "probe.*Tm"),
            ("primer gc low", "primer_design:\n  primer:\n    min_gc_percent: -0.1\n", "primer.*min_gc_percent.*0.*100"),
            ("product range", "primer_design:\n  product_size_min: 201\n", "product_size_min"),
            ("candidate region short", "primer_design:\n  candidate_region_length: 150\n  product_size_max: 200\n", "candidate_region_length"),
            ("enabled without conservation", "primer_design:\n  enabled: true\n", "requires enabled conservation"),
        )
        base_yaml = (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
        )
        for name, primer_design_yaml, message in invalid_cases:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                yaml = base_yaml
                if name == "enabled without conservation":
                    yaml = base_yaml.replace("conservation:\n  enabled: true\n", "")
                self._load_yaml(yaml + primer_design_yaml)

    def test_rejects_invalid_direct_primer_design_configuration(self):
        invalid = (
            PrimerDesignConfig(max_candidate_regions=True),
            PrimerDesignConfig(max_region_overlap_fraction=float("nan")),
            PrimerDesignConfig(primer=OligoConstraints(21, 20, 25, 58.0, 60.0, 62.0, 40.0, 60.0)),
            PrimerDesignConfig(probe=OligoConstraints(18, 25, 30, 71.0, 70.0, 72.0, 30.0, 80.0)),
            PrimerDesignConfig(candidate_region_length=150, product_size_max=200),
        )
        for primer_design in invalid:
            with self.subTest(primer_design=primer_design), self.assertRaises(ValueError):
                validate_primer_design_config(primer_design)

    def test_accepts_two_bit_primer_design_entropy_limit(self):
        validate_primer_design_config(PrimerDesignConfig(max_mean_entropy_bits=2.0))

    def test_rejects_invalid_conservation_yaml_configuration(self):
        invalid_cases = (
            ("non-mapping", "conservation: true\n", "conservation.*mapping"),
            (
                "unknown key",
                "conservation:\n  extra: true\n",
                "conservation.*extra.*unrecognized",
            ),
            (
                "enabled integer",
                "conservation:\n  enabled: 1\n",
                "enabled.*boolean",
            ),
            (
                "window bool",
                "conservation:\n  window_size: true\n",
                "window_size.*integer",
            ),
            (
                "window fractional",
                "conservation:\n  window_size: 1.5\n",
                "window_size.*integer",
            ),
            (
                "window zero",
                "conservation:\n  window_size: 0\n",
                "window_size.*1.*1000000",
            ),
            (
                "step high",
                "conservation:\n  step_size: 1000001\n",
                "step_size.*1.*1000000",
            ),
            (
                "step greater than window",
                "conservation:\n  window_size: 10\n  step_size: 11\n",
                "step_size.*cannot exceed.*window_size",
            ),
        )
        for name, conservation_yaml, message in invalid_cases:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                self._load_yaml(
                    "target:\n  name: target\n"
                    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
                    + conservation_yaml
                )

    def test_selected_input_rejects_invalid_direct_conservation_configuration(self):
        invalid = (
            ConservationConfig(enabled=1),
            ConservationConfig(window_size=True),
            ConservationConfig(window_size=1.5),
            ConservationConfig(window_size=0),
            ConservationConfig(step_size=1_000_001),
            ConservationConfig(window_size=10, step_size=11),
        )
        for conservation in invalid:
            with self.subTest(conservation=conservation):
                config = PipelineConfig(
                    target_name="target",
                    input_fasta=FIXTURE_FASTA,
                    conservation=conservation,
                )
                with self.assertRaises(ValueError):
                    _ = config.selected_input

        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            alignment=AlignmentConfig(enabled=False),
            conservation=ConservationConfig(enabled=True),
        )
        with self.assertRaisesRegex(ValueError, "requires enabled alignment"):
            _ = config.selected_input

    def test_rejects_invalid_alignment_yaml_configuration(self):
        invalid_cases = (
            ("unknown key", "  extra: true\n", "alignment.*unrecognized"),
            ("enabled integer", "  enabled: 1\n", "enabled.*boolean"),
            ("threads bool", "  threads: true\n", "threads.*integer"),
            ("threads fractional", "  threads: 1.5\n", "threads.*integer"),
            ("threads zero", "  threads: 0\n", "threads.*1.*256"),
            ("threads high", "  threads: 257\n", "threads.*1.*256"),
            ("numeric reference", "  reference_id: 7\n", "reference_id.*non-blank string"),
            ("empty reference", "  reference_id: \"\"\n", "reference_id.*non-blank string"),
            ("whitespace reference", "  reference_id: \" \"\n", "reference_id.*non-blank string"),
        )
        for name, alignment_yaml, message in invalid_cases:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                self._load_yaml(
                    "target:\n  name: target\n"
                    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
                    "alignment:\n" + alignment_yaml
                )

    def test_selected_input_rejects_invalid_direct_alignment_configuration(self):
        for alignment in (
            AlignmentConfig(enabled=1),
            AlignmentConfig(threads=True),
            AlignmentConfig(threads=0),
            AlignmentConfig(reference_id=7),
            AlignmentConfig(reference_id=" "),
        ):
            with self.subTest(alignment=alignment):
                config = PipelineConfig(
                    target_name="target",
                    input_fasta=FIXTURE_FASTA,
                    alignment=alignment,
                )
                with self.assertRaises(ValueError):
                    _ = config.selected_input

    def test_rejects_invalid_clustering_yaml_configuration(self):
        invalid_cases = (
            ("unknown key", "  extra: true\n", "clustering.*unrecognized"),
            ("enabled", "  enabled: 1\n", "enabled.*boolean"),
            ("identity bool", "  identity: true\n", "identity.*number"),
            ("identity string", "  identity: nope\n", "identity.*number"),
            ("identity low", "  identity: 0.799\n", "identity.*0.80.*1.0"),
            ("identity high", "  identity: 1.01\n", "identity.*0.80.*1.0"),
            ("threads bool", "  threads: true\n", "threads.*integer"),
            ("threads fractional", "  threads: 1.5\n", "threads.*integer"),
            ("threads low", "  threads: 0\n", "threads.*1.*256"),
            ("threads high", "  threads: 257\n", "threads.*1.*256"),
            ("memory bool", "  memory_mb: true\n", "memory_mb.*positive integer"),
            ("memory fractional", "  memory_mb: 1.5\n", "memory_mb.*positive integer"),
            ("memory non-positive", "  memory_mb: 0\n", "memory_mb.*positive integer"),
        )
        for name, clustering_yaml, message in invalid_cases:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                self._load_yaml(
                    "target:\n  name: target\n"
                    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
                    "clustering:\n" + clustering_yaml
                )

    def test_selected_input_rejects_invalid_direct_clustering_configuration(self):
        for clustering in (
            ClusteringConfig(enabled=1),
            ClusteringConfig(identity=True),
            ClusteringConfig(identity=0.799),
            ClusteringConfig(identity=1.01),
            ClusteringConfig(threads=0),
            ClusteringConfig(memory_mb=0),
        ):
            with self.subTest(clustering=clustering):
                config = PipelineConfig(
                    target_name="target",
                    input_fasta=FIXTURE_FASTA,
                    clustering=clustering,
                )
                with self.assertRaises(ValueError):
                    _ = config.selected_input


if __name__ == "__main__":
    unittest.main()
