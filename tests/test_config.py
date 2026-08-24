import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import ClusteringConfig, NcbiInputConfig, PipelineConfig, load_config


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

    def test_rejects_invalid_clustering_yaml_configuration(self):
        invalid_cases = (
            ("unknown key", "  extra: true\n", "clustering.*unrecognized"),
            ("enabled", "  enabled: 1\n", "enabled.*boolean"),
            ("identity bool", "  identity: true\n", "identity.*number"),
            ("identity string", "  identity: nope\n", "identity.*number"),
            ("identity low", "  identity: 0.74\n", "identity.*0.75.*1.0"),
            ("identity high", "  identity: 1.01\n", "identity.*0.75.*1.0"),
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
            ClusteringConfig(identity=0.74),
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
