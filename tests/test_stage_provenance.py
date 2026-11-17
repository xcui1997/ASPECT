"""Independent stages reject stale or missing upstream classification inputs."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import pandas as pd

from aspect.pipeline import _load_test_sequences, _migration_test_targets, _write_classification_provenance


class StageProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.text_dir = Path(self.temp.name)
        self.cfg = {'migration_test_params': {'analyze_type': 'All'},
                    'classification_params': {'threshold': 1},
                    'analysis_settings': {}, 'nna_params': {}}
        (self.text_dir / 'refined_manifest.json').write_text('["refined_example.csv"]')
        (self.text_dir / 'refined_example.csv').write_text('refined_cluster,mag\n0,1\n')
        pd.DataFrame([{'file_name': 'refined_example.csv', 'refined_cluster': 0,
                       'sequence_type': 'Swarm'}]).to_csv(
                           self.text_dir / 'classification_results.csv', index=False)
        _write_classification_provenance(self.cfg, self.text_dir)

    def tearDown(self):
        self.temp.cleanup()

    def test_unmodified_inputs_are_accepted(self):
        self.assertEqual(len(_migration_test_targets(self.cfg, self.text_dir)), 1)

    def test_changed_refined_data_is_rejected(self):
        (self.text_dir / 'refined_example.csv').write_text('refined_cluster,mag\n0,2\n')
        with self.assertRaisesRegex(ValueError, 'Classification is stale'):
            _migration_test_targets(self.cfg, self.text_dir)

    def test_missing_refined_data_is_not_silently_skipped(self):
        (self.text_dir / 'refined_example.csv').unlink()
        with self.assertRaisesRegex(FileNotFoundError, 'Refined file missing'):
            _migration_test_targets(self.cfg, self.text_dir)
        targets = pd.DataFrame([{'file_name': 'refined_example.csv', 'refined_cluster': 0}])
        with self.assertRaisesRegex(FileNotFoundError, 'Refined file missing'):
            list(_load_test_sequences(self.cfg, self.text_dir, targets))

    def test_absent_cluster_is_not_silently_skipped(self):
        pd.DataFrame({'refined_cluster': [0], 'stime': [0.], 'x': [0.],
                      'y': [0.], 'dep': [1.], 'mag': [1.]}).to_csv(
                          self.text_dir / 'refined_example.csv', index=False)
        targets = pd.DataFrame([{'file_name': 'refined_example.csv', 'refined_cluster': 9}])
        with self.assertRaisesRegex(ValueError, 'Cluster 9 is missing'):
            list(_load_test_sequences(self.cfg, self.text_dir, targets))

    def test_changed_classification_settings_are_rejected(self):
        self.cfg['classification_params']['threshold'] = 2
        with self.assertRaisesRegex(ValueError, 'Classification is stale'):
            _migration_test_targets(self.cfg, self.text_dir)


if __name__ == '__main__':
    unittest.main()
