"""The 50-event gate applies before classification."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from aspect.pipeline import _classification_minimum_size, run_classify

class _FakeClassifier:
    def __init__(self, cfg): pass
    def calculate_metrics(self, cluster): return {'count': len(cluster)}
    def classify(self, frame): return frame
    def plot_diagnostic(self, results, path): pass

class ClassificationGateTest(unittest.TestCase):
    def test_fixed_and_fractional_requirements(self):
        cfg = {'analysis_settings': {'min_events_per_cluster': 50}, 'classification_params': {'min_parent_fraction': .1}}
        self.assertEqual(_classification_minimum_size(cfg, 400), 50)
        self.assertTrue(49 < _classification_minimum_size(cfg, 400) <= 50)
        self.assertEqual(_classification_minimum_size(cfg, 800), 80)

    def test_short_child_never_reaches_classification(self):
        cfg = {'analysis_settings': {'min_events_per_cluster': 50}, 'classification_params': {'min_parent_fraction': .1}, 'dbscan_params': {'min_cluster_size': 30}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            # Both are valid DBSCAN outputs, but only 50 meets the paper gate.
            df = pd.DataFrame({'refined_cluster': [0]*49 + [1]*50, 'mag': [2.]*99, 'lat': [35.]*99, 'lon': [-117.]*99})
            df.to_csv(path/'refined_test.csv', index=False)
            # A stale output file is ignored when absent from the current manifest.
            df.to_csv(path/'refined_stale.csv', index=False)
            (path/'refined_manifest.json').write_text('["refined_test.csv"]')
            with patch('aspect.pipeline.ensure_output_dirs', return_value=(path, path, path)), patch('aspect.pipeline.ClusterClassifier', _FakeClassifier), patch('aspect.pipeline._prepare_sequence_dataframe', side_effect=lambda df, *a, **k: df):
                run_classify(cfg)
            result = pd.read_csv(path/'classification_results.csv')
            self.assertEqual(result['refined_cluster'].tolist(), [1])
            self.assertEqual(result['count'].tolist(), [50])

if __name__ == '__main__': unittest.main()
