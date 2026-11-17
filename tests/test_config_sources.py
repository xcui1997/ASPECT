"""Configuration overrides must not retain the bundled example catalog."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from aspect.config import load_config


class ConfigSourceTests(unittest.TestCase):
    def test_source_mapping_is_replaced(self):
        with TemporaryDirectory() as directory:
            config_path = Path(directory) / 'config.yaml'
            config_path.write_text('data_sources:\n  MY_CATALOG: input.csv\n')
            cfg = load_config(config_path)
            self.assertEqual(set(cfg['data_sources']), {'MY_CATALOG'})
            self.assertEqual(cfg['data_sources']['MY_CATALOG'], str(Path(directory).resolve() / 'input.csv'))


if __name__ == '__main__':
    unittest.main()
