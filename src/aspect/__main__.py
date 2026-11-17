"""Run the ASPECT workflow using a YAML configuration."""
from pathlib import Path
import argparse
from .pipeline import load_and_run
from .config import DEFAULT_CONFIG_PATH


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG_PATH,
                        help='YAML configuration path; defaults to the project config.yaml.')
    parser.add_argument('--steps', default=None,
                        help='Comma-separated workflow stages, e.g. nna,dbscan,classify,migration_test.')
    args = parser.parse_args()
    steps = [item.strip() for item in args.steps.split(',') if item.strip()] if args.steps else None
    load_and_run(args.config.resolve(), steps=steps)


if __name__ == '__main__':
    main()
