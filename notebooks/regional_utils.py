"""Shared, auditable stage execution for the San Jacinto regional notebooks."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


def project_root():
    here = Path.cwd().resolve()
    for path in (here, *here.parents):
        if (path / 'run_aspect.py').is_file() and (path / 'src/aspect').is_dir():
            return path
    raise RuntimeError('Open this notebook from within the ASPECT repository.')


ROOT = project_root()
sys.path.insert(0, str(ROOT / 'src'))
from aspect.config import load_config  # noqa: E402

CONFIG_PATH = ROOT / 'examples/config_san_jacinto_region.yaml'
if not CONFIG_PATH.is_file():
    raise FileNotFoundError(f'Regional config not found: {CONFIG_PATH}')
CFG = load_config(CONFIG_PATH)
OUT = Path(CFG['output_dir'])
TEXT = OUT / 'text'
MIGRATION = OUT / 'migration_test'
SESSION = OUT / 'regional_notebook_session.json'
STAGES = ('nna', 'dbscan', 'classify', 'migration_test')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def input_fingerprint():
    paths = [CONFIG_PATH, ROOT / 'config.yaml', ROOT / 'run_aspect.py']
    paths += sorted((ROOT / 'src/aspect').glob('*.py'))
    paths += [Path(p) for p in CFG['data_sources'].values()]
    return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha256(p) for p in paths}


def product_paths(stage):
    if stage == 'nna':
        return sorted(TEXT / f'nna_{name}.csv' for name in CFG['data_sources'])
    if stage == 'dbscan':
        manifest = TEXT / 'refined_manifest.json'
        return [manifest] + [TEXT / name for name in json.loads(manifest.read_text())]
    if stage == 'classify':
        return [TEXT / 'classification_results.csv', TEXT / 'classification_provenance.json']
    if stage == 'migration_test':
        return [MIGRATION / 'migration_test_summary.csv', MIGRATION / 'migration_test_status.csv',
                MIGRATION / 'migration_test_runs.csv.gz']
    raise ValueError(f'Unknown stage: {stage}')


def stage_digest(stage):
    return {str(p.relative_to(OUT)): sha256(p) for p in product_paths(stage)}


def validated_session(previous_stage):
    if not SESSION.is_file():
        raise RuntimeError('Run 01_nna.ipynb before this notebook.')
    state = json.loads(SESSION.read_text())
    if state.get('input_fingerprint') != input_fingerprint():
        raise RuntimeError('Input catalog, configuration, or core code has changed. Rerun from 01_nna.ipynb.')
    index = STAGES.index(previous_stage)
    for stage in STAGES[:index + 1]:
        if stage not in state.get('stages', {}):
            raise RuntimeError(f'Stage {stage} is incomplete. Run the notebooks in order.')
        if state['stages'][stage] != stage_digest(stage):
            raise RuntimeError(f'Outputs for stage {stage} have changed. Rerun from that stage.')
    return state


def run_stage(stage):
    index = STAGES.index(stage)
    fingerprint = input_fingerprint()
    if index:
        state = validated_session(STAGES[index - 1])
    else:
        state = {'input_fingerprint': fingerprint, 'stages': {}}
    # Invalidate this stage and downstream stages before execution.
    for downstream in STAGES[index:]:
        state['stages'].pop(downstream, None)
    OUT.mkdir(parents=True, exist_ok=True)
    SESSION.write_text(json.dumps(state, indent=2) + '\n')
    command = [sys.executable, str(ROOT / 'run_aspect.py'), '--config', str(CONFIG_PATH), '--steps', stage]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    print(result.stdout)
    if result.returncode:
        print(result.stderr)
        raise RuntimeError(f'Stage {stage} failed with status {result.returncode}')
    if fingerprint != input_fingerprint():
        raise RuntimeError('Inputs changed during execution; rerun from NNA.')
    state['stages'][stage] = stage_digest(stage)
    SESSION.write_text(json.dumps(state, indent=2) + '\n')
    print(f'Completed {stage}; inputs and products recorded in {SESSION.name}.')


def refined_files():
    validated_session('dbscan')
    return [TEXT / name for name in json.loads((TEXT / 'refined_manifest.json').read_text())]
