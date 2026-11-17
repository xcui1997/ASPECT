<p align="center"><img src="assets/aspect_logo.png" alt="ASPECT logo" width="600"></p>

# ASPECT

[![Tests](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml/badge.svg)](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml) [![Python](https://img.shields.io/badge/python-%3E%3D3.10-blue)](pyproject.toml) [![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE) [![Developed with OpenAI Codex](https://img.shields.io/badge/Developed%20with-OpenAI%20Codex-087f8c)](https://openai.com/codex/)

**English** | [中文](README_zh.md)

[Source code](src/aspect/) | [Parameter reference](docs/parameters.md) | [San Jacinto example](notebooks/)

ASPECT is a Python workflow for identifying and classifying earthquake sequences and testing their spatiotemporal migration from an earthquake catalog. Its modular algorithms and YAML configuration support both complete runs and step-by-step exploration.

## Features

- Extract preliminary sequences with nearest-neighbor analysis (NNA), then refine them with spatial Area and temporal Knee clustering.
- Classify sequences using interpretable timing, magnitude and tree-structure features.
- Construct a rolling P70 distance front and compare growth models with non-migration references.
- Assess detection stability using the full sequence and 18 random 90% subsets; report individual runs and jointly passing velocity estimates.
- Run the workflow from a shell script or explore the San Jacinto example in a Jupyter notebook.

Spatial clustering uses projected horizontal coordinates; migration distances use three-dimensional coordinates including depth. Automatic UTM projection is intended for a geographically local catalog; broad regions, polar data, or data crossing the antimeridian require an appropriate projection strategy. This general workflow processes the supplied catalog without cross-catalog deduplication.

## Install

Python 3.10 or newer is required. Install the core workflow from the project directory:

```bash
python -m pip install -e .
```

For the regional notebooks and accelerated NNA, install the optional dependencies:

```bash
python -m pip install -e '.[notebook,acceleration]'
```

`numba` is strongly recommended for the regional example; runtime depends on the machine and catalog.

## San Jacinto regional example

The regional input includes **292,468 events across all available magnitudes** for 2008–2017 in the rectangular area 33.0–34.0°N and 117.4–116.4°W. The example NNA configuration applies M ≥ 0.3, retaining **72,061 events** for analysis; this is a magnitude threshold, not a selection of particular cases. The rectangle is a geographic selection, not an exact fault-zone boundary, and includes background and nearby-fault events. The source is a locally prepared QTM-derived catalog. See [data selection and provenance](examples/SAN_JACINTO_REGION.md).

Run the full regional workflow:

```bash
bash run_all.sh --config examples/config_san_jacinto_region.yaml
```

Or open the executed notebooks in order. Each notebook runs its stage and displays the complete stage results; later notebooks require outputs from earlier stages and check their provenance.

1. [Regional catalog](notebooks/00_regional_catalog.ipynb)
2. [NNA sequence extraction](notebooks/01_nna.ipynb)
3. [Spatial and temporal refinement](notebooks/02_refinement.ipynb)
4. [Sequence classification](notebooks/03_classification.ipynb)
5. [Migration detection](notebooks/04_migration.ipynb)

```bash
jupyter lab notebooks/
```

The regional configuration stores results under `results/san_jacinto_region/`. The original `config.yaml` remains a smaller 136-event smoke test; `bash run_all.sh` without arguments runs that smaller example. To use another catalog, copy a YAML, set `data_sources`, adjust parameters and run `bash run_all.sh --config /path/to/config.yaml`. Input requires `datetime`, `lat`, `lon`, `dep`, `mag`; time is UTC, coordinates are in degrees and depth in km. Paths are relative to the YAML file. See the [parameter reference](docs/parameters.md).

## Model comparison

The example sets the step and logarithmic reference margins to 0 and 4 IC units, respectively. The reported `dic` is `min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_growth`; therefore it is a margin-adjusted comparison, not an unadjusted IC difference against every reference. Both margins are configurable. A run jointly passes when `dic > 0` and the migration-phase R² exceeds 0.5 under the example configuration.

## Output

- `<output_dir>/text/nna_*.csv`: preliminary event assignments and nearest-neighbor information.
- `<output_dir>/text/refined_*.csv`: refined event groups and spatial-exclusion flags.
- `<output_dir>/text/classification_results.csv`: sequence features and classes.
- `<output_dir>/migration_test/migration_test_summary.csv`: model comparison, migration decision and velocity estimates.
- `<output_dir>/migration_test/migration_test_status.csv` and `migration_test_runs.csv.gz`: processing status and individual test runs.

## Example data and reproducibility

The regional and small examples are extracted from a **locally prepared QTM-derived** catalog. Earlier preprocessing of the project source has not been independently verified; selection boundaries, magnitude threshold and data references are recorded in [examples/README.md](examples/SAN_JACINTO_REGION.md). Neither example reproduces the full California study: the paper's 246-sequence sample includes additional review and case-specific handling not automatically reproduced by this general workflow.

The QTM catalog is described by [Ross et al. (2019)](https://doi.org/10.1126/science.aaw6888); see the [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html) for catalog access and data acknowledgments.

## Tests and coverage

The [Tests workflow](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml) runs unit tests and the small end-to-end example on Python 3.10–3.12. Coverage is reported in the workflow and provided as an XML artifact for each Python version.

## Citation

If you use ASPECT, please cite the software and the earthquake catalogs used in your analysis:

> Cui, X. (2026). *ASPECT: Earthquake sequence identification and migration analysis* [Computer software]. GitHub. https://github.com/xcui1997/ASPECT

```bibtex
@misc{cui2026aspect,
  author = {Cui, Xin},
  title  = {ASPECT: Earthquake sequence identification and migration analysis},
  year   = {2026},
  url    = {https://github.com/xcui1997/ASPECT}
}
```

For reproducibility, record the commit hash used in your analysis. Citation metadata is also available in [CITATION.cff](CITATION.cff).

## Development acknowledgment

ASPECT was developed by Xin Cui in collaboration with [OpenAI Codex](https://openai.com/codex/), which assisted with implementation, testing, and documentation. Individual Codex-assisted commits may also identify Codex as a co-author in their commit metadata.

## License and attribution

ASPECT is released under the MIT License; see [LICENSE](LICENSE). The magnitude-completeness and b-value helper `McB.py` was adapted from [sachalapins/bvalues](https://github.com/sachalapins/bvalues); its original license is preserved in the project attribution files.
