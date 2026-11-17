<p align="center"><img src="assets/aspect_logo.png" alt="ASPECT logo" width="600"></p>

# ASPECT

**English** | [中文](README_zh.md)

ASPECT is a Python workflow for identifying and classifying earthquake sequences and testing their spatiotemporal migration from an earthquake catalog. Its modular algorithms and YAML configuration support both complete runs and step-by-step exploration.

## Features

- Extract preliminary sequences with nearest-neighbor analysis (NNA), then refine them with spatial Area and temporal Knee clustering.
- Classify sequences using interpretable timing, magnitude and tree-structure features.
- Construct a rolling P70 distance front and compare growth models with non-migration references.
- Assess detection stability using the full sequence and 18 random 90% subsets; report individual runs and jointly passing velocity estimates.
- Run the workflow from a shell script or explore the San Jacinto example in a Jupyter notebook.

Spatial clustering uses projected horizontal coordinates; migration distances use three-dimensional coordinates including depth. Automatic UTM projection is intended for a geographically local catalog; broad regions, polar data, or data crossing the antimeridian require an appropriate projection strategy. This general workflow processes the supplied catalog without cross-catalog deduplication.

## Install

Python 3.10 or newer is required. From the project directory:

```bash
python -m pip install -e .
```

For the notebook, also install JupyterLab:

```bash
python -m pip install jupyterlab
```

## Quick start

The default configuration runs a small QTM-derived example near the northwestern San Jacinto fault zone:

```bash
bash run_all.sh
```

Results are saved under `results/`. The command runs NNA, spatial/temporal refinement, classification and migration testing. To run an individual stage after its required inputs have been produced:

```bash
python run_aspect.py --config config.yaml --steps nna
python run_aspect.py --config config.yaml --steps dbscan
python run_aspect.py --config config.yaml --steps classify
python run_aspect.py --config config.yaml --steps migration_test
```

To inspect each stage and its figures, open [notebooks/san_jacinto_demo.ipynb](notebooks/san_jacinto_demo.ipynb):

```bash
jupyter lab notebooks/san_jacinto_demo.ipynb
```

For your own catalog, copy and edit `config.yaml`, then run `bash run_all.sh --config /path/to/config.yaml`. The input CSV must contain `datetime`, `lat`, `lon`, `dep` and `mag`: times are UTC, latitude and longitude are in degrees, and depth is in km. Paths in a configuration file are resolved relative to that file. Review magnitude thresholds, clustering and migration parameters for the target catalog; the example defaults are not universal.

## Model comparison

The example sets the step and logarithmic reference margins to 0 and 4 IC units, respectively. The reported `dic` is `min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_growth`; therefore it is a margin-adjusted comparison, not an unadjusted IC difference against every reference. Both margins are configurable. A run jointly passes when `dic > 0` and the migration-phase R² exceeds 0.5 under the example configuration.

## Output

- `results/text/nna_*.csv`: preliminary event assignments and nearest-neighbor information.
- `results/text/refined_*.csv`: refined event groups and spatial-exclusion flags.
- `results/text/classification_results.csv`: sequence features and classes.
- `results/migration_test/migration_test_summary.csv`: model comparison, migration decision and velocity estimates.
- `results/migration_test/migration_test_status.csv` and `migration_test_runs.csv.gz`: processing status and individual test runs.

## Example data and reproducibility

The example contains 136 earthquakes from a short time and spatial window in a **locally prepared QTM-derived** catalog. See [examples/README.md](examples/README.md) for the exact selection, source provenance and references. Earlier preprocessing of the project source has not been independently verified. The example demonstrates the workflow; it is not a full-catalog or full-paper reproduction. In particular, the paper's reported 246-sequence sample includes additional review and case-specific handling that are not automatically reproduced by this general workflow.

The QTM catalog is described by [Ross et al. (2019)](https://doi.org/10.1126/science.aaw6888); see the [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html) for catalog access and data acknowledgments.

## License and attribution

ASPECT is released under the MIT License; see [LICENSE](LICENSE). The magnitude-completeness and b-value helper `McB.py` was adapted from [sachalapins/bvalues](https://github.com/sachalapins/bvalues); its original license is preserved in the project attribution files.
