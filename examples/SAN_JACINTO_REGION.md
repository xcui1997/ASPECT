# San Jacinto regional QTM example

`san_jacinto_region_qtm.csv` contains **292,468 events at all magnitudes**, selected from every available date in the project's locally prepared QTM-derived catalog. The inclusive rectangular selection is:

- Longitude: **−117.4° to −116.4°**.
- Latitude: **33.0° to 34.0°**.
- Time: all dates in the source, from **2008-01-01 00:04:18.539 UTC** to **2017-12-31 23:44:19.261 UTC** for the selected events.

No specific sequence, magnitude class, or time window is selected before analysis. The 2008–2017 interval is a limitation of the locally available source file; it is not the full historical or current QTM archive. The rectangle contains parts of the San Jacinto region, other nearby structures, and background earthquakes; it is a geographic window, not a fault-association classification, and does not cover the entire fault system.

The CSV contains all 292,468 source events inside the rectangle. The regional NNA configuration applies **M ≥ 0.3**, leaving **72,061 events** for analysis; **220,407** smaller events do not enter the subsequent workflow. This configurable magnitude threshold is not evidence of catalog completeness throughout the region and interval. Assess completeness for any scientific analysis.

## Reproduce the selection

Given the same project-preprocessed input file, run:

```bash
python examples/extract_san_jacinto_region.py /path/to/qtm_relo.csv
```

The source file contains 639,723 records and has SHA-256 `3c861f5beb09800b7e12b9679189558669393848c5f2a3b7146c6d5eb3fe5b22`. The resulting CSV has SHA-256 `b7bf50fc7f68dcf4565e5f22cc58d9c441ecbeff5284d57a20e26e59d9f73d84`; independent extractions produce the same bytes. It is a locally prepared QTM-derived catalog, not a fresh extraction directly from an SCEDC distribution. Its earlier preprocessing and relocation history have not been independently established. The retained columns are `datetime` (UTC), `event_id`, `lat`, `lon` (degrees), `dep` (km), `mag`, and the source's `rel` flag. The extraction does not reinterpret or filter `rel`.

To run the regional analysis, install the package with optional acceleration and supply the regional configuration:

```bash
python -m pip install -e '.[acceleration,notebook]'
bash run_all.sh --config examples/config_san_jacinto_region.yaml
```

Results go to `results/san_jacinto_region/`. NNA uses an exact spatial pruning algorithm; numba is strongly recommended at this scale. Runtime depends on hardware and the Python environment. Clustering and classification apply the configured thresholds; they are not equivalent to mapping all physical earthquake sequences, and their results depend on the selected spatial and magnitude bounds.

## Data citation

The QTM catalog is described by Ross, Z. E., Trugman, D. T., Hauksson, E., & Shearer, P. M. (2019), Searching for hidden earthquakes in Southern California, *Science*, 364, 767–771, https://doi.org/10.1126/science.aaw6888. The [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html) provides catalog information. Its data acknowledgement identifies the Southern California Seismic Network (https://doi.org/10.7914/SN/CI) and SCEDC (https://doi.org/10.7909/C3WD3xH1).
