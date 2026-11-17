# San Jacinto regional QTM example

`san_jacinto_region_qtm.csv` contains **288,129 events at all magnitudes**, selected from every available date in the project's locally prepared QTM-derived catalog, from 2008-01-01 through 2017-12-31 (UTC). Every input event is at most 15 km from the geographic selector described below. No sequences, magnitude classes, or time windows were manually selected. The earliest selected event is 2008-01-01 00:04:12.531 UTC and the latest is 2017-12-31 23:44:19.261 UTC. The 2008–2017 temporal coverage is limited to that of the source file, not the full historical or current QTM archive.

This is a regional demonstration of the San Jacinto area, **not an authoritative fault-zone polygon**. The selector is a deliberately simplified line connecting (longitude, latitude, degrees):

```
(-117.50, 34.15), (-117.17, 33.93), (-116.98, 33.77),
(-116.70, 33.55), (-116.38, 33.25), (-115.95, 33.00)
```

Distances to its segments are measured in UTM zone 11N (EPSG:32611). All events within 15 km on either side or around the line's endpoints are included, irrespective of later cluster assignment. The selector is approximate; it can include nearby off-fault and background events, and it can exclude fault-zone events beyond the corridor. A mapped fault-trace dataset and its spatial uncertainty would be needed for a more formal fault-association study. The [Southern California Earthquake Data Center](https://scedc.caltech.edu/earthquake/sanjacinto.html) describes the broader multi-strand San Jacinto fault system.

The CSV contains all 288,129 source events in the corridor. The regional configuration applies **M >= 0.3** at the NNA stage, leaving 80,574 events for analysis and excluding 207,555 smaller events from the subsequent workflow. The magnitude threshold is a workflow setting, not a demonstration that the catalog is complete throughout the whole area and interval. Assess magnitude completeness for a scientific analysis.

## Reproduce the selection

Given the same project-preprocessed input file, run:

```bash
python examples/extract_san_jacinto_region.py /path/to/qtm_relo.csv
```

The source file has 639,723 records and SHA-256 `3c861f5beb09800b7e12b9679189558669393848c5f2a3b7146c6d5eb3fe5b22`. The resulting selection has SHA-256 `c0881882f765484244d7b94374793d13630c7174eee885ef36b694e6a571f1ac`. It is a locally prepared QTM-derived catalog, not a fresh extraction directly from an SCEDC distribution. Its earlier preprocessing and relocation history have not been independently established. The retained columns are `datetime` (UTC), `event_id`, `lat`, `lon` (degrees), `dep` (km), `mag`, and the source's `rel` flag. The extraction does not reinterpret or filter `rel`.

To run the regional analysis, install the package with optional acceleration and supply the regional configuration:

```bash
python -m pip install -e '.[acceleration]'
bash run_all.sh --config examples/config_san_jacinto_region.yaml
```

Results go to `results/san_jacinto_region/`. NNA uses an exact spatial pruning algorithm; numba is strongly recommended at this scale. Runtime depends on hardware and the Python environment. Clustering and classification apply the thresholds in the configuration; they are not equivalent to mapping all physical earthquake sequences, and their numerical results depend on the selected spatial and magnitude bounds.

## Data citation

The QTM catalog is described by Ross, Z. E., Trugman, D. T., Hauksson, E., & Shearer, P. M. (2019), Searching for hidden earthquakes in Southern California, *Science*, 364, 767–771, https://doi.org/10.1126/science.aaw6888. The [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html) provides catalog information. Its data acknowledgement identifies the Southern California Seismic Network (https://doi.org/10.7914/SN/CI) and SCEDC (https://doi.org/10.7909/C3WD3xH1).
