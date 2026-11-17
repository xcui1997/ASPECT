# San Jacinto QTM examples

For the full regional example and selection details, see [SAN_JACINTO_REGION.md](SAN_JACINTO_REGION.md).

## Small smoke-test example

`san_jacinto_qtm.csv` is a small example of earthquakes near the northwestern San Jacinto fault zone, extracted from the QTM-derived, locally prepared relocation catalog used by the ASPECT project. It contains **136 events** selected with:

- Origin time: 2013-06-01 00:00:00 UTC (inclusive) to 2013-07-01 00:00:00 UTC (exclusive).
- Latitude: 34.00°–34.15° N (inclusive).
- Longitude: −117.55° to −117.42° (inclusive).
- Magnitude: M ≥ 0.30.

The demonstration intentionally uses a short time span and geographic window; its outcomes are not a full-catalog reproduction and may depend on these boundaries. The earliest selected event is 2013-06-14 03:09:11.461 UTC and the latest is 2013-06-29 14:42:27.675 UTC.

Columns are `datetime` (UTC origin time), `event_id` (source identifier), `lat` and `lon` (degrees), `dep` (depth in km), `mag` (catalog magnitude), and `rel` (flag retained unchanged from the project input). No interpretation of `rel` is required by the example workflow.

The source file used for this subset is the existing project input `qtm_relo.csv` (639,723 records; SHA-256 `3c861f5beb09800b7e12b9679189558669393848c5f2a3b7146c6d5eb3fe5b22`). It is a locally prepared QTM-derived input, not a fresh extraction directly from an SCEDC distribution. Its earlier preprocessing and relocation provenance have not been independently established here; the example reproduces the records as present in this project input. The subset SHA-256 is `d5104f38d1a9df81169abb7a0a12680501a1766f2da027fc4113495e4138fab3`.

The extraction can be repeated from that prepared source with `python examples/extract_qtm_example.py /path/to/qtm_relo.csv /path/to/output.csv`. The extraction script reproduces the checked-in sample byte-for-byte from the source identified above.

## Data citation

The QTM catalog is described by Ross, Z. E., Trugman, D. T., Hauksson, E., & Shearer, P. M. (2019), Searching for hidden earthquakes in Southern California, *Science*, 364, 767–771, https://doi.org/10.1126/science.aaw6888. The Southern California Earthquake Data Center provides the catalog at https://scedc.caltech.edu/data/qtm-catalog.html. Its data acknowledgement identifies the Southern California Seismic Network (https://doi.org/10.7914/SN/CI) and SCEDC (https://doi.org/10.7909/C3WD3xH1).
