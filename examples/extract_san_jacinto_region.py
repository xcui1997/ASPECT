"""Extract all events in a rectangular San Jacinto region.

Example: python examples/extract_san_jacinto_region.py /path/to/qtm_relo.csv
No magnitude or temporal cut is imposed; all available source events whose
coordinates lie inside the inclusive rectangle are retained.
"""
import argparse
from pathlib import Path
import pandas as pd

LON_MIN, LON_MAX = -117.4, -116.4
LAT_MIN, LAT_MAX = 33.0, 34.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Locally prepared QTM-derived CSV")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("san_jacinto_region_qtm.csv"))
    args = parser.parse_args()
    df = pd.read_csv(args.source)
    required = {"datetime", "event_id", "lat", "lon", "dep", "mag", "rel"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")
    mask = df.lon.between(LON_MIN, LON_MAX, inclusive="both") & df.lat.between(LAT_MIN, LAT_MAX, inclusive="both")
    selected = df.loc[mask].copy()
    selected.sort_values(["datetime", "event_id"], kind="stable", inplace=True)
    selected.to_csv(args.output, index=False)
    print(f"Wrote {len(selected):,} events to {args.output}")


if __name__ == "__main__":
    main()
