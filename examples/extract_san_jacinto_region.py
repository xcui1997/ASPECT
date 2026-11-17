"""Extract a regional QTM-derived catalog without preselecting sequences.

Example: python examples/extract_san_jacinto_region.py /path/to/qtm_relo.csv
The output contains every source event in a 15-km-wide (each side) corridor
around a simplified regional line, at all magnitudes and all available dates.
The line is an approximate geographic selector, not a mapped fault trace.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import Transformer

# Approximate NW-to-SE selector from San Bernardino towards the Salton region.
SELECTOR_LON_LAT = [(-117.50, 34.15), (-117.17, 33.93), (-116.98, 33.77),
                    (-116.70, 33.55), (-116.38, 33.25), (-115.95, 33.00)]
CORRIDOR_KM = 15.0


def corridor_distance_km(lon, lat):
    """Return distance to the nearest segment in projected UTM 11N coordinates."""
    transform = Transformer.from_crs(4326, 32611, always_xy=True)
    vertices = np.asarray([transform.transform(*pair) for pair in SELECTOR_LON_LAT]) / 1000.0
    east, north = transform.transform(np.asarray(lon), np.asarray(lat))
    points = np.column_stack([east, north]) / 1000.0
    distance = np.full(len(points), np.inf)
    for start, end in zip(vertices[:-1], vertices[1:]):
        vector = end - start
        fraction = np.clip((points - start) @ vector / (vector @ vector), 0, 1)
        closest = start + fraction[:, None] * vector
        distance = np.minimum(distance, np.linalg.norm(points - closest, axis=1))
    return distance


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
    mask = df.lat.notna() & df.lon.notna()
    eligible = df.loc[mask].copy()
    selected = eligible.loc[corridor_distance_km(eligible.lon, eligible.lat) <= CORRIDOR_KM].copy()
    selected.sort_values(["datetime", "event_id"], kind="stable", inplace=True)
    selected.to_csv(args.output, index=False)
    print(f"Wrote {len(selected):,} events to {args.output}")


if __name__ == "__main__":
    main()
