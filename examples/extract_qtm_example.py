"""Extract the small San Jacinto demonstration from the project's QTM CSV.

Usage: python examples/extract_qtm_example.py /path/to/qtm_relo.csv output.csv
The source is a locally prepared catalog; see README.md for provenance.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Path to the project QTM CSV.")
    parser.add_argument("output", type=Path, help="Path for the extracted CSV.")
    args = parser.parse_args()
    events = pd.read_csv(args.source)
    times = pd.to_datetime(events["datetime"], utc=True, errors="raise")
    mask = (
        (times >= pd.Timestamp("2013-06-01", tz="UTC"))
        & (times < pd.Timestamp("2013-07-01", tz="UTC"))
        & events["lat"].between(34.00, 34.15, inclusive="both")
        & events["lon"].between(-117.55, -117.42, inclusive="both")
        & (events["mag"] >= 0.3)
    )
    events.loc[mask].sort_values("datetime").to_csv(args.output, index=False)
    print(f"Wrote {mask.sum()} events to {args.output}")


if __name__ == "__main__":
    main()
