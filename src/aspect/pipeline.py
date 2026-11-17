"""Workflow scheduling, input normalization, classification gates, and CSV export.

Projected coordinates and depth are in km; stime is in seconds and migration time is in days.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import zlib
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .classifier import ClusterClassifier
from .config import ensure_output_dirs, load_config
from .migration_test import MigrationTester
from .migration_front import trim_early_events
from .nna_cluster import BPTree
from .st_dbscan import SeismicClusterAnalyzer
from .utils import convert_gps_to_utm, get_refined_b_value, select_region


def _parse_single_datetime(value):
    if pd.isna(value):
        return pd.NaT
    if isinstance(value, pd.Timestamp):
        return value.tz_localize("UTC") if value.tzinfo is None else value.tz_convert("UTC")
    if isinstance(value, datetime):
        ts = pd.Timestamp(value)
        return ts.tz_localize("UTC") if ts.tzinfo is None else ts.tz_convert("UTC")

    text = str(value).strip()
    if not text:
        return pd.NaT

    normalized = text.replace("/", "-")
    normalized = re.sub(r"\s+", " ", normalized)
    normalized = re.sub(r"(\d)\.\d+Z$", r"\g<0>", normalized)
    try:
        parsed = pd.to_datetime(normalized, format="mixed", utc=True, errors="coerce")
    except TypeError:
        parsed = pd.to_datetime(normalized, utc=True, errors="coerce")
    if not pd.isna(parsed):
        return parsed

    for fmt in (
        "%Y-%jT%H:%M:%S.%f",
        "%Y-%jT%H:%M:%S",
        "%Y-%j %H:%M:%S.%f",
        "%Y-%j %H:%M:%S",
        "%Y%j%H%M%S.%f",
        "%Y%j%H%M%S",
        "%Y%j",
        "%Y%m%dT%H%M%S.%f",
        "%Y%m%dT%H%M%S",
        "%Y%m%d %H%M%S.%f",
        "%Y%m%d %H%M%S",
        "%Y%m%d%H%M%S.%f",
        "%Y%m%d%H%M%S",
        "%Y%m%d",
    ):
        try:
            return pd.Timestamp(datetime.strptime(normalized, fmt), tz="UTC")
        except ValueError:
            continue
    raise ValueError(f"Unsupported datetime format: {value!r}")


def parse_datetime_series(values):
    if isinstance(values, pd.Series):
        return values.apply(_parse_single_datetime)
    return _parse_single_datetime(values)


def preprocess_catalog(df: pd.DataFrame, stime_origin: str = "1900-01-01") -> pd.DataFrame:
    df = df.copy()
    df["datetime"] = parse_datetime_series(df["datetime"])
    origin = parse_datetime_series(stime_origin)
    df["stime"] = (df["datetime"] - origin).dt.total_seconds()
    return df.sort_values(by="stime").reset_index(drop=True)


def _safe_file_stem(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_") or "sequence"


def _prepare_sequence_dataframe(
    df: pd.DataFrame,
    cfg: dict,
    source_name: str,
    refined_cluster: int | None = None,
) -> pd.DataFrame:
    out = df.copy()
    if "datetime" in out.columns:
        out["datetime"] = parse_datetime_series(out["datetime"])
        out = out.dropna(subset=["datetime"]).copy()
        if "stime" not in out.columns:
            out = preprocess_catalog(out)
        else:
            out = out.sort_values("stime").reset_index(drop=True)
    elif "stime" in out.columns:
        out = out.sort_values("stime").reset_index(drop=True)
    elif "time" in out.columns:
        out["stime"] = out["time"]
        out = out.sort_values("stime").reset_index(drop=True)
    else:
        raise KeyError("Sequence data must contain one of: datetime, stime, or time.")

    if "time_day" not in out.columns:
        out["time_day"] = (out["stime"] - out["stime"].min()) / 86400.0

    if ("x" not in out.columns or "y" not in out.columns) and {"lat", "lon"}.issubset(out.columns):
        dbscan_cfg = cfg.get("dbscan_params", {})
        out["x"], out["y"], _ = convert_gps_to_utm(
            out["lat"].values,
            out["lon"].values,
            zone=dbscan_cfg.get("zone_code"),
        )

    if "cluster_id" not in out.columns:
        out["cluster_id"] = 0
    if refined_cluster is not None and "refined_cluster" not in out.columns:
        out["refined_cluster"] = refined_cluster
    if "source_name" not in out.columns:
        out["source_name"] = source_name
    return out




def _write_single_refined_file(
    df: pd.DataFrame,
    cfg: dict,
    source_name: str,
    text_dir: Path,
) -> tuple[Path, pd.DataFrame]:
    out = _prepare_sequence_dataframe(df, cfg, source_name, refined_cluster=0)
    out["refined_cluster"] = out.get("refined_cluster", 0)
    out_name = f"refined_{_safe_file_stem(source_name)}_single.csv"
    out_path = text_dir / out_name
    out.to_csv(out_path, index=False)
    return out_path, out

def run_nna(cfg: dict) -> None:
    _, text_dir, png_dir = ensure_output_dirs(cfg)
    nna_cfg = cfg.get("nna_params", {})

    for source_name, path in cfg.get("data_sources", {}).items():
        print(f"--- Processing NNA for {source_name} ---")
        df_raw = pd.read_csv(path)
        required = {"datetime", "lat", "lon", "dep", "mag"}
        missing = required - set(df_raw.columns)
        if missing:
            raise ValueError(f"Catalog {source_name} is missing required columns: {sorted(missing)}")
        if df_raw["datetime"].isna().any():
            raise ValueError(f"Catalog {source_name} contains missing event times")
        for column in ("lat", "lon", "dep", "mag"):
            values = pd.to_numeric(df_raw[column], errors="coerce")
            if not np.isfinite(values).all():
                raise ValueError(f"Catalog {source_name} contains invalid {column} values")
            df_raw[column] = values
        if not df_raw["lat"].between(-90, 90).all() or not df_raw["lon"].between(-180, 180).all():
            raise ValueError(f"Catalog {source_name} contains invalid latitude or longitude")
        mag_min = nna_cfg.get("mag_thresholds", {}).get(source_name, nna_cfg.get("minimum_magnitude"))
        df = df_raw[df_raw["mag"] >= mag_min].copy() if mag_min is not None else df_raw.copy()
        if len(df) < 2:
            raise ValueError(f"Catalog {source_name} must contain at least two selected events")
        df = preprocess_catalog(df)

        tree = BPTree(
            catalog=df,
            name=source_name,
            b=nna_cfg.get("b_value", 1.0),
            df=nna_cfg.get("df_value", 1.6),
            p=nna_cfg.get("p_value", 0.5),
            lag=nna_cfg.get("lag"),
            bp_mode=nna_cfg.get("bp_mode", "exact"),
            bp_screening_margin=nna_cfg.get("bp_screening_margin", 1.2),
            bp_bootstrap_candidates=nna_cfg.get("bp_bootstrap_candidates", 32),
            bp_verify_fast_result=nna_cfg.get("bp_verify_fast_result", False),
            bp_progress_every=nna_cfg.get("bp_progress_every"),
        )
        tree.grow()
        tree.prune(c=nna_cfg.get("prune_c"))
        # Replace the destination entry atomically: reused NNA inputs can be
        # symlinks to a frozen run, which must never be overwritten in place.
        nna_path = text_dir / f"nna_{source_name}.csv"
        temporary = nna_path.with_suffix(".csv.tmp")
        tree.save_results(temporary)
        temporary.replace(nna_path)
        tree.plot_heatmap(png_dir / f"{source_name}_heatmap.pdf")


def run_dbscan(cfg: dict) -> None:
    _, text_dir, _ = ensure_output_dirs(cfg)
    analyzer = SeismicClusterAnalyzer(cfg)
    do_physics = cfg.get("calc_b_value", False)
    dbscan_cfg = cfg.get("dbscan_params", {})
    analysis_cfg = cfg.get("analysis_settings", {})
    fixed_mc = analysis_cfg.get("fixed_mc")
    refined_paths = []

    for source_name in cfg.get("data_sources", {}):
        file_name = text_dir / f"nna_{source_name}.csv"
        if not file_name.exists():
            raise FileNotFoundError(f"NNA output missing: {file_name}. Run the nna stage first.")
        print(f"--- Processing ST-DBSCAN for {source_name} ---")
        full_df = _prepare_sequence_dataframe(pd.read_csv(file_name), cfg, source_name)

        for cid, cluster in full_df.groupby("cluster_id"):
            if len(cluster) < cfg.get("min_prelim_size", 5):
                continue
            if fixed_mc is not None:
                selected = cluster[cluster["mag"] >= float(fixed_mc)].copy()
                suffix = f"_Mc{float(fixed_mc):.1f}"
            elif do_physics:
                bg = select_region(
                    full_df,
                    cluster,
                    radius_km=analysis_cfg.get("radius_km", 100.0),
                    period_years=analysis_cfg.get("period_years", 5),
                )
                b_val, mc = get_refined_b_value(bg["mag"].values, mbin=0.1)
                selected = cluster[cluster["mag"] >= mc].copy()
                suffix = f"_b{b_val:.1f}_Mc{mc:.1f}"
            else:
                selected = cluster.copy()
                suffix = "_Direct"

            if len(selected) < cfg.get("min_prelim_size", 5):
                continue

            main_ev = selected.loc[selected["mag"].idxmax()]
            if "datetime" in selected.columns:
                main_ev_time = pd.to_datetime(main_ev["datetime"], utc=True, errors="coerce")
                date_tag = main_ev_time.strftime("%Y%m%d") if not pd.isna(main_ev_time) else "undated"
            else:
                date_tag = "undated"
            out_name = f"{source_name}_C{cid}_{date_tag}_M{main_ev['mag']:.1f}{suffix}"
            clustered = analyzer.run_sequential_dbscan(selected, file_tag=out_name)
            final = analyzer.refine_clusters(clustered)
            output_path = text_dir / f"refined_{out_name}.csv"
            final.to_csv(output_path, index=False)
            refined_paths.append(output_path.name)
            if dbscan_cfg.get("plot_clusters", True):
                analyzer.plot_clusters(final, file_name=out_name)

    # Record the current stage outputs so old files cannot enter classification.
    manifest_path = text_dir / "refined_manifest.json"
    manifest_path.write_text(json.dumps(sorted(refined_paths), indent=2) + "\n")


def _classification_minimum_size(cfg: dict, parent_size: int) -> float:
    """Minimum retained members: fixed threshold and parent fraction both apply."""
    minimum = cfg.get("analysis_settings", {}).get("min_events_per_cluster", 50)
    fraction = cfg.get("classification_params", {}).get("min_parent_fraction", 0.1)
    return max(minimum, fraction * parent_size)


def _sha256(path: Path) -> str:
    """Hash a file without loading the entire file into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _refined_snapshot(text_dir: Path) -> dict:
    """Identify the manifest and exact refined files used for classification."""
    manifest = text_dir / "refined_manifest.json"
    if not manifest.exists():
        raise FileNotFoundError("Refinement manifest missing. Run the dbscan stage first.")
    names = json.loads(manifest.read_text())
    if not isinstance(names, list) or not all(isinstance(name, str) for name in names):
        raise ValueError("Invalid refinement manifest. Run the dbscan stage again.")
    files = {}
    for name in names:
        if Path(name).name != name:
            raise ValueError(f"Invalid refined filename in manifest: {name}")
        path = text_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"Refined file missing: {path}. Run dbscan and classify again.")
        files[name] = _sha256(path)
    return {"manifest_sha256": _sha256(manifest), "refined_sha256": files}


def _classification_settings(cfg: dict) -> dict:
    """Retain settings that influence classification and its eligibility gate."""
    nna = cfg.get("nna_params", {})
    return {"classification_params": cfg.get("classification_params", {}),
            "analysis_settings": cfg.get("analysis_settings", {}),
            "tree_defaults": {"b_value": nna.get("b_value"), "df_value": nna.get("df_value")}}


def _write_classification_provenance(cfg: dict, text_dir: Path) -> None:
    snapshot = _refined_snapshot(text_dir)
    snapshot.update({"classification_sha256": _sha256(text_dir / "classification_results.csv"),
                     "classification_settings": _classification_settings(cfg)})
    (text_dir / "classification_provenance.json").write_text(
        json.dumps(snapshot, sort_keys=True, indent=2) + "\n")


def _validate_classification_provenance(cfg: dict, text_dir: Path) -> None:
    provenance_path = text_dir / "classification_provenance.json"
    if not provenance_path.is_file():
        raise FileNotFoundError("Classification provenance missing. Run the classify stage first.")
    try:
        saved = json.loads(provenance_path.read_text())
    except (ValueError, OSError) as error:
        raise ValueError("Invalid classification provenance. Run classify again.") from error
    current = _refined_snapshot(text_dir)
    current.update({"classification_sha256": _sha256(text_dir / "classification_results.csv"),
                    "classification_settings": _classification_settings(cfg)})
    if saved != current:
        raise ValueError("Classification is stale: refined files or classification settings changed. "
                         "Run classify again before migration_test.")

def run_classify(cfg: dict) -> None:
    _, text_dir, png_dir = ensure_output_dirs(cfg)
    classifier = ClusterClassifier(cfg)
    manifest_path = text_dir / "refined_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError("Refinement manifest missing. Run the dbscan stage first.")
    snapshot = _refined_snapshot(text_dir)
    refined_files = [str(text_dir / name) for name in snapshot["refined_sha256"]]
    if not refined_files:
        raise ValueError("DBSCAN found no sequences for classification.")

    all_metrics = []
    for fpath in refined_files:
        df = _prepare_sequence_dataframe(
            pd.read_csv(fpath),
            cfg,
            Path(fpath).stem.replace("refined_", "", 1),
            refined_cluster=0,
        )
        if "refined_cluster" not in df.columns:
            df["refined_cluster"] = 0
        cluster_data = df[df["refined_cluster"] != -1]
        for cid, cluster in cluster_data.groupby("refined_cluster"):
            # Paper gate: >=50 events and >=10% of the parent sequence;
            # DBSCAN's min_cluster_size is a distinct refinement parameter.
            min_size = _classification_minimum_size(cfg, len(df))
            if len(cluster) < min_size:
                continue
            metrics = classifier.calculate_metrics(cluster)
            if metrics:
                metrics["file_name"] = Path(fpath).name
                metrics["refined_cluster"] = cid
                metrics["lat"] = cluster["lat"].mean() if "lat" in cluster.columns else np.nan
                metrics["lon"] = cluster["lon"].mean() if "lon" in cluster.columns else np.nan
                max_mag_idx = cluster["mag"].idxmax()
                metrics["datetime"] = (
                    cluster.loc[max_mag_idx, "datetime"]
                    if "datetime" in cluster.columns
                    else np.nan
                )
                metrics["mag"] = cluster.loc[max_mag_idx, "mag"]
                metrics["main_lat"] = cluster.loc[max_mag_idx, "lat"] if "lat" in cluster.columns else np.nan
                metrics["main_lon"] = cluster.loc[max_mag_idx, "lon"] if "lon" in cluster.columns else np.nan
                all_metrics.append(metrics)

    if not all_metrics:
        raise ValueError("No refined sequences satisfy the classification requirements.")
    results = classifier.classify(pd.DataFrame(all_metrics))
    results.to_csv(text_dir / "classification_results.csv", index=False)
    _write_classification_provenance(cfg, text_dir)
    classifier.plot_diagnostic(results, png_dir / "classification_diagnostic.png")


def _select_refined_cluster(data: pd.DataFrame, refined_cluster) -> pd.DataFrame:
    cluster = data[data["refined_cluster"] == refined_cluster]
    if not cluster.empty:
        return cluster

    label = str(refined_cluster)
    try:
        value = float(refined_cluster)
        if value.is_integer():
            label = str(int(value))
    except (TypeError, ValueError):
        pass
    return data[data["refined_cluster"].astype(str) == label]




def _migration_test_targets(cfg: dict, text_dir: Path) -> pd.DataFrame:
    """Select classification rows after verifying their upstream provenance."""
    class_res_path = text_dir / "classification_results.csv"
    if class_res_path.exists():
        _validate_classification_provenance(cfg, text_dir)
        targets = pd.read_csv(class_res_path)
        analyze_type = cfg["migration_test_params"].get("analyze_type", "All")
        if analyze_type != "All":
            targets = targets[targets["sequence_type"] == analyze_type]
        return targets
    raise FileNotFoundError("Classification results missing. Run the classify stage first.")


def _load_test_sequences(cfg: dict, text_dir: Path, targets: pd.DataFrame):
    """Yield (row, coords, time_days, mag, seed) for each target sequence (events sorted by time)."""
    params = cfg["migration_test_params"]
    cols = params.get("coordinate_columns", ["x", "y", "dep"])
    base_seed = params.get("perturbation", {}).get("base_seed", 4242)
    for _, row in targets.iterrows():
        f_path = text_dir / row["file_name"]
        if not f_path.exists():
            raise FileNotFoundError(f"Refined file missing: {f_path}. Run dbscan and classify again.")
        data = _prepare_sequence_dataframe(pd.read_csv(f_path), cfg, Path(f_path).stem.replace("refined_", "", 1), refined_cluster=0)
        cluster = _select_refined_cluster(data, row["refined_cluster"]).sort_values("stime")
        if cluster.empty:
            raise ValueError(f"Cluster {row['refined_cluster']} is missing from {f_path}. Run classify again.")
        if len(cluster) < params.get("min_events_required", 50):
            raise ValueError(f"Cluster {row['refined_cluster']} in {f_path} has {len(cluster)} events, "
                             "below the migration min_events_required threshold.")
        # reproducible seed per sequence
        seed = base_seed + zlib.crc32(f"{row['file_name']}|{row['refined_cluster']}".encode()) % 1000003
        yield row, cluster[cols].to_numpy(float), cluster["stime"].to_numpy(float) / 86400.0, cluster["mag"].to_numpy(float), seed


def run_migration_test(cfg: dict) -> None:
    """Revision migration test: binary migration / no migration, reason, parameters and velocity stability from the same 19 runs."""
    out, text_dir, _ = ensure_output_dirs(cfg)
    mig_dir = out / "migration_test"
    (mig_dir / "png").mkdir(parents=True, exist_ok=True)
    params = cfg["migration_test_params"]
    tester = MigrationTester(cfg)
    rows = []
    statuses = []
    run_rows = []
    for row, coords, time, mag, seed in _load_test_sequences(cfg, text_dir, _migration_test_targets(cfg, text_dir)):
        entry = {"file_name": row["file_name"], "refined_cluster": row["refined_cluster"],
                 "sequence_type": row.get("sequence_type", np.nan),
                 "seed": seed}
        res = tester.analyze(coords, time, mag, seed=seed)
        if res is None:
            early = params.get("early_screen", {})
            if early.get("enabled", True):
                _, screened_time, _, _ = trim_early_events(
                    coords, time, mag, early.get("lower_percentile", 10),
                    early.get("upper_percentile", 90), early.get("multiplier", 1.0))
                n_events = len(screened_time)
            else:
                n_events = len(time)
            statuses.append({**entry, "status": "insufficient_after_early_screen",
                             "n_input_events": len(time), "n_events": n_events,
                             "migration": pd.NA})
            continue
        # Save every run to make each joint pass and the speed spread traceable.
        run_rows.extend({"file_name": row["file_name"], "refined_cluster": row["refined_cluster"],
                         "seed": seed, **run} for run in res["runs"])
        entry.update(tester.summary_row(res))
        rows.append(entry)
        statuses.append({**{k: entry[k] for k in ("file_name", "refined_cluster", "sequence_type", "seed")},
                         "status": "analyzed", "n_input_events": len(time),
                         "n_events": entry["n_events"], "migration": entry["migration"]})
        if (res["migration"] and params.get("plot_migration", True)) or (not res["migration"] and params.get("plot_none", False)):
            name = f"{Path(row['file_name']).stem}_cluster{row['refined_cluster']}"
            tester.plot_report(res, title=name, save_path=mig_dir / "png" / f"migration_test_{name}.png")
    pd.DataFrame(run_rows).to_csv(mig_dir / "migration_test_runs.csv.gz", index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(mig_dir / "migration_test_summary.csv", index=False)
    pd.DataFrame(statuses).to_csv(mig_dir / "migration_test_status.csv", index=False)
    if len(summary):
        print(f"Migration test: {int(summary['migration'].sum())} migrations out of {len(summary)} sequences.")
        print(summary["reason"].value_counts().to_string())







def run_steps(cfg: dict, steps: Iterable[str]) -> None:
    """Validate stage names before running them in the configured order."""
    stages = {"nna": run_nna, "dbscan": run_dbscan, "classify": run_classify,
              "migration_test": run_migration_test}
    steps = list(steps)
    unknown = set(steps) - stages.keys()
    if unknown:
        raise ValueError(f"Unknown workflow steps: {sorted(unknown)}")
    for step in steps:
        stages[step](cfg)


def load_and_run(config: str | Path, steps: Iterable[str] | None = None) -> None:
    """Run the configured stages or an explicit subset of stages."""
    cfg = load_config(config_path=config)
    run_steps(cfg, steps if steps is not None else cfg["workflow_steps"])
