"""Merge YAML configuration and resolve input and output paths.

Projected coordinates and depth are in km; stime is in seconds and migration time is in days.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config.yaml"


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_config(config_path: str | Path | None = None, region: str | None = None) -> dict[str, Any]:
    default_cfg = load_yaml(DEFAULT_CONFIG_PATH)

    if config_path is None:
        if region is None:
            raise ValueError("Either --region or --config must be provided.")
        config_path = PROJECT_ROOT / "regions" / region / "config.yaml"
    else:
        config_path = Path(config_path)

    config_path = Path(config_path).resolve()
    region_dir = config_path.parent
    overrides = load_yaml(config_path)
    cfg = deep_merge(default_cfg, overrides)
    # A user-specified catalog mapping replaces the example mapping entirely.
    if "data_sources" in overrides:
        cfg["data_sources"] = deepcopy(overrides["data_sources"])
    for flag in ("merge_clusters", "plot_diagnostics"):
        if not isinstance(cfg.get("dbscan_params", {}).get(flag, True), bool):
            raise ValueError(f"dbscan_params.{flag} must be a boolean")
    cfg["_config_path"] = str(config_path)
    cfg["_region_dir"] = str(region_dir)

    paths = cfg.setdefault("paths", {})
    output_dir = Path(paths.get("output_dir", cfg.get("output_dir", "results")))
    if not output_dir.is_absolute():
        output_dir = region_dir / output_dir
    paths["output_dir"] = str(output_dir)
    cfg["output_dir"] = str(output_dir)

    raw_sources = cfg.get("data_sources", {})
    if isinstance(raw_sources, (str, Path)):
        raw_sources = {"single_sequence": raw_sources}
    elif isinstance(raw_sources, list):
        raw_sources = {
            Path(source_path).stem or f"source_{i}": source_path
            for i, source_path in enumerate(raw_sources)
        }

    data_sources = {}
    for source_name, source_path in raw_sources.items():
        p = Path(source_path)
        if not p.is_absolute():
            p = region_dir / p
        data_sources[source_name] = str(p)
    cfg["data_sources"] = data_sources

    return cfg


def ensure_output_dirs(cfg: dict[str, Any]) -> tuple[Path, Path, Path]:
    out = Path(cfg["output_dir"])
    text_dir = out / "text"
    png_dir = out / "png"
    text_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    return out, text_dir, png_dir
