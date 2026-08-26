# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project state

This repository is a freshly scaffolded `uv` Python project (created via `uv init`). It currently contains only
the package skeleton — no application logic, tests, or dependencies have been added yet. Expect to build most
things from scratch rather than finding existing conventions to follow.

## Purpose

Based on the project and data-directory names, this repo is for evaluating WRF (Weather Research and Forecasting
model) simulations at the BNF site (DOE ARM Bankhead National Forest site) against reference/forcing datasets.

## Environment and commands

Package management is via `uv` (build backend `uv_build`, package layout under `src/`).

- Install/sync dependencies: `uv sync`
- Add a dependency: `uv add <package>`
- Run the CLI entry point: `uv run bnf-wrf-eval` (maps to `bnf_wrf_eval:main`)
- Run any script/module in the project's environment: `uv run python -m <module>`
- Build the package: `uv build`

Python version is pinned to 3.11 (`.python-version`, `requires-python = ">=3.11"` in `pyproject.toml`).

No test runner, linter, or formatter is configured yet — check `pyproject.toml` for `[dependency-groups]`/`[tool.*]`
sections before assuming a tool (e.g. pytest, ruff) is in use.

## Structure

- `src/bnf_wrf_eval/` — the only Python package in the repo (`src` layout). `__init__.py` currently just defines
  `main()`, the entry point wired up in `pyproject.toml`.
- `satoshi_forcing_data/` and `satoshi_testruns/` — symlinks into shared project storage
  (`/gpfs/wolf2/arm/cli120/proj-shared/sey/bnf/wrf/...`), **not part of the git repo** (untracked, and large).
  - `satoshi_forcing_data/` holds reference/forcing datasets (`era5/`, `era5rda/`, `hrrr/`).
  - `satoshi_testruns/` holds WRF LASSO BNF test-run output directories, each named for a config/date
    (e.g. `20250502lassobnfwrfera5ml3`, `20250917lassobnfwrfhrrr2`). These are read-only reference data owned by
    another user (`sey`) — do not modify files under these symlinked paths.
