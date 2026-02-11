# PRODUCT

## Purpose

Glyph Soup exists to test whether structured complexity signals can emerge from simple reaction rules in artificial chemistry, and to provide reproducible evidence for or against that hypothesis.

## Target users

- Computational researchers exploring emergence and assembly-like metrics
- Engineers building reproducible simulation pipelines
- Contributors extending reaction models and experiments

## Core value

- Reproducible experiment runs via explicit config + seeded RNG
- Clear, inspectable outputs (trace CSV + summary JSON + batch calibration)
- Extensible architecture for additional experiments and reaction policies

## Current features

- Binary-tree molecule representation
- Stochastic bond and break dynamics
- Incremental `a_total` tracking with optional verification
- Single-seed and batch Experiment A execution
- Batch-level analysis, including transition acceleration detection
- Golden trace regression for deterministic behavior checking

## Business / research objectives

- Establish a reliable baseline (null model) for random chemistry
- Calibrate thresholds from simulated distributions rather than imported constants
- Enable rapid iteration on new hypotheses without breaking reproducibility

## Out of scope (for now)

- Production web app / UI productization
- Distributed compute orchestration
- Cross-language runtime ports as default path
