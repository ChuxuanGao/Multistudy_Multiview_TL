# AGENTS.md

Guidelines for AI agents (Codex) working on this repository.

## Project Overview

This repository implements extensions of the **ptLasso framework** for transfer learning and multi-study modeling.

The current development goal is to extend ptLasso to support **multistudy multiview learning** 
using the the extended version of `cv.multiview` framework called `cvar.multiview` which enables alpha and rho tuning.

Primary languages: **R**

Core modeling paradigms:

- Penalized regression
- Transfer learning
- Multi-study learning
- Multi-view data integration

The implementation should remain consistent with the **original ptLasso design philosophy**.

---

# Key Design Principles

1. **Preserve ptLasso semantics**

Existing ptLasso workflows should not change unless absolutely necessary.

2. **Minimal API changes**

Function interfaces should remain consistent with existing ptLasso and glmnet conventions.

3. **Follow R modeling conventions**

Functions should follow common R package patterns:

- `model()`
- `cv.model()`
- `predict.model()`
- `predict.cv.model()`

4. **Readable statistical code**

Code should prioritize:

- clarity
- reproducibility
- explicit handling of model arguments.

---

# Core Functions

The following functions define the multiview ptLasso pipeline:
