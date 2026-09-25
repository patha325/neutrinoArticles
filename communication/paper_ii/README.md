# Paper II working package

**Title:** *End-to-End Simulation of a Neutrino Communication Channel Through the Earth*

This package now contains a runnable DUNE benchmark, seven generated figures, CSV/JSON outputs, and a revised article draft. The benchmark uses public DUNE TDR inputs. Its rate is source-folded for a 40 kt far detector, but detector migration, spill timing, communication backgrounds, and synchronization acquisition are not modeled. The rounded headline is a conditional benchmark, not a detector prediction: no physical uncertainties are propagated, the reconstructed-energy efficiency is used against true energy, and the flux-plane/baseline mapping is unconfirmed. Results should be read with the assumptions in the article and `data/simulation/simulation_summary.json`.

## Contents

- [Article draft](manuscript/paper_ii_article_draft.md) — v0.1 with methods, benchmark results, limits, and next steps.
- [Research proposal](proposal.md) — linked Paper II and future few-event capacity/coding study.
- [Methods draft](manuscript/paper_ii_draft.md) — the prior model plan and definitions.
- `analysis/neutrino_channel.py` — PREM and three-flavor evolution, binned rate fold, Poisson channels, and packet simulation.
- `analysis/run_simulation.py` — reproducible CSV/JSON and figure generation.
- `analysis/test_neutrino_channel.py` — numerical invariants and limiting-case checks.
- `data/inputs/` — public DUNE text inputs and [provenance](data/inputs/PROVENANCE.md).
- `data/simulation/` — folded spectrum, slot scale, equal-event-budget OOK/PPM comparison, quadrature/profile-step convergence, independent matrix-exponential cross-check, finite packet outcomes, and assumption summary.
- `figures/` — SVG figures for the manuscript and PNG previews.

## Reproduce

```bash
python -m pip install -r communication/paper_ii/analysis/requirements.txt
python communication/paper_ii/analysis/run_simulation.py
python -m unittest discover -s communication/paper_ii/analysis -p 'test_*.py'
```

The source is a DUNE TDR FHC far-detector flux histogram. PREM propagation is integrated within each flux-energy bin using Gaussian quadrature and checked against higher orders and a separate SciPy matrix-exponential evolution. DUNE's published selection efficiency is applied directly against true energy as a proxy because this pass does not fold the full migration matrix. Calendar-average exposure is used for slot-duration illustrations, with no explicit spill schedule or second uptime factor. Background scenarios are fixed per-slot assumptions, not predictions. The coding comparison reports a conditional incident-proton beam-energy equivalent; accelerator electrical energy and achieved communications rate are not estimated.

Paper II uses a simple physical channel and decoder baseline. *The Neutrino Channel: Capacity and Coding in the Few-Event Regime* remains a distinct follow-on study that can reuse this channel model for its more detailed coding analysis.
