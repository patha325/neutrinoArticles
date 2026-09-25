# Paper II working package

**Title:** *End-to-End Simulation of a Neutrino Communication Channel Through the Earth*

This folder contains the new source-to-decoder follow-up to Paper I. The manuscript is currently a methods draft, not a submission-ready paper. It deliberately contains no long-baseline numerical claim because the source configuration and detector response have not yet been fixed and validated.

## Scope

The analysis will replace Paper I's conditional event-rate input with:

1. an explicit energy- and angle-resolved accelerator-neutrino source;
2. geometric transport and three-flavor matter propagation through a stated Earth density profile;
3. interaction cross sections folded with detector geometry, material, and selection efficiency;
4. time-dependent signal and background counts in actual communication slots;
5. synchronization, framing, finite-block decoding, and complete-message outcomes;
6. facility-input energy per successfully delivered payload bit.

The measured NuMI–MINERvA communication demonstration is a validation reference only where the published inputs permit. It will not be described as a reproduction of its full decoder without the required event, timing, and protocol inputs.

## Manuscript

- [paper_ii_draft.md](manuscript/paper_ii_draft.md) — current manuscript draft and model definitions.
- [proposal.md](proposal.md) — research proposal connecting Paper II to the few-event capacity and coding study.

## Result gate

No result may be described as source-derived until the repository contains a versioned source flux input and provenance, a detector response definition, a propagation implementation with numerical checks, background assumptions, and a decoder run that emits complete-message metrics. All rate, energy, and power terms must retain the distinctions defined in the manuscript.
