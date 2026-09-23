# Response to peer review (working checklist)

Version 0.5 is framed as a technical scoping note. It does not claim a
long-baseline design or a new source-derived feasibility result. The reviewer
recommendations are divided into implemented editorial/analysis changes and
the source and detector work still needed for a stronger feasibility paper.

| Review issue | Implemented in version 0.5 | Still required |
| --- | --- | --- |
| Title implies empirical CRC packets | Retitled to “NuMI-Calibrated Simulation and Far-Field Sensitivity”; abstract, methods, tables and conclusion call the CRC results a simulation | None for terminology |
| Far-field source divergence is free | Far-field figures and captions consistently call out conditional geometric sensitivity; no source requirement is claimed | Published/simulated energy-angle flux and useful-neutrino power for a chosen facility |
| Detector represented by mass only | Limitations remain prominent; T2K iron cross section is only a comparison point, not a generic response | Actual projected geometry, material, efficiency, containment, backgrounds and flux fold |
| Propagation omitted | Labeled vacuum two-flavor estimate retained only as a sensitivity illustration; reference [9] now has complete bibliographic details | Three-flavor matter propagation through a specified Earth profile and spectral source |
| CRC packet results do not use measured event records | Wording corrected throughout; packet simulation uses the published count mean and pulse schedule but is not described as measured or as a decoder reproduction | Reconstruct observed counts/frame sync/code or obtain event-level data |
| Input-mean uncertainty omitted | Re-estimate `lambda = 2×1402/3454`; propagate approximate Poisson 95% limits through repeated-message simulation; distinguish this from Monte Carlo uncertainty | Selection/systematic uncertainty and observed-data comparison |
| Power/rate conventions unclear | Define slot duration/rate, raw OOK bit rate, peak on-slot power, full-duty average and input-weighted mean in Section 2; standardize equations, captions and plots | Apply a physical facility model to convert beam/facility input energy |
| Pulse schedule averaged without phase | Simulate 25 pulses at 2.2 s spacing plus the stated supercycle gap; report latency range over uniform arrival phase | Validate schedule interpretation against original timing records |
| Coding/throughput is assumed | Rename as nominal bookkeeping; add Stancil's 0.37 bit/slot zero-background capacity as a theoretical reference; no FEC gain is claimed | Compare concrete finite-block codes at matched average energy and slot schedule |
| Compression ratio is hypothetical | Keep 4:1 example explicitly secondary and unmeasured | Test compressor on an identified corpus |
| Solar System prior art/units | Cite planetary-blockage discussion in [1]; Appendix reports peak one-slot power and explains long integration intervals | Body transmission, geometry and ephemeris simulation |
| Reproducibility | README includes commands for scripts, tables, figures and tests; the draft PR head is an immutable Git commit that identifies this revision | Archive/release if an archival citation is needed |

## Prioritized next analysis

1. Choose a source with published energy- and angle-resolved flux and a proton
   or muon input budget. Compute delivered flux at the receiver instead of
   assigning divergence independently.
2. Fold that spectrum through three-flavor oscillations and a stated Earth
   density model, target material and detector response; include backgrounds
   and uncertainties.
3. Simulate complete messages on the real source timing at matched average
   facility energy. Compare OOK/repetition with a finite-block code or sparse
   pulse scheme, including acquisition and packet failure.
4. Validate against observed detector counts where accessible. Keep submarine,
   finance and occultation use cases out of feasibility claims until each has
   an appropriate source and receiver model.

## Packet simulation scope

`empirical_message.py` samples uniform five-byte payloads with CRC-8/ATM and
one-sided missed detections derived from a Poisson mean. Packet timing uses
the published 25-pulse supercycle and uniformly distributed issue phase.
Receivers receive slot synchronization externally. It records exact correct
acceptance, CRC acceptance, undetected errors, slot latency and proton-beam
energy incident on the target. The input-mean interval reflects only Poisson
counting statistics in the event count used by Stancil et al.; it does not
cover selection systematics. This simulation is not the experiment's decoder.
