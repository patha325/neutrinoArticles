# Response to methodological review (working checklist)

The manuscript remains a scoping paper. The table separates revisions made
from the additional evidence required for a full engineering feasibility claim.

| Review issue | Revision in version 0.4 | Still required |
| --- | --- | --- |
| Undefined long-baseline source | Empirical NuMI benchmark isolated from hypothetical 3 GeV far-field geometry; 120 GeV on-pulse proton energy and published supercycle now tracked separately | Measured or published spectral angular flux, pointing and conversion from facility energy for a chosen long-baseline source |
| Detector/cross section/propagation | Explicit T2K cross-section measurement and a labelled illustrative vacuum oscillation sensitivity; removed implication of physical lower bound | Three-flavor matter propagation with flux-weighted detector-specific cross section, selection and background uncertainties |
| Unsupported payload number | 1.64 nominal bits/s clearly identified as bookkeeping; separate finite-message repeated-OOK/CRC-8 protocol simulated at measured NuMI link | Measured or source-specific finite-message performance with synchronization, realistic coding, message classes, dropouts and channel feedback |
| Circular benchmark simulation | Figure 2 called an implementation check; published observed frame BER described and distinguished from generated counts | Digitize observed BER points with uncertainties or obtain run-level data to test the full decoder |
| Planetary prior art | Appendix A now cites the 2012 paper's mention of planetary blockage and explains full on-slot duration in rate inversions | Beam/receiver and body-transmission simulation for an actual ephemeris |
| Divergent applications | Submarine, finance and occultation treated as stress tests with named limitations | Independent source and detector design and end-to-end requirements for each use case |

## Prioritized next analysis

1. Select a neutrino source with an angular and energy-resolved flux and a
   published proton or muon input budget. Compute how many useful neutrinos
   reach a stated detector, rather than assigning divergence independently.
2. Fold three-flavor flavor evolution in a specified Earth density profile
   with detector material, geometry, event selection and backgrounds.
3. Simulate a complete source pulse schedule and finite messages. Compare OOK,
   repetition and one sparse-pulse option at equal **average** facility energy,
   recording packet reliability and acquisition time.
4. Validate the model against experimental event distributions where accessible,
   provide confidence intervals and repeat for a single use case.

## What the current packet example establishes

`empirical_message.py` samples uniform five-byte payloads with CRC-8/ATM.
The receiver knows slot boundaries externally and declares on when any event
appears across repeated on-symbol slots. It records undetected wrong payloads
as well as correctly accepted packets. The on-pulse energy column is **proton
beam energy to the target**; it is neither neutrino-carried energy nor facility
electricity. The original 2012 study used a convolutional code and 64-bit
synchronization word and reported a different decoded rate [1 in manuscript].
