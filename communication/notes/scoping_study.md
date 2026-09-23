# Direct through-Earth neutrino communication: reproducible scoping study

23 September 2026. Research note, **not a submission-ready manuscript**.

## Verified relevant publications and technical sources

| Work | What it contributes | URL |
| --- | --- | --- |
| Hallsjö (2026), *Locating nuclear-powered submarines with antineutrinos* | Companion study: count statistics, benchmark significance, detector deployment and geometry; passive low-energy reactor antineutrinos, not a communication beam. | https://arxiv.org/abs/2605.15642 |
| Hallsjö (2018), *Charged current quasi-elastic muon neutrino interactions in the Baby MIND detector* (PhD thesis) | Beam interaction reconstruction, detector efficiency, magnetic spectrometer and nuSTORM context. | https://theses.gla.ac.uk/41123/ |
| Stancil et al. (2012), *Demonstration of Communication using Neutrinos*, Mod. Phys. Lett. A 27, 1250077 | Actual one-way NuMI–MINERvA link; 0.81 events/on pulse, 0.1 decoded bit/s demonstration, 1.035 km including 240 m rock. | https://arxiv.org/abs/1203.2847 |
| Huber (2010), *Submarine neutrino communication*, Phys. Lett. B 692, 268–271 | Muon-storage-ring submarine downlink proposal; projected, not demonstrated. | https://arxiv.org/abs/0909.4554 |
| Learned, Pakvasa & Zee (2008), *Galactic Neutrino Communication* | Very long distance communication concepts and energy/beam choices, not an Earth-link benchmark. | https://arxiv.org/abs/0805.2429 |
| Fidalgo Prieto et al. (2022), *Submarine Navigation using Neutrinos* | Adjacent navigation and timing application, different message-delivery objective. | https://arxiv.org/abs/2207.09231 |
| DUNE Collaboration (2020), *Far Detector Technical Design Report, Volume I* | Multi-kt fixed detector reference; does not establish a deployable communication receiver. | https://arxiv.org/abs/2002.02967 |
| Hyper-Kamiokande Proto-Collaboration (2018), *Design Report* | Alternative large water-Cherenkov receiver design reference. | https://arxiv.org/abs/1805.04163 |
| nuSTORM Collaboration (2025), *Neutrinos from Stored Muons (nuSTORM)* | Controlled source and flux concepts; no assumed communication facility performance. | https://arxiv.org/abs/2505.06137 |
| Particle Data Group (2025), *Neutrino Cross Section Measurements* | Energy-dependent interaction inputs needed for a physics-grade spectrum integration. | https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf |
| Earth-core oscillation study (2021), *Neutrino Oscillations through the Earth's Core* | Flags matter/oscillation effects excluded from this first calculation. | https://arxiv.org/abs/2110.01148 |

This is a verified *working bibliography* from searches, not a claim to have exhaustively identified all literature. Bibliographies and forward citations of the first four works should be screened before submission. An earlier conversation incorrectly described Stancil as a Physical Review Letters paper; the publisher record linked from arXiv identifies Modern Physics Letters A.

## Research question and contribution

Can a fixed transmitter and receiver communicate directly through Earth at a stated decoded rate and error probability within defensible source-power and detector-size constraints? The novelty test is a common accounting framework that separates the measured demonstration, projected neutrino energy intercepted in an idealized beam, actual accelerator wall power, and message-level reliability. The companion Hallsjö (2026) paper studies an uncontrolled reactor antineutrino source; the source physics and optimization cannot simply be transferred to a controlled muon-neutrino beam.

## Reproducible preliminary model

Run `python model.py`. The outputs are `benchmark.csv` and `sensitivity.csv`. Python standard library only; fixed RNG seed.

The reproduction uses the *reported* 0.81 selected events per on-pulse from Stancil et al., not raw experimental events. For equally likely bits, zero background, perfect slot synchronization and one-or-more-event threshold, uncoded BER after pooling `n` independent pulses per bit is `0.5 exp(-0.81 n)`. The Monte Carlo simulates this Bernoulli consequence of a Poisson event count; it does **not** recreate Stancil's convolutional decoder, 64-bit sync sequence, 25-pulse supercycle or measured 0.1 decoded bit/s. The five-pulse prediction is 0.00871 BER, consistent in order with the paper's 99% correct bit statement.

For a *separate sensitivity exercise*, assume a 3 GeV monoenergetic collimated muon-neutrino beam, full-angle parameter represented by far-field half-angle theta=1 mrad, a 10 kt fiducial detector that is small relative to the beam footprint, selected efficiency 0.5, no background, no attenuation, no oscillation, perfect time alignment, equal 0/1 probabilities, and a representative charged-current cross section `0.67×10^-38 (E/GeV) cm²` per nucleon. This approximate cross section is **not** a precision value for 3 GeV. Neutrino kinetic energy is used as an idealized floor; the actual proton-beam and wall-plug power must exceed it and are not estimated here.

`A = pi (theta L)^2`, `p_selected = (M/m_u) sigma epsilon / A`, `lambda_req = -ln(2 BER)`, and `P_nu = R_bit lambda_req E_nu / p_selected`. The model intentionally takes flavor survival as unity. A real source has an energy spectrum and angle-energy correlations, and Earth oscillations may change selected charged-current rates. These optimistic results must **not** be described as real facility specifications or a universal physical lower bound.

At 1 raw bit/s and target uncoded BER 1%, the calculated beam-carried neutrino energies per second are:

| Chord | 10 kt, 1 mrad case | Status |
| --- | ---: | --- |
| 1,000 km | 0.976 MW | optimistic illustrative beam-energy requirement |
| 5,000 km | 24.401 MW | optimistic illustrative beam-energy requirement |
| 12,000 km | 140.549 MW | optimistic illustrative beam-energy requirement |

At fixed mass and divergence the result grows with the *square* of baseline. A tenfold narrower divergence reduces the idealized result by 100; whether that is feasible at the selected energy requires a beam-optics study. A tenfold larger detector reduces it by ten. The `sensitivity.csv` provides the 3 × 3 × 3 sweep. None of these figures include coding overhead, background, beam duty cycle, pointing loss, beam production efficiency, construction cost or receiver geometry; each tends to weaken the practical link relative to this optimistic case.

## Article draft: proposed structure and figures

**Provisional title:** *Direct Neutrino Communication Through the Earth: A Quantitative Feasibility Map*.

**Abstract (working):** We evaluate an idealized one-way through-Earth neutrino communication link using a Poisson counting receiver. A benchmark derived from the reported NuMI–MINERvA selected-event rate reproduces the expected uncoded on–off-keying error relation but does not reproduce its full decoding chain. A separate monoenergetic, far-field sensitivity model maps the neutrino beam energy required for an illustrative 1 bit/s, 1% uncoded-error link as a function of chord, beam divergence and detector mass. The large variation across these parameters motivates source-specific spectra, oscillation, detector-response and full-message studies before claims about feasible installations can be made.

1. Introduction: distinction among beam detection, raw bits, decoded payload and practical infrastructure; connect to prior Hallsjö passive-source analysis.
2. Prior work: Huber projection versus Stancil demonstration, adjacent navigation and detector studies.
3. Link model: chord geometry, beam spectrum, flavor propagation, detector acceptance, backgrounds and power definitions.
4. Decoder: OOK and synchronized Poisson slots; exact false-positive/false-negative curves, then sync and error-correction overhead.
5. Validation: published 0.81 event/pulse statistic, five-pulse error prediction, limits of that validation.
6. Sensitivity results: heatmaps of 1% BER beam energy over baseline/divergence and mass; data in `sensitivity.csv`.
7. Discussion: accelerator duty cycle, pointing, Earth matter and oscillations, engineering scale, comparison with alternatives for each use case.
8. Conclusion: explicit feasible or excluded regions only after a source-specific calculation.

**Figure plan:** (1) beam–Earth–detector geometry with chord length and beam footprint; (2) analytical versus Monte Carlo BER at the published selected-event operating point; (3) required neutrino-carried MW versus chord for three divergence values and several masses; (4) credible facility power versus detector mass once source conversion and backgrounds are modelled. Figure 4 cannot be populated responsibly from the current inputs.

## Remaining work before submission

1. Extract beam energy spectra, protons-on-target, duty cycle and actual far-field angular distributions from specific facility design reports. Define a physically achievable transmitter at each chord.
2. Integrate an energy-dependent cross section from PDG or a validated interaction model; propagate oscillations through an Earth density model and quantify beam absorption where relevant.
3. Specify receiver material/shape and detector projected area; apply energy-dependent efficiency and backgrounds. The simple mass/beam-area expression is valid only when the beam footprint is appreciably larger than the detector.
4. Implement nonzero-background optimal likelihood threshold, time acquisition, realistic synchronization, message framing, coding, outages and confidence intervals. Benchmark complete payload throughput against Stancil's published protocol, avoiding confusion between raw and decoded rates.
5. Sensitivity and uncertainty: beam energy spectrum, divergence, flavor survival, reconstruction efficiency, backgrounds, power conversion, 1/10/40 kt masses and multiple rates. Compare fixed site links and mobile receivers separately.
6. Review all references forward and backward, deposit model code and parameter manifest, and have the physics and communications calculations independently checked before any publication claim.

## Revision 0.2: detector, submarine and latency extensions

`analysis/tradeoffs.py` and the CSV files `data/power_detector_tradeoff.csv` and `data/latency_example.csv` invert the same conditional model. At a 5,000 km chord and 1 mrad, a one-second on symbol has `P_neutrino(MW) × M(kt) ≈ 244 MW·kt` at ideal 1% raw BER; 10 kt and 1 ms imply 24.4 GW carried by neutrinos. The conventional-fiber comparison uses a hypothetical geometry and Corning 1550 nm group index 1.4677, not a measured financial route. See the manuscript for the full assumptions and literature context.
