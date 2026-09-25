# Research proposal draft

## From sparse neutrino detections to reliable through-Earth communication

### Proposed linked studies

1. **Paper II:** *End-to-End Simulation of a Neutrino Communication Channel Through the Earth*
2. **Follow-on information-theory study:** *The Neutrino Channel: Capacity and Coding in the Few-Event Regime*

**Status:** Research proposal draft. The source and receiver configuration for Paper II remains to be selected; no new source-derived performance result is claimed here.

## 1. Summary

This proposal develops a physics-grounded model of direct neutrino communication through the Earth, then uses that model to study how much information can be delivered when a receiver observes only a few neutrino interactions per symbol. The work is motivated by two complementary strands of prior neutrino research: a study of passive antineutrino detection from nuclear-powered submarines, and a scoping study of accelerator-generated neutrino communication through the Earth.

The prior detection work examined whether sparse, background-limited events from an uncontrolled reactor source could support a detection task. The communication scoping paper considered a controlled accelerator beam and message transmission, but its long-baseline rate calculation was conditional on simplified assumptions. The proposed work connects these questions through a common event-level framework while preserving their physical differences. Reactor-antineutrino rates and detector assumptions will not be transferred to the accelerator communication case.

Paper II will derive signal events from an explicit source spectrum, Earth propagation, interaction cross sections, detector geometry and response, and backgrounds. It will pass the resulting time-dependent event stream through a declared, reproducible communication protocol. A separate follow-on study will use that physics-derived channel to compare channel capacity and finite-block coding strategies in the few-event regime. Together, the studies will distinguish a beam that produces detectable events from a link that delivers reliable payload at a stated rate and energy cost.

## 2. Prior work and rationale

### 2.1 Passive antineutrino detection

In *Locating nuclear-powered submarines with antineutrinos*, Hallsjö studied passive detection of an uncontrolled reactor-antineutrino source in a sparse-event setting. Its analysis emphasizes geometry, event counts, background, operating points, detector mass, and deployment configuration. It is a useful methodological predecessor for handling low counts and decision thresholds. It addresses detection of a source, however, rather than transmission of a chosen message with a controlled beam. Its reactor spectra, rates, and detector assumptions therefore cannot be used as inputs to the proposed communication link. The work is available as arXiv:2605.15642.

### 2.2 Accelerator-generated neutrino communication

The working manuscript *Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity* separates a NuMI–MINERvA-calibrated packet simulation from an idealized far-field sensitivity calculation. The packet exercise uses published event-rate and timing inputs but is not a reproduction of the experiment's decoder or a measurement of those simulated packet outcomes. The far-field result is conditional on simplified source, propagation, receiver, and background assumptions. It identifies the central next step: replace the assumed event rate with a source- and detector-derived rate.

Stancil *et al.* demonstrated a neutrino communication link over a 1.035 km baseline, including 240 m of earth. This is an experimental proof of principle at its reported configuration. It does not establish a long-baseline source–receiver design or validate the proposed end-to-end model.

### 2.3 Research gap

The missing bridge is a documented mapping from accelerator operation to useful, decoded information:

\[
\text{source spectrum}
\rightarrow \text{Earth propagation}
\rightarrow \text{selected event stream}
\rightarrow \text{message reliability and delivery cost}.
\]

Without that bridge, a Poisson channel can be studied mathematically, but its signal and background rates remain assumptions. Conversely, a detailed event-rate calculation alone does not determine which modulation and coding scheme delivers the most useful information under timing, energy, and reliability constraints.

## 3. Research objectives and questions

The project has five objectives:

1. Select one accelerator-source and detector configuration with flux, geometry, operating, and response information sufficient for a reproducible benchmark.
2. Compute the energy- and time-dependent neutrino flux at the receiver, including three-flavor propagation through a stated Earth density and electron-fraction profile.
3. Derive selected signal and background events using specified interaction channels, detector material, projected geometry, efficiency, and live-time gates.
4. Evaluate end-to-end message performance with explicit synchronization, framing, decoding, and failure accounting.
5. Characterize the few-event channel's capacity and compare practical finite-block codes under matched energy, pulse-rate, and reliability budgets.

The central questions are:

- How does a specified accelerator source produce selected event rates across a stated Earth chord and detector configuration?
- How do oscillations, energy-angle correlations, detector selection, and background alter the event statistics available to a decoder?
- At those physics-derived event rates, what payload rate, message success probability, latency distribution, and facility-input energy per delivered bit can a concrete protocol achieve?
- How do simple OOK, sparse-pulse or pulse-position modulation, and finite-block codes compare when the receiver observes only a few events per symbol?

## 4. Proposed work

### Work package 1: Source and benchmark definition

Choose a source with a traceable energy- and angle-resolved flux or a documented simulation input, a stated beam timing and duty cycle, and a defensible primary-beam or facility-input budget. Define the transmitter and receiver coordinates, chord length, source normalization, flux surface, detector material, geometry, fiducial mass, and event-selection response.

A source data sufficiency gate will be applied before numerical results are described as source-derived. If only a near-detector flux is available, any transformation to the communication geometry must be documented and validated. A free beam-divergence parameter will not be presented as an achieved accelerator specification.

### Work package 2: Propagation and interaction model

Propagate the source spectrum along the specified chord using three-flavor matter evolution with a named Earth density model, such as PREM, and an explicit electron-fraction profile. Test vacuum, constant-density, probability-normalization, and numerical convergence limits. Include energy-changing interactions during transit if relevant for the selected spectrum; otherwise quantify why they can be neglected.

Fold the propagated flux through energy-dependent cross sections and a receiver response. The geometry and target convention will prevent double counting of acceptance. Versioned software, parameter choices, and uncertainty inputs will be recorded.

### Work package 3: Event stream and baseline protocol

Generate signal and background events in actual acquisition windows using the accelerator time structure. Specify whether timing is externally supplied or acquired from the neutrino stream. Implement a transparent baseline protocol with modulation, synchronization, framing, error detection, and decoder behavior.

Report intermediate event counts separately from complete-message outcomes. Core metrics will include false acquisition, packet success, detected and undetected errors, useful payload bits per wall-clock second, latency quantiles, and facility joules per successfully delivered bit. Neutrino-carried energy and accelerator electrical/input energy will remain separate quantities.

### Work package 4: Few-event capacity and coding study

Use the channel characterized in Paper II as the physical input to the follow-on study. Begin with a Poisson count model,

\[
K_i\mid x_i\sim\operatorname{Poisson}(b_i+s_i x_i),
\]

where \(s_i\) and \(b_i\) are taken from the source-, propagation-, detector-, and background model or its uncertainty ensemble. Calculate information-theoretic reference limits for declared input constraints, then compare finite-length OOK, repetition with soft combining, pulse-position modulation, and selected error-correcting codes.

Compare schemes under matched average source or facility energy and the same timing, background, acquisition, and packet-success constraints. Include synchronization and framing overhead in delivered payload metrics. Capacity bounds will be labeled as theoretical bounds, while simulation results will be labeled by the decoder and finite-message protocol actually implemented.

### Work package 5: Uncertainty, validation, and reproducibility

Separate uncertainty from source flux, oscillation parameters and Earth density, cross sections, detector response, backgrounds, timing, and Monte Carlo sampling. Retain correlations where available and identify any unsupported uncertainty component.

Validation will proceed from numerical limiting cases to source-flux reproduction, detector-rate comparisons for matching configurations, analytical communication-channel checks, and comparison with published NuMI–MINERvA quantities where reconstructible. The published demonstration will be a reference, not proof that the full proposed system has been validated.

Code, input-data provenance, environment and software versions, run configurations, generated outputs, and manuscript sources will be stored with the project. Reproduction instructions and fixed random seeds will be provided for stochastic calculations.

## 5. Relationship between the papers

| Study | Main question | Main contribution | Boundary |
| --- | --- | --- | --- |
| *Locating nuclear-powered submarines with antineutrinos* | Can a sparse, uncontrolled reactor-antineutrino signal support detection under stated geometry and background? | Passive-source detection framework and operating-point analysis | Does not model controlled message transmission |
| *Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity* | What does the experimental benchmark and an idealized sensitivity map imply for through-Earth communication? | Scoping calculations and explicit identification of missing source, propagation, receiver, and protocol inputs | Long-baseline rates remain conditional |
| Paper II, *End-to-End Simulation…* | What event stream follows from a specified accelerator source, Earth path, and receiver? | Physics-grounded channel model and complete baseline-message simulation | Does not claim optimal coding or facility feasibility without validated inputs |
| Follow-on, *The Neutrino Channel…* | How much information can be delivered over that few-event channel? | Capacity analysis and matched-budget comparison of finite-block coding and modulation | Uses Paper II's channel model; does not replace its source and detector calculations |

The two earlier studies provide complementary experience with sparse-event modeling: passive detection from reactor antineutrinos and controlled communication with accelerator neutrinos. Their source physics and task objectives remain distinct. The proposed sequence joins them at the statistical and reproducibility level, not by treating their rates as interchangeable.

## 6. Expected outcomes

The project is expected to produce:

- a versioned source-to-receiver event-rate calculation with traceable inputs;
- a reproducible end-to-end baseline decoder with message-level metrics;
- uncertainty-aware channel parameters for specified source, path, detector, and background conditions;
- capacity bounds and finite-block coding comparisons for a physically grounded few-event channel;
- tables and figures that clearly distinguish measured benchmarks, model outputs, theoretical limits, and conditional sensitivities.

The project will report numerical feasibility only for the declared benchmark and assumptions. A computed receiver event rate does not itself establish a practical installation. A measured short-baseline link does not establish a useful global link.

## 7. Risks and responses

**Insufficient public source information.** A facility may not publish a suitable far-field energy-angle distribution. The project will document this limitation and either select a source with sufficient inputs or clearly identify any required source simulation; it will not substitute a free divergence assumption.

**Detector response is unavailable.** A design report may lack the efficiency curve or background for the proposed operating mode. Use a transparent parameterized response only as a sensitivity study, and label results as conditional until a defensible response exists.

**Computational model scope expands.** Full beamline, transport, and detector simulation may exceed the initial study. Use staged validation, begin with flux tables and response parametrizations, and add complexity only when it changes the inference.

**Potential overlap between manuscripts.** Keep Paper II centered on deriving the physical event channel and demonstrating one explicit end-to-end protocol. Keep capacity and coding theory in the follow-on study, while ensuring Paper II exports the channel distributions and uncertainty inputs that the coding study needs.

## 8. Draft publication sequence

1. Complete and review Paper II's source, propagation, detector, background, and baseline protocol model.
2. Publish or archive the reproducible Paper II inputs, code, and outputs with explicit limitations.
3. Use that model as the physical channel for the capacity and finite-block coding analysis.
4. Consider a later network study only after a single-link channel has defensible source-derived performance.

## References

1. S.-P. Hallsjö, “Locating nuclear-powered submarines with antineutrinos,” arXiv:2605.15642 (2026). https://arxiv.org/abs/2605.15642
2. S.-P. Hallsjö, “Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity,” working manuscript v0.7, September 2026, repository snapshot at commit 773faadeb00af5b8e0f162feb3d3491327bb22af: https://github.com/patha325/neutrinoArticles/tree/773faadeb00af5b8e0f162feb3d3491327bb22af/communication/manuscript
3. D. D. Stancil *et al.*, “Demonstration of Communication using Neutrinos,” *Modern Physics Letters A* **27**, 1250077 (2012). https://doi.org/10.1142/S0217732312500770 ; https://arxiv.org/abs/1203.2847
4. A. M. Dziewonski and D. L. Anderson, “Preliminary reference Earth model,” *Physics of the Earth and Planetary Interiors* **25**, 297–356 (1981). https://doi.org/10.1016/0031-9201(81)90046-7
5. C. Andreopoulos *et al.*, “The GENIE Neutrino Monte Carlo Generator,” *Nuclear Instruments and Methods in Physics Research A* **614**, 87–104 (2010). https://doi.org/10.1016/j.nima.2009.12.009
