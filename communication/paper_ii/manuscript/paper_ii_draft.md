# Paper II: End-to-End Simulation of a Neutrino Communication Channel Through the Earth

*Working manuscript draft; author details omitted from this repository copy.*

> **Status: methods and analysis plan.** Numerical long-baseline results are intentionally withheld until source flux, propagation, detector response, background, and decoder inputs are fixed and validated. Nothing in this draft establishes engineering feasibility.

## Abstract

The first paper in this sequence separated a NuMI-calibrated finite-message simulation from an idealized far-field sensitivity calculation. Its long-baseline event rate was conditional on a monoenergetic beam, a free angular-divergence parameter, unit flavor survival, zero background, and a simplified mass-only receiver. Here we formulate an end-to-end simulation in which the expected selected-event rate is derived from an energy- and angle-resolved source yield, propagation through the Earth, interaction probabilities in a stated detector material and geometry, and event-selection response. The resulting time-tagged signal and background events are passed through a finite-message protocol with explicit framing, synchronization, and coding. The calculation will report selected events per transmitted slot, packet success and undetected-error probabilities, decoded payload rate, latency distributions, and facility-input energy per successfully delivered payload bit. The primary calculation will use a declared accelerator-source configuration and documented flux inputs; the NuMI–MINERvA demonstration will serve as a validation reference only to the extent that published inputs permit. Three-flavor matter evolution will use a stated density and electron-fraction profile, while interaction and detector uncertainties will be propagated separately. This draft specifies the model and validation gates; source-specific results remain pending.

**Keywords:** neutrino communication; Earth propagation; accelerator neutrino flux; detector response; finite-block coding; Poisson channel; reproducible simulation.

## 1. Introduction

A neutrino beam offers a direct path through matter, but weak interaction makes reception inefficient. The relevant question is not whether a neutrino can cross the Earth; it is whether a specified transmitter and receiver can deliver a finite message with stated reliability, latency, and energy cost. An assumed detector event mean can make a channel calculation mathematically complete while leaving its most consequential physical inputs unspecified.

Paper I, *Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity*, established two deliberately separate reference calculations. One used the selected-event mean and timing reported for the NuMI–MINERvA communication demonstration to simulate a different repeated-OOK/CRC protocol. The other mapped conditional far-field sensitivity under a monoenergetic approximation and ideal receiver assumptions. Neither calculation supplied a source-derived long-baseline spectrum, Earth-propagated flavor composition, detector response, or background model. Its power figures were conditional neutrino-carried power during on-slots, not accelerator electrical power or a facility requirement.

This paper addresses that gap by making the event-rate parameter an output of a chain of physical and communication models. The chain begins with an accelerator and focusing configuration and its energy-angle neutrino yield. It follows the beam through a specified Earth chord, applies flavor evolution and any material transport effects, folds the resulting flux through interaction cross sections and a geometrical detector response, and generates time-tagged signal and background events. A receiver then performs acquisition, synchronization, framing, decoding, and error detection. The final reported quantity is successful decoded payload per elapsed time, with message-level reliability and energy accounting.

The objective is methodological and quantitative. The free beam-divergence parameter in Paper I will not be treated as a source specification. Accelerator wall power will not be inferred from neutrino-carried energy, and event-count capacity will not be equated with successfully decoded payload. Results will be conditional on the declared source, detector, path, and protocol, and their uncertainties will be reported with those conditions.

Stancil *et al.* reported a one-way message transmitted over a 1.035 km baseline, including 240 m of earth, using the NuMI beam and MINERvA detector [1]. That result is a useful experimental reference, not a long-baseline validation dataset. The published selected-event mean does not by itself validate a new source model, Earth-propagation calculation, or decoder.

## 2. Scope and reported quantities

The system boundary is

\[
\text{accelerator and target}
\rightarrow \text{focusing and decay region}
\rightarrow \frac{d^2N_{\nu_\alpha}}{dE\,d\Omega}
\rightarrow \text{Earth propagation}
\rightarrow \text{detector interactions and selection}
\rightarrow \text{event record}
\rightarrow \text{message decoder}.
\]

Each benchmark will state the source facility and operating mode; primary-beam energy and intensity; usable exposure and duty cycle; flux table and normalization; source-detector baseline and chord; oscillation and Earth model; detector material, active mass, projected area and position; selection response; background and acquisition gate; modulation, frame, code and decoder; and facility-input energy definition.

We distinguish:

1. **Incident neutrino fluence:** neutrinos crossing the receiver plane per area, energy, flavor, and time.
2. **Selected-event mean:** expected accepted signal events in a declared communication slot.
3. **Raw symbol rate:** modulation symbols transmitted per elapsed second, including inactive source periods as specified.
4. **Decoded payload rate:** application payload bits in successfully accepted messages divided by wall-clock acquisition time.
5. **Facility energy per delivered bit:** accelerator and beam-facility input energy divided by successfully delivered payload bits. Neutrino-energy flux is not a substitute.
6. **Message reliability:** probability that a complete framed payload is accepted correctly, including detected and undetected errors, synchronization loss, and acquisition failure.

The Shannon capacity of a simplified counting channel is a theoretical reference, not a measured or simulated finite-message throughput.

## 3. Source spectrum and geometric transport

Let source configuration \(s\) specify the primary-beam species and energy, target, focusing system, decay region, and operating mode. The source model supplies

\[
Y_{\alpha}(E,\Omega\mid s)
= \frac{d^2N_{\nu_\alpha}}{dE\,d\Omega}
\quad [\mathrm{neutrinos\;POT^{-1}\,GeV^{-1}\,sr^{-1}}],
\]

or an equivalent flux table at a documented reference surface. The table must identify flavor, polarity, normalization, binning, geometric surface, and whether oscillations or detector effects are already applied. Conversion from protons on target to facility exposure is recorded separately.

The input spectrum will not be replaced by a freely chosen monoenergetic beam or divergence. If public source data provide only a near-detector flux rather than an angle-resolved emission yield, transformation to the target baseline must be documented and validated; a near-detector spectrum cannot be treated as a universal far-field angular distribution.

For a point source and a receiver small compared with its distance, the unoscillated fluence contains the geometric factor \(1/L^2\). For an extended decay region or finite receiver plane, the calculation integrates over source and receiver coordinates:

\[
d\Phi_{\alpha}^{0}
= \int_{V_s}d^3x_s\,q_{\alpha}(E,\Omega;\mathbf{x}_s)
\frac{\cos\vartheta}{r^2}\,dA .
\]

Here \(q_{\alpha}\) is the emitted yield density under the selected normalization, \(r\) is the source-element to receiver-element distance, and \(\vartheta\) is incidence angle on the receiver plane. The implementation will state when it uses the point-source approximation and verify the inverse-square limit.

The source is a decision gate for numerical results. A baseline configuration will be selected from a facility with documented flux inputs suitable for the chosen geometry. Input tables and provenance will be recorded before any result is labeled source-derived.

## 4. Propagation through the Earth

### 4.1 Chord and density profile

Endpoint coordinates define a straight-line chord through a spherical Earth with a stated radius. The path coordinate is \(x\in[0,L]\); radial distance determines local density. The initial reference model will be PREM [2], with an explicit electron fraction \(Y_e(r)\). Crust composition near the endpoints will be stated because a radial average does not describe every local site.

For an incident flavor \(\alpha\), let \(P_{\alpha\beta}(E,L)\) denote the probability of detection as flavor \(\beta\). In natural units, coherent three-flavor evolution is

\[
i\frac{d}{dx}\nu_f(x)=
\left[
\frac{1}{2E}U
\begin{pmatrix}
0&0&0\\
0&\Delta m^2_{21}&0\\
0&0&\Delta m^2_{31}
\end{pmatrix}
U^\dagger+
\begin{pmatrix}
V(x)&0&0\\
0&0&0\\
0&0&0
\end{pmatrix}
\right]\nu_f(x),
\qquad V(x)=\sqrt{2}G_FN_e(x).
\]

For antineutrinos, the matter potential changes sign and the mixing matrix is complex-conjugated. Numerical integration will use a documented set of oscillation parameters and will distinguish parameter uncertainty from numerical integration error. The evolution operator will be checked for unitarity in the absence of absorption and against vacuum and constant-density limits.

At energies where interactions during transit are non-negligible, flavor evolution alone is insufficient. The transport stage will include energy-changing neutral-current scattering and relevant charged-current absorption or regeneration through a transport kernel, or quantitatively show that these effects are negligible over the chosen energy and chord range. The approximation will be selected from the benchmark spectrum rather than presumed in advance.

## 5. Detector response and selected event rate

Let \(\Phi_\beta(E,\mathbf{x},t)\) be the flux at the receiver after geometric transport and propagation. The selected-event rate for interaction channel \(c\) is

\[
\frac{dN_{\mathrm{sig}}}{dt}
=
\sum_{\beta,c}\int dE\int_{A_{\mathrm{proj}}}dA\,
\Phi_\beta(E,\mathbf{x},t)
\,n_{T,c}(\mathbf{x})
\,\sigma_{\beta c}(E)
\,\epsilon_{\beta c}(E,\mathbf{x},\mathrm{event}),
\]

where \(n_{T,c}\) is the target density appropriate to the cross-section convention, \(\sigma_{\beta c}\) is the interaction cross section, and \(\epsilon\) includes trigger, reconstruction, containment, and selection. Geometry must not be counted twice: flux is either integrated over the actual detector volume or combined with an effective-area/volume response.

The detector benchmark will name material, geometry, fiducial volume, and readout model. A mass-only receiver is insufficient when beam footprint, projected area, containment, or event selection depends on energy or interaction topology. The initial interaction calculation may use tabulated cross sections with published uncertainties. A generator such as GENIE may be used for final-state and detector-level response after its version, tune, target model, and configuration are pinned [3]. Generator output is not experimental data.

The expected selected signal in slot \(j\) is

\[
\mu_{s,j}=\int_{t_j}^{t_j+\Delta t_j}dt\,
\frac{dN_{\mathrm{sig}}}{dt}.
\]

Slot means need not be equal: accelerator pulse structure, beam uptime, timing phase, and source interruptions will be represented at appropriate granularity.

## 6. Backgrounds and event stream

For slot \(j\), an initial count model is

\[
N_j\sim\operatorname{Poisson}(b_j+x_j\mu_{s,j}),
\qquad x_j\in\{0,1\},
\]

where \(b_j\) is expected background in the same detector, selection, and live-time gate as the signal. Background components will be listed separately when timing or event distributions differ. If event-level data are available, decoding will use event times and reconstructed observables rather than reducing all information to a count; otherwise the limitation of count-only decoding will be explicit.

The background model must correspond to a specified receiver location and veto, shielding, depth, and acquisition window. Zero background may be used only as a labeled diagnostic limit. Signal and background parameter uncertainties will be reported separately from finite Monte Carlo uncertainty.

## 7. Modulation, synchronization, and decoding

The transmitter will use a declared modulation and source schedule. The first benchmark will compare OOK with at least one energy-matched sparse-pulse or pulse-position scheme: equal peak power does not imply equal average source energy. A receiver without external slot timing must acquire synchronization from the received stream or from a declared independent timing channel.

For count-only OOK with known slot timing, the count likelihood ratio is

\[
\Lambda_j(n_j)=
\frac{P(n_j\mid x_j=1)}{P(n_j\mid x_j=0)}
=\exp(-\mu_{s,j})
\left(1+\frac{\mu_{s,j}}{b_j}\right)^{n_j},
\]

for \(b_j>0\), with the zero-background limit handled separately. This is a symbol detector, not a complete protocol. The implemented receiver will account for preamble detection, frame boundaries, payload and check bits, decoder failure, undetected errors, retransmission or acknowledgement policy if used, and elapsed idle time.

Principal outputs are complete-message success probability; detected and undetected error probabilities; useful payload bits per wall-clock second; latency quantiles; and facility joules per successful payload bit. Raw event and symbol rates are intermediate quantities only.

## 8. Uncertainty, validation, and reproducibility

Uncertainty sources will be partitioned into source-flux normalization and shape; parent production and focusing where available; oscillation parameters and Earth profile; interaction cross sections and final-state modeling; detector efficiency and selection; background; beam timing and uptime; and Monte Carlo sampling. Correlations will be retained where supported by source inputs. Results will identify uncertainties not represented.

Validation proceeds in layers:

1. **Numerics:** vacuum and constant-density oscillation limits; probability normalization; stability under path-step and energy-bin refinement; geometric inverse-square behavior where applicable.
2. **Source:** reproduce published flux projections or tables at their documented reference geometry before transforming to the receiver geometry.
3. **Interactions and detector:** compare folded rates with published predictions for a matching beam, target, and selection where available.
4. **Communication:** reproduce analytical Poisson OOK error probabilities and test framing, synchronization, and decoding on deterministic fixtures.
5. **Experimental reference:** compare only NuMI–MINERvA observables that can be reconstructed from published information. This is not a reproduction of the original experiment unless its event selection, timing, and decoder are actually reproduced.

Code, input manifest, software versions, run configuration, generated tables, and figures will be committed alongside the manuscript. Each result will be regenerable from a clean checkout using documented commands and fixed random seeds for stochastic sampling.

## 9. Results (pending model lock)

No source-derived end-to-end numerical result is reported in this draft. This section will be populated only after source and detector inputs are selected, transport and response stages are implemented, and validation gates in Section 8 pass.

Planned outputs are:

- energy and angle distributions at source and receiver, before and after flavor propagation;
- selected signal and background events per slot for each transmitted symbol;
- message-level reliability, decoded payload rate, latency distribution, and facility energy per successful payload bit;
- uncertainty intervals and sensitivity contours for named operating points;
- modulation and coding comparisons at matched source or facility energy.

A table of required numerical fields and provenance will accompany the implementation. Missing inputs will be marked unavailable rather than replaced with an unreferenced nominal value.

## 10. Discussion and limitations (pending results)

The eventual discussion will separate simulation limitations from installation limitations. A source simulation does not establish that a beamline can be constructed or operated in the required orientation. A large underground detector is not automatically a communication receiver, and the NuMI–MINERvA demonstration does not validate a different baseline or protocol. Facility power, duty cycle, pointing, receiver deployment, backgrounds, and message deadlines must be addressed before a use case can be called feasible.

This draft defines the calculation and evidentiary standard rather than reporting a feasibility conclusion. Application examples from Paper I will not be upgraded into feasibility claims by applying a propagation factor to prior conditional numbers.

## 11. Conclusion (working)

Paper II is designed to replace a conditional event-rate parameter with a reproducible source-to-decoder calculation. Its central object is the selected, time-dependent event distribution derived from source spectrum, Earth propagation, detector interactions and response, and background. The intended endpoint is complete-message delivery performance and facility-input energy per successful payload bit. Numerical conclusions await selection and validation of a concrete source and receiver configuration.

## References

[1] D. D. Stancil *et al.*, “Demonstration of Communication using Neutrinos,” *Modern Physics Letters A* **27**, 1250077 (2012). https://doi.org/10.1142/S0217732312500770 ; https://arxiv.org/abs/1203.2847

[2] A. M. Dziewonski and D. L. Anderson, “Preliminary reference Earth model,” *Physics of the Earth and Planetary Interiors* **25**, 297–356 (1981). https://doi.org/10.1016/0031-9201(81)90046-7

[3] C. Andreopoulos *et al.*, “The GENIE Neutrino Monte Carlo Generator,” *Nuclear Instruments and Methods in Physics Research A* **614**, 87–104 (2010). https://doi.org/10.1016/j.nima.2009.12.009 ; https://arxiv.org/abs/0905.2517

[4] T. Ohlsson, “Neutrino oscillations with three flavors in matter: Applications to neutrinos traversing the Earth,” *Physics Letters B* **486**, 19–26 (2000). https://doi.org/10.1016/S0370-2693(00)00741-2 ; https://arxiv.org/abs/hep-ph/9912295

[5] P. B. Denton and R. Pestes, “Neutrino Oscillations through the Earth's Core,” *Physical Review D* **104**, 113007 (2021). https://doi.org/10.1103/PhysRevD.104.113007 ; https://arxiv.org/abs/2110.01148

[6] Particle Data Group, “Neutrino Cross Section Measurements,” *Review of Particle Physics* (2025). https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf
