# Direct Neutrino Communication Through the Earth: An Initial Quantitative Feasibility Map

**Sven-Patrik Hallsjö**  
**Working manuscript, version 0.1 — 23 September 2026**

> **Status.** This draft combines a published experimental benchmark with a deliberately simplified sensitivity calculation. It has not established the feasibility of a particular accelerator–detector installation. Values described as “beam power” are energy carried by neutrinos; they are neither proton beam power nor facility electrical power. The source-specific and propagation calculations identified in Section 6 are necessary before submission.

## Abstract

Neutrinos can cross substantial amounts of matter, making a direct communication path through the Earth physically possible. Their small interaction probability, however, transfers the difficulty to the transmitter and receiver. We construct an initial quantitative map of that trade-off. First, we reproduce the uncoded on–off-keying error probability implied by the NuMI–MINERvA communication demonstration, using its reported mean of 0.81 selected signal events per transmitted pulse. Pooling five ideal independent pulses gives a predicted bit-error probability of 0.00871; a seeded 200,000-bit simulation gives 0.00893. Second, we examine an independent, idealized long-baseline model with a monoenergetic 3 GeV beam, specified angular spread, detector mass and event-selection efficiency. For a 10 kt receiver, 1 mrad beam half-angle and 50% efficiency, the neutrino-carried power associated with 1 raw bit s⁻¹ at 1% uncoded bit error is approximately 0.98, 24.4 and 140.5 MW for 1,000, 5,000 and 12,000 km chords, respectively. These are sensitivity results conditional on optimistic assumptions, not predictions for an existing facility. The strongest dependence is on angular spread and distance; the remaining determinants include the achievable beam spectrum, flavor evolution, detector response, background, timing and source conversion efficiency. A practical feasibility claim requires integrating these effects and measuring decoded payload throughput.

**Keywords:** neutrino communication; through-Earth links; Poisson channel; accelerator beam; detector sensitivity; feasibility.

## 1. Introduction

Communication through rock or across a large terrestrial chord ordinarily requires an indirect route. Neutrino beams offer an unusual direct path because most neutrinos traverse matter without interacting. The same property makes reception difficult: very few transmitted particles generate identifiable events in a finite detector. Consequently, a beam that can be detected statistically need not carry a useful message at an acceptable error rate or cost.

The possibility has passed an experimental proof-of-principle test. Stancil *et al.* transmitted an encoded message using the NuMI beam and MINERvA detector over 1.035 km, including 240 m of earth. They report a decoded rate of 0.1 bit s⁻¹ and a 1% bit error rate [1]. That result establishes a link at its particular source, geometry, detector and decoding protocol; it does not answer whether a regional or global direct link could be operated at a useful rate.

Huber studied a one-way link to a submerged submarine using a high-energy neutrino beam from a muon storage ring [2]. Learned, Pakvasa and Zee examined signalling on galactic scales [3]. A later submarine-navigation proposal addresses position, navigation and timing rather than reliable delivery of an arbitrary message [4]. These studies motivate a clear separation between demonstrated transmission and projections under particular technologies.

This study also follows Hallsjö's analysis of locating nuclear-powered submarines by their emitted antineutrinos [5]. That work treated an uncontrolled, comparatively low-energy reactor source as a detection problem and examined sparse counts, background, geometry and detector deployment. Here the question is different: for a controlled and modulated source, what signal count is needed to decode bits, and how does an idealized source requirement change with baseline, angular spread and receiver size? The detector and interaction expertise documented in Hallsjö's Baby MIND thesis informs the next, detector-specific stage of this work [6]. Reactor-antineutrino event rates or detector assumptions from [5] are not transferred to the high-energy communication example.

Our immediate aim is a transparent *scoping calculation*. We report what the stated assumptions imply and identify which missing inputs prevent an engineering conclusion.

## 2. Link definition and communication metric

Consider a stationary source that sends on–off-keyed symbols toward a stationary receiver through a chord of length `L`. In an ideal synchronized slot, an on symbol yields a selected signal count with Poisson mean `s`; an off symbol has no signal. A background count with mean `b` may be present in either slot. Conditional on transmitted symbol `x ∈ {0,1}`, the count model is

`K | x ~ Poisson(b + s x)`.                                                    (1)

For equally probable input symbols, `b = 0`, and a threshold of one observed event, false-positive errors vanish and an on symbol is missed with probability `exp(−s)`. Therefore

`BER_uncoded = ½ exp(−s)`.                                                     (2)

This equation measures uncoded, synchronized *raw* bits. A delivered payload rate must additionally account for symbol slots, synchronization, framing, coding, source duty cycle and any failed frames. An actual background `b > 0` introduces both false-positive and false-negative errors; neither Eq. (2) nor a threshold of one should then be assumed optimal.

We choose a target `BER_uncoded = 0.01`, so the mean signal requirement under the stated ideal conditions is

`s_required = −ln(2 × 0.01) = 3.912 selected events per on symbol`.            (3)

For 1 raw bit s⁻¹ and equal on/off probabilities, the source emits an on pulse in half the symbol slots on average. To avoid ambiguity, the calculation below treats the *on-slot* signal requirement as a continuous per-bit source budget: it effectively specifies the peak beam energy during an on symbol, or a conservative full-duty equivalent at 1 symbol s⁻¹. The **time-averaged transmitted energy for equiprobable data would be half of the tabulated full-duty-equivalent values**, before adding framing and other overhead. This convention must be held fixed when comparing source power. 

## 3. Methods

### 3.1 Published experimental benchmark

Stancil *et al.* report an average of approximately `0.81` selected muon events for a beam-on pulse in the reduced-intensity communication run [1]. Their measurement includes charged-current interactions in upstream rock whose muons enter MINERvA, as well as a smaller in-detector component. Consequently, it is incorrect to turn that count directly into a fiducial-target interaction probability for the detector alone.

For the limited purpose of validating Eq. (2), we assume independent pulses and zero background, and pool `n` pulses for each bit, giving `s = 0.81 n`. We compare the analytical error probability with a fixed-seed simulation of 200,000 equally likely bits for each `n`. A simulated on bit is missed with probability `exp(−0.81n)`; an off bit is always decoded correctly under this model. This is a simulation of the Poisson zero-count event, not a reconstruction of MINERvA data. Stancil *et al.* used a synchronization sequence, convolutional error correction and repeated frames; this study does not reimplement those procedures [1].

### 3.2 Illustrative far-field receiver model

For a separate link sensitivity study, we assume an approximately circular far-field neutrino footprint with characteristic half-angle `θ`. At distance `L`, its area is `A = π(θL)²`. If the footprint is appreciably larger than the receiver and illumination is sufficiently uniform, an approximate selected-event probability per neutrino directed into that footprint is

`p_sel ≃ (M/m_u) σ_CC(E) ε / [π(θL)²]`,                                     (4)

where `M` is fiducial target mass, `m_u` is the atomic mass unit expressed in grams, `σ_CC` is the per-nucleon charged-current cross section and `ε` subsumes trigger and event-selection efficiency. This expression assumes a target small compared with the beam footprint and does not specify a receiver shape. It cannot be extrapolated to an arbitrarily narrow beam or a detector larger than the footprint.

We set `E = 3 GeV`, `σ_CC = 0.67 × 10⁻³⁸(E/GeV) cm²` per nucleon as an approximate illustrative parametrization, `ε = 0.5`, and vary `M = 1, 10, 40 kt`, `θ = 0.1, 1, 10 mrad`, and `L = 1,000, 5,000, 12,000 km`. The cross-section approximation is not a precision calculation at 3 GeV; a subsequent version must use energy-dependent inputs and uncertainties from the Particle Data Group [7]. We assume a single selected flavor at the receiver, no oscillation loss, no absorption, no off-slot background, perfect pointing and timing, and an unconstrained beam of the specified divergence. In particular, a 0.1 mrad beam at 3 GeV is a sensitivity parameter, **not an established accelerator capability**.

The energy carried by 3 GeV neutrinos to deliver `R` raw on-symbol opportunities per second with the required expected count is represented by

`P_ν,full = R × s_required × E / p_sel`.                                       (5)

This is an idealized **full-duty-equivalent neutrino energy rate** for on symbols, in watts when `E` is in joules. With random equally probable OOK symbols, the corresponding time-average is `P_ν,full/2`; real coding and protocol overhead may change the on-symbol fraction. Eq. (5) does not compute the intensity or electrical power required to *produce* and *focus* those neutrinos. For physically achievable sources, total input power will depend on source conversion, energy distribution, duty cycle and beam optics.

The numerical calculation uses `m_u = 1.66053906660 × 10⁻²⁴ g` and `1 GeV = 1.602176634 × 10⁻¹⁰ J`. It is implemented in the accompanying `model.py`, with fixed random seed `20260923`; the scenario grid and benchmark are supplied as CSV files [11].

## 4. Results

### 4.1 Communication benchmark

| Pulses pooled per raw bit | Mean selected on-bit events | Eq. (2) BER | Simulated BER, 200,000 bits |
| ---: | ---: | ---: | ---: |
| 1 | 0.81 | 0.22243 | 0.22318 |
| 5 | 4.05 | 0.00871 | 0.00893 |
| 9 | 7.29 | 0.000341 | 0.000355 |

The simulation agrees with the analytical Bernoulli zero-count prediction at the shown precision. Five pulses give approximately 99.1% correct *uncoded synchronized bits* under these assumptions. This is consistent in scale with Stancil *et al.*'s statement that pooling five frames allowed 99% of transmitted bits to be decoded correctly, but frames and pulses are not interchangeable and this table is **not** an independent reproduction of their decoded 0.1 bit s⁻¹ performance [1].

### 4.2 Illustrative long-baseline sensitivity

Table 2 gives the full-duty-equivalent neutrino energy rate calculated from Eqs. (3)–(5), for `R = 1` raw symbol s⁻¹ and `θ = 1 mrad`. All scenarios target 1% uncoded bit error with zero background. Units are MW of neutrino-carried energy **during on-symbol-equivalent operation**, not MW of accelerator electricity.

| Chord length | 1 kt receiver | 10 kt receiver | 40 kt receiver |
| ---: | ---: | ---: | ---: |
| 1,000 km | 9.76 MW | 0.976 MW | 0.244 MW |
| 5,000 km | 244 MW | 24.4 MW | 6.10 MW |
| 12,000 km | 1,405 MW | 140.5 MW | 35.1 MW |

At fixed target mass and divergence, Eq. (4) implies `P_ν,full ∝ L² θ²/M`. The 10 kt, 5,000 km reference case is 24.4 MW at 1 mrad; setting the hypothetical divergence to 0.1 mrad changes this illustrative value to 0.244 MW, while 10 mrad changes it to 2,440 MW. Across all 27 calculated combinations, the smallest and largest entries are 0.00244 MW and approximately 140,549 MW, respectively. These extremes primarily reveal the leverage and danger of treating beam divergence as an unconstrained parameter. The complete grid is in [11].

The displayed values are conditional, optimistic estimates. They are **not** a claim that the desired rates are feasible, nor a universal lower bound on all possible neutrino communication systems. A beam energy spectrum, realistic source geometry or another detection channel changes the result; backgrounds, pointing losses, flavor oscillations and synchronization generally increase the requirements of this particular link design.

## 5. Discussion

### 5.1 Relation to earlier work

The Fermilab demonstration establishes message transfer and supplies a measured count distribution against which a simple Poisson decoder can be checked [1]. Huber's submarine analysis presents a specific proposed source and receiver concept [2]. The present long-baseline sweep does not supersede either result: it supplies a common, explicitly simplified accounting of how detected counts constrain communication and how one assumed geometry scales.

The relationship to [5] is methodological. Both studies require the analyst to separate a statistical detection score from a concrete operating point and to make detector geometry explicit. Passive submarine detection in [5] and active communication here have different energies, source control, backgrounds and receiver objectives. Direct citation to [5] should introduce this change of question and any legitimately reused methods, rather than imply that its antineutrino sensitivity numbers validate Eq. (5). Hallsjö's detector thesis [6] likewise motivates scrutiny of selection efficiency and interaction reconstruction but is not evidence that a 50% selected efficiency applies to the receiver postulated here.

### 5.2 Limits of the first calculation

**Source feasibility.** The strongest unverified assumption is the angular distribution of useful neutrinos at the selected energy. A genuine accelerator study must link pion or stored-muon production, focusing, decay geometry, proton or muon current, duty cycle, beam losses, and resulting energy–angle correlations. The generated flux must be integrated over the actual receiver footprint. Neutrino-carried energy is only one component of that calculation. The nuSTORM study [8] is a relevant source-design starting point.

**Propagation and interaction.** The calculation holds flavor survival at one. On long terrestrial chords, oscillations, including matter effects, can alter a selected muon-flavor charged-current sample [9]. A narrow 3 GeV beam also does not represent the energy spread of a real source. The appropriate rate is an energy and angle integral over source flux, oscillation probability, Earth transmission, cross section and detector efficiency; the PDG review [7] supplies cross-section context. Errors in these factors should be propagated together.

**Detector and receiver geometry.** A mass-only model hides projected area, target depth, channel choice, cosmic and atmospheric backgrounds, event containment, trigger and timing. DUNE and Hyper-Kamiokande designs [10,12] establish examples of large stationary detector technologies, but neither design report demonstrates a dedicated communication receiver. A mobile submarine receiver requires its own volume, power, timing and mechanical analysis, and should not be grouped with fixed underground detectors in an engineering claim.

**Information delivery.** Eq. (2) assumes known symbol boundaries and zero background. A useful article must select the actual symbol duration and accelerator pulse schedule, calculate nonzero-background false-positive and false-negative probabilities with an appropriate likelihood threshold, then simulate synchronization, framing, coding, dropped slots and message completion. The experimental 0.1 bit s⁻¹ decoded rate cannot be inferred from the 0.81 events per pulse alone [1]. The present `1 raw symbol s⁻¹` row is a comparison convention, not a payload-throughput prediction.

**Use-case comparison.** Any claim of practicality also requires a specified use case and competing communication path. A fixed global direct link, a submerged mobile receiver and a planetary occultation link have different requirements. The value of direct propagation cannot be inferred solely from the achievable event rate.

## 6. Work required for a submission-quality article

The next revision should (i) select one specific source design and obtain its published flux and power inputs; (ii) implement spectral transport through a stated Earth density profile and an energy-dependent interaction model; (iii) choose a detector geometry and response curve; (iv) incorporate measured or defensibly projected background in synchronized time windows; and (v) simulate complete coded messages and compute delivered information per elapsed second and per joule of facility input. The 2012 experimental setup should then be reproduced as far as its published inputs permit, explicitly recording any inaccessible raw data. Finally, the resulting source–mass–rate feasibility contours should include uncertainty ranges and assumptions sufficient for independent reproduction.

These are material missing analyses, not editorial refinements. Until completed, the numerical power values in Table 2 should be presented only as a sensitivity map.

## 7. Conclusion

A verified experimental communication link exists, and the Poisson zero-count model explains why pooling sparse events can sharply reduce an idealized raw bit-error probability. For longer through-Earth chords, a simple finite-detector model predicts a steep dependence on distance and assumed beam divergence. In the illustrative 10 kt, 1 mrad case, meeting a 1% uncoded error target at one raw symbol per second corresponds to full-duty-equivalent neutrino-carried energy rates of approximately 0.98 MW over 1,000 km and 140.5 MW over 12,000 km. The calculations isolate the scale of a technical problem; they do not establish the power, cost, error performance or feasibility of a real accelerator link. A source-specific, spectrum-integrated and protocol-level analysis is required to answer the article's motivating question.

## References

[1] D. D. Stancil *et al.*, “Demonstration of Communication using Neutrinos,” *Modern Physics Letters A* **27**, 1250077 (2012). https://doi.org/10.1142/S0217732312500770 ; https://arxiv.org/abs/1203.2847

[2] P. Huber, “Submarine neutrino communication,” *Physics Letters B* **692**, 268–271 (2010). https://doi.org/10.1016/j.physletb.2010.08.003 ; https://arxiv.org/abs/0909.4554

[3] J. G. Learned, S. Pakvasa and A. Zee, “Galactic Neutrino Communication” (2008). https://arxiv.org/abs/0805.2429

[4] J. Fidalgo Prieto *et al.*, “Submarine Navigation using Neutrinos” (2022). https://arxiv.org/abs/2207.09231

[5] S.-P. Hallsjö, “Locating nuclear-powered submarines with antineutrinos” (2026). https://arxiv.org/abs/2605.15642

[6] S.-P. Hallsjö, *Charged current quasi-elastic muon neutrino interactions in the Baby MIND detector*, PhD thesis, University of Glasgow (2018). https://theses.gla.ac.uk/41123/

[7] Particle Data Group, “Neutrino Cross Section Measurements,” *Review of Particle Physics* (2025). https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf

[8] nuSTORM Collaboration, “Neutrinos from Stored Muons (nuSTORM)” (2025). https://arxiv.org/abs/2505.06137

[9] “Neutrino Oscillations through the Earth's Core” (2021). https://arxiv.org/abs/2110.01148

[10] DUNE Collaboration, *Deep Underground Neutrino Experiment, Far Detector Technical Design Report, Volume I: Introduction to DUNE* (2020). https://arxiv.org/abs/2002.02967

[11] S.-P. Hallsjö, *Direct through-Earth neutrino communication: reproducible scoping study*, accompanying `model.py`, `benchmark.csv` and `sensitivity.csv` (working files, 2026).

[12] Hyper-Kamiokande Proto-Collaboration, *Hyper-Kamiokande Design Report* (2018). https://arxiv.org/abs/1805.04163
