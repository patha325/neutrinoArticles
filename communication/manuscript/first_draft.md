# Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity

**Sven-Patrik Hallsjö**  
**Working manuscript, version 0.6 — 23 September 2026**

> **Status.** This technical scoping study combines a packet simulation calibrated to the NuMI–MINERvA selected-event mean with a separate, deliberately simplified geometric sensitivity calculation. The simulated packets were not experimentally transmitted. The study has not established a source–detector design for a long-baseline link. Far-field figures use energy carried by neutrinos, **not proton beam power or facility electricity**; the NuMI-derived proton energy estimate describes the published target run and must not be substituted into the far-field model.

## Abstract

Neutrinos can pass through substantial matter, potentially enabling direct links where electromagnetic paths are obstructed, but their weak interactions make reception difficult. We present two separate calculations. First, a repeated OOK/CRC simulation calibrated to the published NuMI–MINERvA selected-event mean correctly accepts 65.7% of simulated 40-bit packets after five repetitions, with about 587 s mean latency under assumed external synchronization. Varying the input mean across its approximate Poisson counting interval gives 60–72% correct acceptance; this is a sensitivity range, not a joint confidence interval. Second, a distinct idealized far-field model gives 24.4 MW peak neutrino-carried on-slot power for a 10 kt receiver at 5,000 km, conditional on 3 GeV neutrinos, 1 mrad divergence, 50% selection efficiency, unit flavor survival, zero background and one-second slots. Neither result establishes a practical long-baseline source or measured CRC performance.

**Keywords:** neutrino communication; through-Earth links; Poisson channel; accelerator beam; detector sensitivity; feasibility.

![Figure 1: To-scale Earth cross-section and a 5,000 km through-Earth chord. The beam width is omitted; the diagram does not depict angular divergence.](../figures/01_earth_chord.png)

**Figure 1.** A fixed-source example used to define chord length. The source and receiver lie on the surface; a real installation would require a beamline oriented along the chord.

## 1. Introduction

Communication through rock or across a large terrestrial chord ordinarily requires an indirect route. Neutrino beams offer an unusual direct path because most neutrinos traverse matter without interacting. The same property makes reception difficult: very few transmitted particles generate identifiable events in a finite detector. Consequently, a beam that can be detected statistically need not carry a useful message at an acceptable error rate or cost.

The possibility has passed an experimental proof-of-principle test. Stancil *et al.* transmitted an encoded message using the NuMI beam and MINERvA detector over 1.035 km, including 240 m of earth. They report a decoded rate of 0.1 bit s⁻¹ and a 1% bit error rate [1]. That result establishes a link at its particular source, geometry, detector and decoding protocol; it does not answer whether a regional or global direct link could be operated at a useful rate.

Huber studied a one-way link to a submerged submarine using a high-energy neutrino beam from a muon storage ring [2]. Learned, Pakvasa and Zee examined signalling on galactic scales [3]. A later submarine-navigation proposal addresses position, navigation and timing rather than reliable delivery of an arbitrary message [4]. These studies motivate a clear separation between demonstrated transmission and projections under particular technologies.

This study also follows Hallsjö's analysis of locating nuclear-powered submarines by their emitted antineutrinos [5]. That work treated an uncontrolled, comparatively low-energy reactor source as a detection problem and examined sparse counts, background, geometry and detector deployment. Here the question is different: for a controlled and modulated source, what signal count is needed to decode bits, and how does an idealized source requirement change with baseline, angular spread and receiver size? The detector and interaction expertise documented in Hallsjö's Baby MIND thesis informs the next, detector-specific stage of this work [6]. Reactor-antineutrino event rates or detector assumptions from [5] are not transferred to the high-energy communication example.

Global neutrino telecommunications were already examined quantitatively by Sáenz *et al.* in 1977 [16]. The present novelty therefore cannot be the existence of that idea. Our aim is a transparent *scoping calculation* that connects an explicit error target to detector size and neutrino-carried energy rate, adds a finite-message example calibrated at the one measured source–receiver geometry, and tests the scaling against proposed uses with distinct receiver and latency constraints. We report what the stated assumptions imply and identify which missing inputs prevent an engineering conclusion.

## 2. Link definition and communication metric

Consider a stationary source that sends on–off-keyed symbols toward a stationary receiver through a chord of length `L`. In an ideal synchronized slot, an on symbol yields a selected signal count with Poisson mean `s`; an off symbol has no signal. A background count with mean `b` may be present in either slot. Conditional on transmitted symbol `x ∈ {0,1}`, the count model is

`K | x ~ Poisson(b + s x)`.                                                    (1)

For equally probable input symbols, `b = 0`, and a threshold of one observed event, false-positive errors vanish and an on symbol is missed with probability `exp(−s)`. Therefore

`BER_uncoded = ½ exp(−s)`.                                                     (2)

This equation measures uncoded, synchronized *raw* bits. A delivered payload rate must additionally account for symbol slots, synchronization, framing, coding, source duty cycle and any failed frames. An actual background `b > 0` introduces both false-positive and false-negative errors; neither Eq. (2) nor a threshold of one should then be assumed optimal.

We choose a target `BER_uncoded = 0.01`, so the mean signal requirement under the stated ideal conditions is

`s_required = −ln(2 × 0.01) = 3.912 selected events per on symbol`.            (3)

| Quantity | Definition used here |
| --- | --- |
| Slot duration, `τ` | Time allocated to one OOK symbol, whether 0 or 1 |
| Slot rate, `R_slot` | `1/τ` slots s⁻¹; for uncoded OOK, this equals raw input bits s⁻¹ |
| Peak on-slot power, `Pν,on` | Neutrino-carried energy required in a 1 slot, divided by `τ` |
| Full-duty average, `Pν,full` | `Pν,on` if every slot transmits a 1 |
| Time-average power, `Pν,avg` | `q₁Pν,on`, where `q₁` is the fraction of slots carrying a 1; `Pν,on/2` for equiprobable OOK |

The far-field rate model sets `τ = 1/R_slot`, so increasing the bit rate at fixed receiver and divergence requires a proportionally greater peak on-slot power. A coded scheme must instead define its coded slot rate, code rate and symbol energy explicitly.

Let `τ` be the duration of one OOK slot, `R_slot = 1/τ` the slot rate, and `q₁` the fraction of slots carrying a one. With one uncoded bit per slot, the raw bit rate is `R_raw = R_slot`. Let `Eν,1 = s Eν/p_sel` be the neutrino-carried energy needed in each on slot. The peak power while transmitting a one is `Pν,on = Eν,1/τ`. If the beam transmits in every slot, this is also its full-duty average power. For random, equiprobable uncoded OOK, the time-averaged neutrino-carried power is `Pν,avg = q₁Pν,on = Pν,on/2`. The numerical tables state **peak on-slot power** unless they explicitly say otherwise. FEC, framing and nonuniform source scheduling change the on-slot fraction and slot rate and must be handled separately.

## 3. Methods

### 3.1 Published experimental benchmark

Stancil *et al.* report an average of approximately `0.81` selected muon events for a beam-on pulse in the reduced-intensity communication run [1]. Their measurement includes charged-current interactions in upstream rock whose muons enter MINERvA, as well as a smaller in-detector component. Consequently, it is incorrect to turn that count directly into a fiducial-target interaction probability for the detector alone.

As a check of Eq. (2), we assume independent pulses and zero background, and pool `n` pulses for each bit, giving `s = 0.81 n`. A fixed-seed simulation of 200,000 equiprobable bits draws misses from `exp(−0.81n)`. Its agreement with Eq. (2) checks the implementation **only**: both use the same assumed probability, and neither independently validates the experimental data. Stancil *et al.* report roughly 78% correctly read uncoded bits in an individual synchronized frame, about 99% with five pooled frames, and no observed errors for nine or more pooled frames; they also report performance after convolutional decoding [1]. These are the published observations against which the predicted bit-error *scale* is compared, without claiming a digitized measurement or an independent dataset. Their published 0.1 bit s⁻¹ decoder-rate estimate is based on their specific frame structure and a limited number of reused frame combinations [1].

For a **simulation calibrated to one measured geometry**, use the published reduced-intensity NuMI run: `2.25 × 10¹³` 120 GeV protons per on pulse and an estimated mean `λ = 0.81` selected muons per on pulse [1]. The paper estimates this mean as `2 × 1402/3454`, assuming half the recorded slots are on. Treating the 1,402 events as Poisson gives `λ̂ = 0.8118` with an approximate 95% counting interval `[0.7693, 0.8543]`; this covers input-counting uncertainty only, not selection systematics. The on-pulse proton energy to the target is about `0.433 MJ`. This is neither facility electrical consumption nor energy delivered as neutrinos. Nor can the 0.81 observed muons be scaled linearly with an arbitrary receiver mass: upstream-rock interactions contribute most of that sample [1].

The accelerator timing is represented explicitly. The published schedule has 25 pulse slots separated by 2.2 s, followed by a 6.267 s interval; our schedule sets the 25 pulse times at offsets 0, 2.2, …, 52.8 s, and the next supercycle begins at 61.267 s [1]. Packet issue time is uniform over the supercycle. Table 1b therefore reports mean and range of modeled wall-clock latency over start phase, rather than multiplying slot count by a mean interval. It assumes the stated repeated schedule continues while a packet is sent.

Our protocol simulation carries five uniformly random payload bytes and a CRC-8/ATM check byte (polynomial 0x07, zero initial state) with `n` repeated on–off slots per bit: `48n` scheduled accelerator slots per packet. Receivers know the packet boundary and slot phase externally. Each on pulse produces independent Poisson counts with mean `λ`; the off-slot mean is zero, approximating the experiment. The decoder combines the `n` counts with a one-event threshold, calculates the CRC and accepts only a matching check. We report correctly accepted packets separately from all accepted packets and CRC-undetected errors. At `n = 3, 5, 9`, we simulate 50,000 random messages at the central `λ̂` estimate and again at both endpoints of its approximate 95% counting interval. The central result has separate finite-Monte-Carlo uncertainty; the endpoint runs provide an input-mean sensitivity interval, not a combined confidence interval. Synchronization is supplied externally, so acquisition and preamble costs are absent. This is not a reproduction of the 2012 transmitted 40-bit word, 92 coded bits, 64-bit synchronization sequence or convolutional decoder [1]. Inputs and outcomes are recorded by `empirical_message.py` and `empirical_messages.csv` [11].

### 3.2 Illustrative far-field receiver model

For a separate link sensitivity study, we assume an approximately circular far-field neutrino footprint with characteristic half-angle `θ`. At distance `L`, its area is `A = π(θL)²`. If the footprint is appreciably larger than the receiver and illumination is sufficiently uniform, an approximate selected-event probability per neutrino directed into that footprint is

`p_sel ≃ (M/m_u) σ_CC(E) ε / [π(θL)²]`,                                     (4)

where `M` is fiducial target mass, `m_u` is the atomic mass unit expressed in grams, `σ_CC` is the per-nucleon charged-current cross section and `ε` subsumes trigger and event-selection efficiency. This expression assumes a target small compared with the beam footprint and does not specify a receiver shape. It cannot be extrapolated to an arbitrarily narrow beam or a detector larger than the footprint.

We set `E = 3 GeV`, `σ_CC = 0.67 × 10⁻³⁸(E/GeV) cm²` per nucleon as an approximate illustrative parametrization, `ε = 0.5`, and vary `M = 1, 10, 40 kt`, `θ = 0.1, 1, 10 mrad`, and `L = 1,000, 5,000, 12,000 km`. The approximation gives `2.01 × 10⁻³⁸ cm²` per nucleon at 3 GeV. It is not a precision calculation: T2K reports an inclusive charged-current measurement on iron of `2.29 ± 0.45 × 10⁻³⁸ cm²` per nucleon at 3.3 GeV [28], illustrating the need for a detector-target-specific cross section and uncertainty [7]. We assume a single selected flavor at the receiver, no oscillation loss, no absorption, no off-slot background, perfect pointing and timing, and an unconstrained beam of the specified divergence. Neither 1 mrad nor 0.1 mrad at 3 GeV has been established here for a source with the modelled intensity. These are **free geometry parameters**, not a demonstrated accelerator beam.

For fixed slot duration `τ` and slot rate `R_slot = 1/τ`, the peak neutrino-carried power needed for an on slot is

`P_ν,on = R_slot × s_required × E / p_sel`.                                     (5)

This is peak power during on slots (equivalently the full-duty average if every slot is on), in watts when `E` is in joules. For equiprobable uncoded OOK, the time-average is half this value. Equation (5) does not compute the intensity or electrical power required to *produce* and *focus* those neutrinos. For physically achievable sources, total input power will depend on source conversion, energy distribution, duty cycle and beam optics.

The numerical calculation uses `m_u = 1.66053906660 × 10⁻²⁴ g` and `1 GeV = 1.602176634 × 10⁻¹⁰ J`. It is implemented in the accompanying `model.py`, with fixed random seed `20260923`; the scenario grid and benchmark are supplied as CSV files [11]. The finite-message example uses the **published NuMI source parameters at one measured receiver**, whereas Eqs. (4)–(5) posit a different, idealized monoenergetic far-field source. Their event rates and powers cannot be joined by rescaling baseline alone.

## 4. Results

### 4.1 Communication benchmark

| Pulses pooled per raw bit | Mean selected on-bit events | Eq. (2) BER | Simulated BER, 200,000 bits |
| ---: | ---: | ---: | ---: |
| 1 | 0.81 | 0.22243 | 0.22318 |
| 5 | 4.05 | 0.00871 | 0.00893 |
| 9 | 7.29 | 0.000341 | 0.000355 |

**Table 1.** Analytic prediction and implementation check for independently pooled on-pulse opportunities under ideal synchronization and zero background. This table uses Stancil et al.'s published rounded mean, `λ = 0.81`; Table 1b uses the event-count estimate `λ̂ = 2 × 1402/3454 = 0.8118`.

![Figure 2: Analytical raw bit-error probability and fixed-seed simulation versus independently pooled pulses.](../figures/02_pooled_pulses_ber.png)

**Figure 2.** Idealized zero-background, synchronized on–off keying at the published selected-event mean, 0.81 per on pulse. The seeded points are drawn from the analytical model itself and are an implementation check, not independent evidence for it.

The simulation agrees with the analytical Bernoulli zero-count prediction by construction. Five pooled on-pulse opportunities give approximately 99.1% correct *uncoded synchronized bits* in this model. Stancil *et al.* report about 99% correct bits from five pooled **frames**; since each frame supplies a corresponding bit opportunity, the count model agrees in scale, but their complete frame synchronization, observed sample and convolutional decoder are not reproduced here [1].

Table 1b instead asks whether a short *simulated* message survives repeated bits and a CRC under the **published experimental source–receiver count model**. The input mean `λ̂ = 0.8118` has approximate 95% Poisson counting limits of 0.7693–0.8543. All rows use 40 payload bits plus an eight-bit CRC, external synchronization, 50,000 random messages per run and the published pulse schedule. “Correct accepted” excludes both CRC-rejected packets and wrongly accepted payloads. At the central mean, the 95% Monte Carlo margin for the correct-acceptance fraction near 0.66 is approximately ±0.004. The low/high input-mean rows in CSV propagate counting uncertainty by simulation at the interval endpoints; they do not include systematics or constitute a joint interval.

| Repeated on slots per bit | Ideal raw bit BER | Correct packets accepted | Wrong payloads accepted | Mean latency (min–max, s) | Correct accepted payload bit/s | Mean proton beam energy to target per attempt |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 0.0438 | 11.47% | 0.090% | 351.5 (345.9–354.4) | 0.0130 | 31.1 MJ |
| 5 | 0.00863 | 65.68% | 0.002% (1 of 50,000) | 586.9 (582.2–590.7) | 0.0448 | 51.8 MJ |
| 9 | 0.000336 | 98.46% | 0 observed in 50,000 | 1,057.6 (1,054.7–1,063.2) | 0.0372 | 93.3 MJ |

**Table 1b.** Simulation calibrated to the published NuMI selected-event mean, from `empirical_messages.csv` [11]. The reported latency range is over the possible starting phases of the modeled supercycle for uniformly timed packet arrivals. “Zero observed” does not establish zero true undetected-error probability. The CRC detects errors; it does not correct them. The proton-energy column counts on-pulse energy sent to the target, averaged over messages; it excludes facility electricity. Rates exclude acquisition and return acknowledgements and cannot be extrapolated to another baseline with the far-field model.

At the lower and upper input-mean limits, simulated correct-acceptance fractions are 8.8–15.3% for three repeats, 60.0–71.9% for five, and 97.65–98.92% for nine. The corresponding accepted-correct payload rates at five repeats are 0.0409–0.0490 bit/s. These intervals show sensitivity to the estimated mean; their width also contains endpoint Monte Carlo noise, and they omit systematics.

At five repetitions, the simulated on-pulse proton energy to the target is 51.8 MJ per attempted packet, or about **2.0 MJ per correctly accepted information bit** after accounting for rejected packets. These numbers neither include power drawn by the accelerator complex nor estimate energy carried by useful neutrinos. A one-packet command has a mean modeled latency of about 587 s in this protocol, with approximately 582–591 s across start phase.

### 4.2 Illustrative long-baseline sensitivity

Table 2 gives peak on-slot neutrino-carried power calculated from Eqs. (3)–(5), for one-second slots (`R_slot = 1` slot/s) and `θ = 1 mrad`. All scenarios target 1% uncoded bit error with zero background, **unit flavor survival** and the same unverified source divergence. For equiprobable uncoded OOK, time-average beam power is half the table entry. The values should not be interpreted as achievable source powers or physical lower bounds once their fixed assumptions are relaxed.

| Chord length | 1 kt receiver | 10 kt receiver | 40 kt receiver |
| ---: | ---: | ---: | ---: |
| 1,000 km | 9.76 MW | 0.976 MW | 0.244 MW |
| 5,000 km | 244 MW | 24.4 MW | 6.10 MW |
| 12,000 km | 1,405 MW | 140.5 MW | 35.1 MW |

![Figure 3: Conditional peak on-slot neutrino-carried power versus chord length for three detector masses.](../figures/03_energy_vs_chord.png)

**Figure 3.** Conditional peak power while transmitting an on slot; the illustrative slot rate is 1 slot/s (one-second slots). For random equiprobable OOK, average neutrino-carried power is half the ordinate. It is not accelerator input power.

![Figure 4: Neutrino-carried energy-rate sensitivity to receiver mass and assumed divergence at 5,000 km.](../figures/04_mass_divergence_grid.png)

**Figure 4.** This matrix visualizes nine of the 27 computed peak on-slot power cases at one slot/s. The 0.1 mrad column is a mathematical sensitivity test; no source capable of that angular spread at 3 GeV has been demonstrated here.

At fixed target mass and divergence, Eq. (4) gives `P_ν,on ∝ L² θ²/(Mτ)`. The 10 kt, 5,000 km, one-second-slot reference case is 24.4 MW peak power during an on slot at 1 mrad; the corresponding time average for equiprobable OOK is 12.2 MW. Setting the hypothetical divergence to 0.1 mrad changes the peak reference value to 0.244 MW, while 10 mrad changes it to 2,440 MW. Across all 27 calculated combinations, the smallest and largest entries are 0.00244 MW and approximately 140,549 MW, respectively. These extremes primarily reveal the leverage and danger of treating beam divergence as an unconstrained parameter. The complete grid is in [11].

![Figure 5: Receiver mass required versus neutrino-carried on-slot power for three chord lengths.](../figures/05_detector_power_tradeoff.png)

**Figure 5.** Inversion of the same idealized model for one-second slots (one slot/s). The ordinate is peak power during an on slot. Each curve is a conditional mass–power trade-off, not a feasible accelerator design.

The model yields a compact operating-point relation. At 1 mrad, 3 GeV and 50% efficiency, for a one-second on symbol and 1% ideal raw BER,

`P_ν,on (MW) × M (kt) = 9.76 × (L/1,000 km)² MW·kt`.                        (6)

More generally, the right-hand side multiplies by `(θ/1 mrad)² × (0.5/ε) × (1 s/τ)` for slot duration `τ`. Equation (6) assumes that the beam footprint is wider than the detector and the other conditions of Section 3.2 hold. It states a **trade-off within one model**, not a minimum cost of any possible communication system.

| Chord | Product for 1 s on symbol | Receiver at 1 MW carried by neutrinos | Receiver at 10 MW | Receiver at 100 MW |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 km | 9.76 MW·kt | 9.76 kt | 0.976 kt | 0.0976 kt |
| 5,000 km | 244 MW·kt | 244 kt | 24.4 kt | 2.44 kt |
| 12,000 km | 1,405 MW·kt | 1,405 kt | 140.5 kt | 14.05 kt |

These values are generated by `tradeoffs.py` and included in `power_detector_tradeoff.csv` [11]. They expose the receiver challenge: a 10 kt detector at 5,000 km requires 24.4 MW carried by the assumed beam during a one-second on symbol; 100 MW carried by neutrinos reduces the model receiver mass to 2.44 kt, still far from a compact installation. Reducing the bit interval by a factor of 1,000 raises the required neutrino-carried power by the same factor at fixed mass. It cannot be compensated by faster neutrino propagation.

**Neutrino energy is distinct from source power.** In this deliberately linear cross-section approximation `σ_CC(E) ∝ E`, the ratio `E/σ_CC(E)` in Eq. (5) cancels when `θ` and `ε` are artificially held fixed. Raising the individual neutrino energy from 3 GeV does **not**, by itself, improve this simplified sensitivity calculation. A physical source changes its energy–angle spectrum and attainable divergence with accelerator design and neutrino energy; interaction cross sections, detector response, Earth propagation and secondary-particle ranges also change. The proton or muon beam energy and wall-plug power are separate quantities. If `η_useful` denotes the total efficiency from facility input to energy carried by neutrinos in the useful beam, the facility input would be `P_input = P_ν/η_useful`, before any other site loads. We have not established `η_useful` for a specified design and therefore cannot give a credible accelerator electrical-power number. For scale only, the published communication run's **all-on proton-pulse energy averaged over its supercycle** is approximately 0.177 MW, whereas this *different source model* requires 24.4 MW carried by useful neutrinos for the 5,000 km, 10 kt, one-slot/s reference. This is an energy-budget comparison across distinct setups, not a way to extrapolate NuMI's angular flux or selected event rate to 5,000 km.

The displayed values are conditional, optimistic estimates. They are **not** a claim that the desired rates are feasible, nor a universal lower bound on all possible neutrino communication systems. A beam energy spectrum, realistic source geometry or another detection channel changes the result; backgrounds, pointing losses, flavor oscillations and synchronization generally increase the requirements of this particular link design.

### 4.3 Raw, coded and compressed throughput

Equation (7) can be inverted into a deliberately normalized *raw slot rate*, holding 3 GeV, 50% efficiency, 1 mrad and 1% ideal uncoded bit error fixed:

`R_raw ≃ 0.1025 [P_ν,on/MW] [M/kt] (1,000 km/L)² slot s⁻¹`.             (7)

For another divergence or efficiency, multiply by `(1 mrad/θ)²(ε/0.5)`. In ordinary mass units, `M/kt = (M/kg)/10⁶`. This normalization is meaningful only for the model's **neutrino-carried peak on-slot power**, receiver mass and raw Poisson decision. `R_raw` counts every OOK slot, each carrying one uncoded bit, including both zero and one slots; it is not a rate per accelerator MW. For equiprobable OOK, time-average power is about `Pν,on/2`, while peak on-slot power is unchanged. Coding can change slot rate and on-slot fraction.

To expose the difference between physical rate and useful rate, define a *bookkeeping example* with forward-error-correction (FEC) code rate `r_c = 1/2` and a frame payload fraction `f = 0.8`. Then `R_payload,illustrative = R_raw r_c f = 0.4 R_raw`, before acquisition, retransmission, acknowledgements or decoder failures. No specific code has been shown to achieve a target post-decoder error rate under this neutrino count distribution; the 1% target refers only to the **uncoded** count decision. Coding can also allow operation at a higher raw error rate and reduce required events per symbol, so a fully optimized coded link cannot be estimated by multiplying one uncoded operating point by 0.4.

For a 5,000 km chord and 10 kt target under these assumptions, the power column is peak on-slot power. A 100 MW entry corresponds to a 4.10 slot/s model rate, or `τ ≈ 0.244 s` slots; at equiprobable uncoded OOK its time-average power would be about 50 MW.

| Peak on-slot neutrino power | Raw slots/s at 1% uncoded BER | Bookkeeping payload bits/s (assumed factors; no decoder) | If a long telemetry stream compresses 4:1: uncompressed-source-equivalent bits/s |
| ---: | ---: | ---: | ---: |
| 1 MW | 0.0410 | 0.0164 | 0.0656 |
| 10 MW | 0.410 | 0.164 | 0.656 |
| 100 MW | **4.10** | **1.64** | **6.56** |

The rightmost column is **not** a faster physical link. The illustrative bookkeeping payload remains 1.64 bit s⁻¹ in the 100 MW row; a hypothetical 4:1 lossless compressor would merely represent a source stream that previously used 6.56 bit s⁻¹. The compression factor is an explicit optimistic assumption for repetitive telemetry, not a measured result and not a plausible default for a tiny already encoded market order. Framing fractions also deteriorate for short messages. The 100 MW example assumes a large fixed detector and an unestablished 100 MW of useful neutrino-carried on-slot power.

A compact figure of merit for the **assumed coded example** is

`Q = R_payload L²/(P_ν,on M) ≃ 0.0410 bit·km²/(s·MW·kg)`,                  (8)

with `L` in km and `M` in kg. Equivalently, at 5,000 km the model yields approximately `0.00164 payload bit s⁻¹/(MW·kt)`, or `1.64×10⁻⁹ payload bit s⁻¹/(MW·kg)`. These normalized numbers must always be quoted with the energy, divergence, selection efficiency, raw BER, assumed code rate and framing fraction. `throughput.py` generates the 27-case `illustrative_throughput.csv` [11].

![Figure 7: Nominal bookkeeping payload rate versus chord for three assumed peak on-slot power levels.](../figures/07_illustrative_payload_rate.png)

**Figure 7.** Nominal bookkeeping illustration with 10 kt mass, half-rate FEC and 80% frame payload fraction. The curves assume the same 1% uncoded event threshold and do not measure performance of an implemented decoder. Power is peak during on slots; the plotted rate counts all raw OOK slots.

## 5. Discussion

### 5.1 Relation to earlier work

The Fermilab demonstration establishes message transfer and supplies a measured count distribution against which a simple Poisson decoder can be checked [1]. Huber's submarine analysis presents a specific proposed source and receiver concept [2]. The present long-baseline sweep does not supersede either result: it supplies a common, explicitly simplified accounting of how detected counts constrain communication and how one assumed geometry scales. The finite-packet exercise uses Stancil *et al.*'s measured source–receiver combination, but its CRC and repetition protocol is our construction; it does not predict event rates for other paths.

The relationship to [5] is methodological. Both studies require the analyst to separate a statistical detection score from a concrete operating point and to make detector geometry explicit. Passive submarine detection in [5] and active communication here have different energies, source control, backgrounds and receiver objectives. Direct citation to [5] should introduce this change of question and any legitimately reused methods, rather than imply that its antineutrino sensitivity numbers validate Eq. (5). Hallsjö's detector thesis [6] likewise motivates scrutiny of selection efficiency and interaction reconstruction but is not evidence that a 50% selected efficiency applies to the receiver postulated here.

### 5.2 Low-rate communication strategies from deep-space research

A useful analogue comes from NASA/JPL's response to Galileo's failed high-gain antenna. Statman describes a revised low-gain link combining compression, arraying of ground antennas, convolutional and variable-redundancy Reed–Solomon coding, decoding feedback and reprocessing of recorded data [13]. These are system-level techniques for recovering useful information from a weak link. Their numerical radio-link gains cannot be transferred to a neutrino beam, but they sharpen the question this paper should ask: how much *verified payload* arrives per unit time and source energy after acquisition and error control?

Moision and Hamkins analyze an optical photon-counting deep-space channel, jointly selecting pulse-position modulation (PPM) order and error-control code rate under average and peak power constraints [14]. Poisson event counting is the useful mathematical connection to Eq. (1). In the neutrino case, pulse timing, allowed accelerator patterns, beam-on energy, background and receiver dead time determine whether PPM or another sparse-pulse code can outperform simple on–off keying. No such gain is assumed in Tables 1–2. NASA/JPL's work on joint decoder and frame synchronization at extremely low data rates also supports treating symbol acquisition as a measured part of the link rather than granting it free of charge [15].

Stancil *et al.* give the zero-background OOK Poisson-channel capacity `C/R_slot = log₂[1 + (1 − e^(−λ)) exp(−λ/(e^λ − 1))]`, approximately **0.37 bits per pulse** at `λ ≈ 0.81` [1]. Our five-repeat CRC simulation yields about `40 × 0.6568/(48 × 5) ≈ 0.109 correctly accepted payload bits per scheduled slot at the central input mean. These are not directly equivalent performance measurements: 0.109 is a finite-packet payload yield per scheduled slot for one assumed protocol, while 0.37 is an asymptotic channel-capacity reference under Stancil et al.'s stated channel model. Neither value establishes an optimized code or the delivered rate of a practical link.

For the next analysis, compare OOK, repetition with soft count combining, and one constrained PPM scheme using the **same** average source-energy budget and an explicit maximum pulse rate. For each scheme report detection, false acquisition, coding and framing overhead; packet success probability; and delivered information bits per second and per joule. Compression should be evaluated only on a specified source message distribution; random or already compressed payloads do not offer a free gain. Combining counts from multiple receivers is a possible analogue of antenna arraying only after their acceptance, timing and independence are modelled. Delay-tolerant networking may help carry an intermittent link's data end to end, but it does not increase the underlying detected-event rate and is outside this one-hop calculation.

### 5.3 Coding and source compression choices

Galileo's low-rate recovery used compression and concatenated error-control codes [13]. For the present link, compare repetition and soft count combining against a finite-blocklength convolutional or LDPC code and a constrained sparse-pulse scheme, with synchronization and energy budgets matched. The CCSDS telemetry synchronization and channel coding standard [21] supplies established design families, not measured coding gains for neutrino counts. Decoder latency and finite frame length are critical in financial messages; a powerful long-block code may deliver fewer bits before a deadline than a short simpler code.

Lossless source coding should be chosen for the message distribution. CCSDS 121.0-B-3 specifies a lossless telemetry compression method [20]; Zstandard with a pre-shared dictionary is another documented general-purpose lossless option [22]. A fixed dictionary of permitted short commands can use fewer transmitted bits than verbose text. Compressors can expand short or nearly random messages after headers; demonstrate a 4:1 factor on an actual representative corpus before using the rightmost column of the throughput table in Section 4.3. Lossy compression may suit some imagery but cannot be applied to exact commands, orders or scientific measurements without an application-specific error budget. Neither source compression nor FEC changes neutrino interaction probability: compression removes source redundancy, while FEC adds channel redundancy to improve recoverability.

### 5.4 Submerged receivers and submarine communication

The submerged-receiver application has a different objective and constraint set from the fixed 10 kt reference. Huber's earlier proposal [2] explicitly considers a one-way high-energy beam from a muon storage ring and reception through muons produced either in the submarine or in surrounding water. That can enlarge the effective interaction volume beyond the on-board apparatus. Our mass-only, 3 GeV charged-current target model does **not** include that mechanism and cannot validate or exclude Huber's proposed rates.

It does show the price of simply miniaturizing the receiver while keeping our fixed-beam assumptions. At a 1,000 km chord, a 100 tonne fiducial target (`0.1 kt`) corresponds to about 97.6 MW carried by neutrinos during a one-second on symbol; a 10 tonne target (`0.01 kt`) corresponds to about 976 MW. At 5,000 km these become approximately 2,440 and 24,400 MW, respectively. These are extrapolations of Eq. (6), **not** numbers for Huber's water-produced-muon receiver. A practical study must model the beam's illuminated area in seawater, muon production and range, optical background, submarine depth and motion, steering and location uncertainty, plus receiver size and energy supply. A one-way downlink can deliver short commands without an onboard accelerator, but acknowledgement or a return data path would require a separate system. A submarine's location also cannot be presumed known precisely enough to keep a narrow beam aligned.

### 5.5 Low-latency financial communication: propagation versus decoding

Short messages between distant financial centres are a plausible reason to examine a direct chord. Research on low-latency trading emphasizes the competition between message reliability and time to decode [17], while measurements of Chicago–New York links document the value of optimized fiber and near-line-of-sight microwave routes [18]. Neither finding establishes a viable neutrino trading link.

Consider a **hypothetical 5,000 km Earth chord**, without naming real trading venues. Neutrino flight time is about `L/c = 16.68 ms`. The shortest surface arc between the same endpoints is about 5,138 km. Using `n_g = 1.4677` as an illustrative group index of a specified conventional optical fiber [19], ideal propagation along that surface arc is `n_g s/c = 25.15 ms`. The maximal propagation-only advantage against this idealized shortest-surface-fiber example is **8.48 ms one-way**. Actual cable routing may lengthen the fiber path, while free-space or alternative fiber technologies may offer much smaller differences. End-to-end comparisons must include source scheduling, message acquisition, decoding, routing and any acknowledgement.

![Figure 6: Ideal count-collection time versus neutrino-carried power compared with an illustrative fiber propagation advantage.](../figures/06_ideal_finance_latency.png)

**Figure 6.** The horizontal line is a propagation-only comparison for one hypothetical geometry. The sloping line is the minimum *modelled on-slot count-collection time* to reach 1% raw BER with a 10 kt receiver. Synchronization, packet overhead and processing would add delay.

In the same idealized 10 kt, 1 mrad receiver, 100 MW carried by neutrinos requires about `244 ms` to collect the target expected signal for one on symbol. Reaching `1 ms` would require about **24.4 GW** carried by neutrinos; reaching the `8.48 ms` propagation-only break-even window requires about **2.88 GW**. At 100 MW, that break-even condition instead implies an approximately **288 kt** receiver. These are conditional inversions of Eq. (6) with a one-bit, 1% raw-error criterion and unit flavor survival; they are especially fragile to the omitted beam and propagation physics. A transaction message, clock acquisition, coding, confirmation and operational redundancy add time or error constraints. The sign of the apparent speed advantage can therefore reverse once decoding is included. A one-way trading signal is also not a complete order-execution and acknowledgement path. The example identifies a research question—whether any physically achievable transmitter–receiver design could meet a specified end-to-end deadline—rather than a positive business case.

### 5.6 Limits of the first calculation

**Source feasibility.** The strongest unverified assumption is the angular distribution of useful neutrinos at the selected energy. A genuine accelerator study must link pion or stored-muon production, focusing, decay geometry, proton or muon current, duty cycle, beam losses, and resulting energy–angle correlations. The generated flux must be integrated over the actual receiver footprint. Neutrino-carried energy is only one component of that calculation. The nuSTORM study [8] is a relevant source-design starting point.

**Propagation and interaction.** The calculation holds flavor survival at one. To expose the size of this omission, a *vacuum, two-flavor illustration only* with `Δm² = 2.5 × 10⁻³ eV²` and maximal mixing yields `P(νμ→νμ) ≈ 1 − sin²[1.267 Δm²(eV²)L(km)/E(GeV)] ≈ 0.29` at 5,000 km and 3 GeV. If that illustrative factor simply multiplied selected events, the corresponding power entry would grow by about 3.5-fold. It is **not** a terrestrial prediction: three-flavor oscillations, matter profile, spectrum and flavor-sensitive efficiency change the answer [9]. A narrow 3 GeV beam also does not represent the published NuMI spectrum, which peaks near 3.2 GeV with a broad width [1]. The appropriate rate is an energy and angle integral over source flux, oscillation probability, Earth transmission, detector-target-specific cross section and efficiency, with propagated uncertainties [7,28]. At these energies absorption may be small compared with oscillation effects for Earth chords, but must still be checked for the chosen geometry and source.

**Detector and receiver geometry.** A mass-only model hides projected area, target depth, channel choice, cosmic and atmospheric backgrounds, event containment, trigger and timing. DUNE and Hyper-Kamiokande designs [10,12] establish examples of large stationary detector technologies, but neither design report demonstrates a dedicated communication receiver. A mobile submarine receiver requires its own volume, power, timing and mechanical analysis, and should not be grouped with fixed underground detectors in an engineering claim.

**Information delivery.** Eq. (2) assumes known slot boundaries and zero background. Table 1b now supplies one finite-message CRC and repetition example with the *published* NuMI pulse schedule, and records both correct accepted and undetected erroneous packets. It still assumes external synchronization, independent stationary event counts, no beam losses and no acknowledgements. For a prospective source, measure nonzero-background false-positive and false-negative probabilities and simulate synchronization, framing, code decoding, dropped slots and total message completion. The experimental 0.1 bit s⁻¹ decoded rate follows the actual published frame and convolutional code [1]; it cannot be inferred from 0.81 events per pulse alone, nor compared as a like-for-like rate with our different protocol. The far-field `1 raw slot s⁻¹` row remains a comparison convention, not a payload-throughput prediction.

**Use-case comparison.** Any claim of practicality also requires a specified use case and competing communication path. A fixed global direct link, a submerged mobile receiver and a planetary occultation link have different requirements. The value of direct propagation cannot be inferred solely from the achievable event rate.

## 6. Work required for a submission-quality article

The empirical NuMI analysis now uses its published proton intensity, time structure and selected-count mean for a *different* finite-packet protocol. The next revision must (i) choose a specific source design with a published **spectral angular flux at the receiver** and source electrical or proton-beam power; (ii) propagate that spectrum through a stated Earth density profile with three-flavor oscillations and an energy-dependent interaction model; (iii) select a receiver material, projected area and response curve, including rock-produced secondaries if relevant; (iv) incorporate measured or defensibly projected background in actual synchronized time windows; and (v) compare complete coded messages at equal average source energy, including a sparse-pulse scheme motivated by [13–15], before reporting delivered information per elapsed second and per joule of facility input. A deeper reproduction of the 2012 experiment would need its actual frame timing, source interruptions, convolutional decoder and preferably the underlying observed counts; Table 1b is **not** that reconstruction. Finally, use uncertainty ranges on flux, propagation, cross section, selection and acquisition in any source–mass–rate feasibility contour.

These are material missing analyses, not editorial refinements. Until completed, the numerical power values in Table 2 should be presented only as a sensitivity map.

## 7. Conclusion

The experimental anchor is Stancil *et al.*'s demonstrated 1.035 km link [1]. In our separate simulation calibrated to its selected-event mean, five repeated slots and CRC correctly accept 65.7% of 40-bit messages on average; modeled mean packet latency is 587 s with a 582–591 s phase range. The estimated Poisson uncertainty in the input mean shifts this acceptance to about 60–72%. These are simulated outcomes under known synchronization, not packet measurements or a recreation of the original decoder.

The independent far-field model remains a geometric sensitivity calculation. Its headline 5,000 km result—**24.4 MW peak neutrino-carried power during on slots** for a 10 kt receiver—is conditional on 3 GeV monoenergetic neutrinos, 1 mrad divergence, 50% selected efficiency, unit flavor survival, zero background, and one-second slots (1 slot/s). For equiprobable uncoded OOK, average power is half the peak. No source–detector installation meeting these assumptions has been established.

Three constraints dominate the interpretation. First, useful source intensity and focusing must be physically demonstrated and converted to facility power. Second, the detector must collect enough *selected* events against background within an actual decoding deadline, including alignment and synchronization. Third, a short decoded bit must be embedded in an end-to-end message protocol whose payload rate and reliability can be measured. Changing the neutrino's energy does not automatically evade the energy–interaction trade-off in our linear cross-section approximation; the coupled spectrum and beam optics require explicit modelling.

The separate long-stream example remains bookkeeping: 100 MW peak on-slot power, 10 kt, 5,000 km, assumed half-rate FEC and 80% frame payload produce 4.10 raw slots/s and 1.64 nominal payload bits/s before any decoder-performance or packet-success estimate. A hypothetical 4:1 compression factor changes only the represented source volume. Neither figure is measured payload throughput. Submarine, finance and Solar System cases remain scoped illustrations because each needs a distinct source, receiver and propagation model.

Huber's submarine concept uses a water-assisted muon receiver and needs a dedicated simulation [2]. In the financial illustration, an ideal 5,000 km chord has an 8.48 ms propagation-only advantage over the stated shortest-surface conventional-fiber comparison; the reference model needs about 2.88 GW peak neutrino-carried on-slot power with a 10 kt receiver just to collect one symbol within that interval, before protocol overhead. Neither application is shown to be practical. A source-specific beam, detector, propagation and protocol model is required for engineering or economic conclusions.

## Appendix A. Occulted links within the Solar System

Neutrino signalling could in principle reach a receiver behind intervening matter when an electromagnetic path is obstructed; Stancil *et al.* explicitly mentioned planetary-body blockage as a possible use [1]. Two distinct cases warrant examination. A lunar far-side installation has no direct view of Earth because the Moon lies between transmitter and receiver; NASA is developing lunar relay services for precisely that coverage gap [23]. During Mars solar conjunction, radio links are disrupted by the solar corona over an interval wider than any literal occultation by the solar disk; NASA commonly limits commanding then [24]. A direct neutrino path near or through the Sun is a separate, much less studied alternative to a radio or optical relay [25]. Neither case has an established neutrino source–receiver design; the new material here is an explicit, conditional distance-scaling calculation.

To put scale on the proposal, extend Eqs. (3)–(5) *only as a geometric sensitivity calculation*. Keep the hypothetical 3 GeV single-flavor beam, 1 mrad far-field half-angle, 10 kt receiver, 50% selected-event efficiency, zero background and 1% uncoded BER for a one-second slot. With no intervening-body absorption or flavor conversion and perfect alignment, peak neutrino-carried on-slot power is

`P_ν,on = 0.976 MW × (L/1,000 km)² × (10 kt/M) × (θ/1 mrad)² × (R_slot/1 slot s⁻¹)`.  (A1)

The corresponding selected-event target remains 3.912 per on slot. The estimates below use one-second slots (one slot/s) and report peak neutrino-carried on-slot power. For equiprobable OOK the time average is half the listed values. The values are computed by `communication/analysis/solar_system_scaling.py` using the same function and constants as the terrestrial model; they are energy flow in neutrinos, **not accelerator electrical power**. Distances use NASA's mean Earth–Moon separation and the defined astronomical unit [26,27]. A far-side surface link would be slightly longer than the centre-to-centre lunar distance, so the lunar row is a scale estimate. The AU rows are illustrative baselines, not claims about any particular Earth–Mars conjunction ephemeris.

| Illustrative separation | Peak neutrino-carried on-slot power at 1 slot/s | One-way vacuum flight time | Interpretation |
| --- | ---: | ---: | --- |
| 5,000 km | 24.4 MW | 0.0167 s | Terrestrial reference, same assumptions |
| 384,400 km | 144,000 MW (144 GW) | 1.28 s | Mean lunar-distance scale |
| 1 AU, 149.6 million km | 21.8 billion MW (21.8 PW) | 8.32 min | Interplanetary sensitivity example |
| 2 AU | 87.4 billion MW (87.4 PW) | 16.6 min | Longer illustrative baseline |

At fixed slot duration, peak power scales linearly with raw slot rate. For example, retaining 100 MW peak on-slot power and a 10 kt target at lunar distance gives about `100/144,222 ≈ 0.00069` raw slots/s. **Each slot would last about 24 minutes; the detector integrates throughout that on interval.** This is not a one-second flash repeated every 24 minutes. At 1 AU, the corresponding slot duration is about 6.9 years at 100 MW peak power. These are inversions of the stipulated event-count model, **not useful operational rates**; coding, background accumulated over such long slots and acquisition further reduce delivered payload. A notional 0.1 mrad divergence would lower each tabulated power by a factor of 100, but this paper has not demonstrated such a source at 3 GeV.

The geometry also changes the physics that the terrestrial toy model omits. For a lunar far-side link, calculate the actual chord through the Moon, its density profile, neutrino transmission and flavor evolution, then test whether a receiver of the assumed mass could be deployed there. At solar conjunction, compute the time-dependent trajectory and impact parameter through the Sun, energy-dependent interactions, oscillations in solar matter and background at the receiver. A path skimming the corona and one crossing the solar interior cannot share a universal transmission factor. Exact pointing, moving endpoints, source duty cycle and the receiver's projected area matter at both scales. Neutrinos travel at essentially the same vacuum speed as an unobstructed electromagnetic signal over these baselines; occultation avoidance does **not** remove lunar or interplanetary light time.

For an engineering comparison, define a message size, deadline and acceptable failure probability. Compare a source-specific neutrino link's delivered bits/s, electrical joules per delivered bit and total receiver mass against a lunar relay or a solar-conjunction relay plus delay-tolerant scheduling [23,25]. The NASA relay study examines architectures to preserve Earth–Mars communications during solar conjunction [25]. Under the present fixed-divergence example, distance-squared dilution alone makes an interplanetary direct neutrino link implausible at useful payload rates; any claim to the contrary needs a demonstrated narrow source, plausible receiver and full propagation and decoding calculation. The possible niche is a very short, high-value signal in a constrained occultation interval, and even that remains unproven.

## References

[1] D. D. Stancil *et al.*, “Demonstration of Communication using Neutrinos,” *Modern Physics Letters A* **27**, 1250077 (2012). https://doi.org/10.1142/S0217732312500770 ; https://arxiv.org/abs/1203.2847

[2] P. Huber, “Submarine neutrino communication,” *Physics Letters B* **692**, 268–271 (2010). https://doi.org/10.1016/j.physletb.2010.08.003 ; https://arxiv.org/abs/0909.4554

[3] J. G. Learned, S. Pakvasa and A. Zee, “Galactic Neutrino Communication” (2008). https://arxiv.org/abs/0805.2429

[4] J. Fidalgo Prieto *et al.*, “Submarine Navigation using Neutrinos” (2022). https://arxiv.org/abs/2207.09231

[5] S.-P. Hallsjö, “Locating nuclear-powered submarines with antineutrinos” (2026). https://arxiv.org/abs/2605.15642

[6] S.-P. Hallsjö, *Charged current quasi-elastic muon neutrino interactions in the Baby MIND detector*, PhD thesis, University of Glasgow (2018). https://theses.gla.ac.uk/41123/

[7] Particle Data Group, “Neutrino Cross Section Measurements,” *Review of Particle Physics* (2025). https://pdg.lbl.gov/2025/reviews/rpp2025-rev-nu-cross-sections.pdf

[8] nuSTORM Collaboration, “Neutrinos from Stored Muons (nuSTORM)” (2025). https://arxiv.org/abs/2505.06137

[9] P. B. Denton and R. Pestes, “Neutrino Oscillations through the Earth's Core,” *Physical Review D* **104**, 113007 (2021). https://doi.org/10.1103/PhysRevD.104.113007 ; https://arxiv.org/abs/2110.01148

[10] DUNE Collaboration, *Deep Underground Neutrino Experiment, Far Detector Technical Design Report, Volume I: Introduction to DUNE* (2020). https://arxiv.org/abs/2002.02967

[11] S.-P. Hallsjö, *Direct through-Earth neutrino communication: reproducible scoping study*, analysis scripts, data and figures archived in the repository snapshot at commit `7dc09e25e9ad65383ff19733ab7f88fd38aa7b9d` (2026). Reproduction checkout and commands are listed in the repository README.

[12] Hyper-Kamiokande Proto-Collaboration, *Hyper-Kamiokande Design Report* (2018). https://arxiv.org/abs/1805.04163

[13] J. I. Statman, “Optimizing the Galileo Space Communication Link,” *Interplanetary Network Progress Report* **42-116**, 114–120 (1994). https://ipnpr.jpl.nasa.gov/progress_report/42-116/116k.html

[14] B. Moision and J. Hamkins, “Deep-Space Optical Communications Downlink Budget: Modulation and Coding,” *Interplanetary Network Progress Report* **42-154**, 1–28 (2003). https://ipnpr.jpl.nasa.gov/progress_report/42-154/154K.html

[15] J. I. Statman, K.-M. Cheung, T. H. Chauvin, J. Rabkin and M. L. Belongie, “Decoder synchronization for deep space missions” (1994), NASA Technical Reports Server ID 19940025161. https://ntrs.nasa.gov/citations/19940025161

[16] A. W. Sáenz, H. Uberall, F. J. Kelly, D. W. Padgett and N. Seeman, “Telecommunication with neutrino beams,” *Science* **198**, 295–297 (1977). https://doi.org/10.1126/science.198.4314.295

[17] M. Karzand and L. R. Varshney, “Communication Strategies for Low-Latency Trading” (2015). https://arxiv.org/abs/1504.07227

[18] G. Laughlin, A. Aguirre and J. Grundfest, “Information Transmission Between Financial Markets in Chicago and New York” (2013). https://arxiv.org/abs/1302.5966

[19] Corning, *SMF-28e+ Photonic Optical Fiber* product information, effective group index 1.4677 at 1550 nm. https://www.corning.com/media/worldwide/csm/documents/Corning%20SMF28e%2B%C2%AE%20Photonic%20Specialty%20Fiber.pdf

[20] Consultative Committee for Space Data Systems, *Lossless Data Compression*, CCSDS 121.0-B-3 (2020). https://ccsds.org/Pubs/121x0b3.pdf

[21] Consultative Committee for Space Data Systems, *TM Synchronization and Channel Coding*, CCSDS 131.0-B-6 (2026). https://ccsds.org/view/bluebooks/entry/4803/

[22] Y. Collet and M. Kucherawy, “Zstandard Compression and the 'application/zstd' Media Type,” RFC 8878 (2021). https://datatracker.ietf.org/doc/html/rfc8878

[23] NASA, “Exploration and Space Communications: Lunar Communications Relay and Navigation Systems.” https://www.nasa.gov/goddard/esc/lcrns/

[24] NASA/JPL, “What's Mars Solar Conjunction, and Why Does It Matter?” (2019). https://www.jpl.nasa.gov/news/whats-mars-solar-conjunction-and-why-does-it-matter/

[25] NASA Technical Reports Server, “Lagrange-Based Options for Relay Satellites to Eliminate Earth-Mars Communications Outages During Solar Superior Conjunctions” (2020). https://ntrs.nasa.gov/citations/20205007788

