# Paper II: End-to-End Simulation of a Neutrino Communication Channel Through the Earth

*Article draft — version 0.1; results are a reproducible benchmark calculation, not an engineering design.*

## Abstract

We build a first source-to-decoder model for a neutrino counting channel using public DUNE far-detector beam inputs. The calculation folds the unoscillated DUNE TDR forward-horn-current (FHC) flux through three-flavor matter evolution along the 1284.9 km GLoBES baseline, charged-current cross sections, and a published DUNE selection-efficiency vector. Twelve-point Gaussian quadrature within the 250 MeV flux bins limits aliasing from rapidly varying oscillation probability. In a 40 kt fiducial-mass benchmark at a nominal exposure of (1.1\times10^{21}) protons on target (POT) per year, the model yields (2.46\times10^{-18}) selected \(\nu_\mu\) charged-current events per POT, or approximately 2703 selected signal events per exposure year. Treating that annual exposure as continuous gives (8.56\times10^{-5}\) events s\(^{-1}\); accumulating 0.1, 0.3, 1, and 3 expected selected events therefore takes approximately 0.32, 0.97, 3.24, and 9.72 h, respectively. We then evaluate ideal Poisson on-off keying (OOK), pulse-position modulation (PPM), repetition-5 soft combining, and a finite 40-bit payload plus CRC-8 packet simulation. Under externally supplied synchronization, 12,000-packet Monte Carlo runs show that packet delivery remains poor at low event means and improves with repetition only by spending proportionally more channel uses. This benchmark does not include reconstructed-energy migration, spill timing, communication-specific backgrounds, synchronization acquisition, or accelerator uptime and electrical-power modeling. It is therefore a physics-informed channel study, not evidence of a practical long-baseline link. The resulting event distributions and transparent decoder baseline provide inputs for a distinct follow-on study of capacity and coding in the few-event regime.

**Keywords:** neutrino communication; DUNE; Earth matter effects; Poisson channel; finite-block coding; pulse-position modulation.

## 1. Introduction

Neutrinos pass through the Earth with little attenuation, but their weak interaction makes them difficult to receive. A communications claim must therefore connect the source to decoded information: a specified beam produces a spectrum and angular distribution; propagation changes flavor; interactions and detector selection create a sparse event record; and a protocol must turn that record into an accepted message. An assumed event mean is enough to study an abstract Poisson channel, but it does not tell us whether a source and receiver produce that mean.

Our first paper, *Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity*, separated a NuMI–MINERvA-calibrated packet exercise from an idealized far-field sensitivity calculation. The former used a reported selected-event mean and timing to simulate a new repeated-OOK/CRC protocol; it did not reproduce the experimental decoder. The far-field calculation exposed how strongly performance depends on assumed divergence, survival probability, detector efficiency, and background. Its neutrino-carried power was not accelerator electrical power. Paper II replaces the free long-baseline event-rate input with a public, source-derived flux fold and then passes the resulting sparse-count channel through explicit baseline decoders.

A second related work, *Locating nuclear-powered submarines with antineutrinos*, studies passive detection of an uncontrolled reactor-antineutrino source. That work motivates careful treatment of sparse counts and background, but it is a detection problem with a different source, spectrum, geometry, and receiver model. None of its event rates or detector assumptions are transferred here. The experimental NuMI–MINERvA demonstration remains a separate proof-of-principle reference: Stancil et al. reported a decoded rate of 0.1 bit s\(^{-1}\) and a 1% bit-error rate over 1.035 km, including 240 m of earth [1].

This paper makes a deliberately limited first step. We use DUNE TDR inputs to construct a reproducible long-baseline \(\nu_\mu\) charged-current signal model and compare basic sparse-event coding choices. A future information-theory paper, *The Neutrino Channel: Capacity and Coding in the Few-Event Regime*, will use the same channel distributions to investigate constrained capacity, practical codes, and finite-blocklength tradeoffs. The present paper supplies its physical channel benchmark; it does not attempt to settle those coding questions.

## 2. Benchmark and model

### 2.1 Source and detector reference

The source spectrum is the public DUNE TDR FHC far-detector flux histogram associated with G4LBNF v3r5p4 and the 2017 optimized-engineered beam configuration. Flux entries are in neutrinos GeV\(^{-1}\) m\(^{-2}\) POT\(^{-1}\), sampled in 250 MeV bins. We use the corresponding public GENIE charged-current cross-section table, GLoBES configuration, and FHC \(\nu_\mu\)-disappearance selection-efficiency vector. The benchmark uses a 40 kt fiducial mass, a 120 GeV proton beam, and (1.1\times10^{21}\) POT per exposure year as stated in the DUNE configuration. The flux file identifies its reference plane as 1297 km from Horn 1, while the GLoBES oscillation baseline is 1284.9 km; we retain and report this 12.1 km source-file/configuration difference rather than silently reconciling it.

This is a DUNE-like simulation benchmark, not a proposed DUNE communications mode. The public far-detector flux already contains beamline and geometric flux prediction. We do not introduce an independent beam-divergence parameter or rescale that flux as if it were a newly designed transmitter.

### 2.2 Earth propagation

We solve three-flavor Schrödinger evolution through matter. For neutrino energy (E), the flavor-basis Hamiltonian in each constant-density segment is

\[
H_f(E,x)=\frac{1}{2E}U\,\mathrm{diag}(0,\Delta m^2_{21},\Delta m^2_{31})U^\dagger
+\mathrm{diag}(V(x),0,0),
\qquad V(x)=\sqrt{2}G_F N_e(x).
\]

The PMNS parameters use the central values configured in the public DUNE GLoBES file (NuFIT 4.0 normal ordering, \(\delta_{CP}=0\) for this benchmark). Density follows the radial PREM profile, with electron fraction fixed to 0.50, and the 1284.9 km chord is split into midpoint segments no longer than 5 km. Layer evolution is multiplied in path order. Probability-matrix normalization and PMNS unitarity are tested numerically.

### 2.3 Selected event rate

For a detector containing (N_T) target nucleons, the selected signal expectation per POT is approximated by

\[
\mu_{\mathrm{sel}}/\mathrm{POT}=N_T\sum_i \Phi_{\nu_\mu,i}\,
\int_{E_i^-}^{E_i^+} dE\,P_{\mu\mu}(E,L)\,
\sigma^{CC}_{\nu_\mu}(E)\,\epsilon(E).
\]

The flux \(\Phi_{\nu_\mu,i}\) is taken constant within its published bin. The integral is evaluated with 12-point Gauss–Legendre quadrature in every bin, because evaluating a rapidly changing oscillation probability only at bin centers can alias the fold. The cross section is interpolated from the DUNE GENIE table. We apply the published selection-efficiency vector directly against true energy as a first-order proxy. In the released GLoBES configuration that vector belongs to a reconstructed-energy response; applying it without the full migration matrix is an approximation and likely the largest avoidable detector-model limitation in this first pass.

The signal model includes only selected \(\nu_\mu\) charged-current events. Neutral-current interactions, other flavors, cosmic activity, detector-specific accidental events, correlated systematic uncertainty, and selection migrations are not included. The 0 and 0.01 event-per-slot backgrounds below are sensitivity scenarios for decoder behavior, not DUNE background predictions.

### 2.4 Symbol and packet model

For binary OOK, the OFF and ON counts are Poisson with means \(b\) and \(b+s\), where \(s\) is the selected signal mean in an ON slot and \(b\) is the background mean per slot. For equiprobable symbols the MAP detector is a count threshold selected from the Poisson likelihood ratio. We report its false-alarm probability, miss probability, and bit-error probability. Repetition-5 uses five slots per payload bit and combines counts before applying the MAP threshold. Capacity in the data files is the maximum mutual information for binary OOK over input probability (0\le p_{on}\le0.5); it is a binary-input, memoryless-channel benchmark, not a system throughput result.

For PPM of order (M\), one of (M\) slots is ON and the largest observed count selects the symbol; ties are randomized. The reported exact symbol-error calculation assumes independent Poisson counts and known symbol timing. The finite packet simulation sends a random 40-bit payload followed by CRC-8/ATM in 48 OOK bits, with optional repetition-5. It uses 12,000 Monte Carlo packets per operating point and a fixed seed. Slot synchronization is supplied externally; acquisition, acknowledgement, retransmission, and synchronization loss are not simulated. Wilson intervals are plotted for packet success. CRC acceptance and undetected wrong payloads are counted separately, since a CRC is an error detector rather than a guarantee of correctness.

### 2.5 Exposure, latency, and energy boundary

For simple scale estimates, the annual POT is divided by calendar seconds. An expected signal mean (s) then corresponds to slot duration \(T_s=s/R_{sel}\). This continuous-exposure conversion is not a spill-level timing simulation: it omits the real bunch/spill structure, shutdown periods, and usable modulation windows. Packet latency multiplies this slot duration by the number of encoded slots. Proton beam energy is tallied separately where shown, using 120 GeV per proton and the expected POT. Accelerator wall-plug energy, target and focusing losses, and facility overhead are not estimated.

## 3. Numerical results

### 3.1 Source fold and propagation

The PREM chord has a path-average density of 3.059 g cm\(^{-3}\); the public GLoBES configuration uses a constant 2.848 g cm\(^{-3}\) reference. In the present rate fold, PREM gives (2.4572\times10^{-18}) selected \(\nu_\mu\) CC events/POT in 40 kt. The constant-density cross-check gives (2.4591\times10^{-18}\) events/POT, a 0.08% difference for this integrated observable under the chosen proxy response. Agreement of these two folds is a check on this one spectrum-integrated rate, not a general bound on matter-profile effects or detector uncertainty.

At (1.1\times10^{21}) POT per exposure year, the result is approximately 2703 selected events per year, or (8.56\times10^{-5}\) events per calendar second. The input spectrum, bin-averaged survival probability, and selected contribution by energy are shown in Fig. 1. Figure 2 shows the PREM chord and compares its oscillation calculation with a constant-density approximation.

![DUNE TDR FHC far flux, bin-averaged survival probability, and selected event contribution.](../figures/01_dune_source_propagation_fold.svg)

*Figure 1. Public DUNE TDR flux and the first-order selected \(\nu_\mu\) CC fold. The bottom panel is the bin-integrated contribution divided by the 250 MeV bin width. The true-energy efficiency proxy is not a full detector migration model.*

![PREM chord and oscillation comparison.](../figures/02_prem_chord_and_oscillation.svg)

*Figure 2. The 1284.9 km chord and three-flavor matter survival probability for the benchmark. The path-averaged PREM density differs from the constant-density GLoBES reference.*

### 3.2 Few-event symbols and finite packets

At an expected 0.1 selected signal events per ON slot, ideal equiprobable OOK with zero background has a BER of 0.452; at one event per ON slot, the BER is 0.184. A background mean of 0.01 events per slot changes the MAP threshold as the signal-to-background ratio changes. Repetition-5 reduces bit errors at a fixed per-slot mean by accumulating five times as much exposure, while reducing uncoded payload rate by the same factor before framing overhead. Figure 3 compares those ideal symbol-level results.

![Few-event OOK BER and soft repetition.](../figures/03_few_event_ber_and_repetition.svg)

*Figure 3. Exact Poisson MAP BER for OOK and count-combined repetition-5, for two assumed background means. These curves omit synchronization, framing, and retransmission.*

The finite packet experiment makes the framing cost visible. For example, at one expected signal event per ON slot and zero background, the simulated 40-bit payload success fraction is (8.3\times10^{-5}\) without repetition and 0.855 with repetition-5. With repetition-5 at three events per ON slot and zero background, all 12,000 packets decode correctly; the 95% Wilson lower bound is about 0.99968, not 1. At one event per slot, the 48-bit frame takes about 155.5 h without repetition and 777.7 h with repetition under the continuous-exposure conversion. The repetition-5 case has a much higher success probability but still a very low successful payload rate because each symbol itself takes hours. These are simulation outcomes under externally supplied timing and the selected rate model, not measured communications performance.

At a 0.01 event-per-slot background, success is lower at the same signal mean. CRC acceptance is nonzero even when no complete payload is correct; the CRC therefore cannot be interpreted as a message-delivery guarantee. Figure 4 and `finite_crc_packet_simulation.csv` report correct packet acceptance, CRC acceptance, and undetected wrong payload outcomes for every simulated operating point.

![Finite packet success with CRC-8 and externally supplied synchronization.](../figures/04_finite_packet_success.svg)

*Figure 4. Correctly accepted 40-bit payload fraction from 12,000 packet trials per point; vertical bars are 95% Wilson intervals. Synchronization is assumed available.*

The event-rate conversion puts the sparse-symbol regime on an operational timescale: 0.1 selected signal events require about 0.324 h, 0.3 events 0.972 h, one event 3.24 h, and three events 9.72 h of nominal exposure. These intervals are per ON slot; an OFF slot and the remaining frame slots add wall-clock time. Figure 5 summarizes the signal accumulation scale.

![Slot duration for selected-event means.](../figures/05_dune_slot_duration.svg)

*Figure 5. Continuous-exposure slot duration required to accumulate the indicated selected signal mean. Real spill timing and uptime are not included.*

## 4. Discussion

The main result is methodological: with public DUNE-like flux and detector inputs, even a 40 kt fiducial receiver yields only about 2700 selected \(\nu_\mu\) CC events per nominal exposure year in this source-rate fold. At the resulting calendar-average rate, a single binary slot with a sub-event mean lasts tens of minutes to hours. Framing multiplies that duration by dozens of encoded bits, and repetition improves reliability only by consuming additional slots. Sparse-event capacity curves are useful for comparing channel strategies, but they do not remove the rate and latency implied by the source and detector.

The numerical result is not a direct prediction of a functioning communication link. The source flux is a far-detector simulation for an oscillation experiment, and the detector response is simplified by using a reconstructed-energy efficiency against true energy. Most importantly for a real decoder, no spill-resolved source schedule, beam-on/beam-off background model, or clock-acquisition procedure has been propagated into the packet simulation. Realistic duty cycle and live-time can make calendar latency longer than the continuous-POT estimate. Conversely, a different source, spectrum, receiver geometry, or coded modulation would define a different benchmark and must be modeled explicitly.

The full migration matrix is the next necessary physics improvement. After that, the model should use spill-level flux and detector live-time, include measured beam-on and beam-off backgrounds, and propagate flux, cross-section, oscillation, and selection uncertainties. A synchronization and packet protocol can then be tested against event timestamps. Facility-input energy requires an independently specified accelerator operating point and cannot be inferred from neutrino-carried energy or POT alone. These improvements are prerequisites before interpreting the outputs as achievable throughput or energy per delivered bit.

Paper I and the present article remain complementary: Paper I records the NuMI-calibrated communication demonstration context and a separate idealized sensitivity exercise; Paper II derives a first source-based event rate and sends sparse counts through explicit baseline decoders. The passive submarine antineutrino study remains methodologically related only through low-count inference. The capacity-and-coding follow-on should inherit the channel law and input constraints produced here while studying the coding problem in more depth.

## 5. Conclusion

We implemented and tested a reproducible DUNE TDR-based source-to-count model and a finite-message Poisson decoder baseline. The model predicts approximately 2703 selected \(\nu_\mu\) CC signal events per 40 kt exposure-year under the stated flux, oscillation, cross-section, and efficiency-proxy assumptions. At nominal continuous exposure, accumulating even one selected event takes about 3.24 h. Packet simulations illustrate the reliability–latency tradeoff between uncoded OOK and repetition-5, but omit synchronization acquisition, spill timing, detector migration, realistic backgrounds, retransmission, and accelerator electrical energy. The code, inputs, tables, tests, and figures are published alongside this draft to make those limits visible and the next model improvements concrete.

## Data and code availability

The analysis is in `communication/paper_ii/analysis/`; public DUNE inputs and their provenance are in `communication/paper_ii/data/inputs/`; generated CSV and JSON outputs are in `communication/paper_ii/data/simulation/`; and vector figures are in `communication/paper_ii/figures/`. Reproduce with:

```bash
python -m pip install -r communication/paper_ii/analysis/requirements.txt
python communication/paper_ii/analysis/run_simulation.py
python -m unittest discover -s communication/paper_ii/analysis -p 'test_*.py'
```

## References

1. D. D. Stancil et al., “Demonstration of Communication using Neutrinos,” *Modern Physics Letters A* **27**, 1250077 (2012), [arXiv:1203.2847](https://arxiv.org/abs/1203.2847).
2. DUNE Collaboration, “Experiment Simulation Configurations Approximating DUNE TDR,” arXiv:2103.04797 (2021), [doi:10.48550/arXiv.2103.04797](https://doi.org/10.48550/arXiv.2103.04797). The public posting includes simulated flux and GLoBES ancillary inputs.
3. A. M. Dziewonski and D. L. Anderson, “Preliminary reference Earth model,” *Physics of the Earth and Planetary Interiors* **25**, 297–356 (1981), [doi:10.1016/0031-9201(81)90046-7](https://doi.org/10.1016/0031-9201(81)90046-7).
4. S.-P. Hallsjö, “Direct Neutrino Communication Through the Earth: NuMI-Calibrated Simulation and Far-Field Sensitivity,” working manuscript, Paper I, repository version.
5. S.-P. Hallsjö, “Locating nuclear-powered submarines with antineutrinos,” arXiv:2605.15642 (2026).
