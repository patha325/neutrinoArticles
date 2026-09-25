"""Generate Paper II benchmark tables and proposal figures."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
from scipy.linalg import expm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from neutrino_channel import (
    DuneBenchmark,
    EARTH_RADIUS_KM,
    KM_TO_EV_INV,
    MATTER_POTENTIAL_EV_PER_RHO_YE,
    OSCILLATION,
    SECONDS_PER_YEAR,
    binary_ook_capacity,
    beam_pot_per_second,
    ook_error_probabilities,
    oscillation_probabilities,
    parse_reco_efficiency,
    ppm_symbol_error,
    proton_beam_energy_j_per_pot,
    prem_density_g_cm3,
    prem_path_segments,
    pmns_matrix,
    read_flux,
    read_xsec,
    repetition_ber,
    selected_spectrum_per_pot,
    simulate_crc_packets,
    slot_duration_for_signal_mean,
)

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
DATA_OUT = PACKAGE / "data" / "simulation"
FIGURES = PACKAGE / "figures"
for directory in (DATA_OUT, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
})


def save_figure(fig, name: str) -> None:
    fig.savefig(FIGURES / f"{name}.svg")
    fig.savefig(FIGURES / f"{name}.png", dpi=180)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    # Preserve the union of columns because the coding sweep includes OOK,
    # repetition, and PPM records with scheme-specific fields.
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def wilson_interval(successes: float, trials: float, z: float = 1.96) -> tuple[float, float]:
    p = successes / trials
    den = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / den
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / den
    low = 0.0 if successes == 0 else max(0.0, center - half)
    high = 1.0 if successes == trials else min(1.0, center + half)
    return low, high


def expm_probability_crosscheck(energies_gev: np.ndarray,
                               baseline_km: float,
                               density_model: str,
                               constant_density_g_cm3: float = 2.848,
                               ye: float = 0.50) -> np.ndarray:
    """Independent SciPy matrix-exponential evolution for validation only.

    The production path diagonalizes each Hermitian Hamiltonian with NumPy.
    This check uses scipy.linalg.expm for each segment and multiplies the
    resulting propagators, so the numerical evolution is implemented by an
    independent algorithm.
    """
    U = pmns_matrix(**{
        "theta12": OSCILLATION["theta12"],
        "theta13": OSCILLATION["theta13"],
        "theta23": OSCILLATION["theta23"],
        "delta_cp": OSCILLATION["delta_cp"],
    })
    mass_sq = np.diag([0.0, OSCILLATION["dm21_ev2"], OSCILLATION["dm31_ev2"]])
    vacuum_mix = U @ mass_sq @ U.conj().T
    if density_model == "prem":
        lengths, densities = prem_path_segments(baseline_km, 5.0)
    elif density_model == "constant":
        lengths = np.array([baseline_km])
        densities = np.array([constant_density_g_cm3])
    else:
        raise ValueError("density_model must be prem or constant")
    output = []
    for energy in np.asarray(energies_gev, dtype=float):
        evolution = np.eye(3, dtype=complex)
        for length, density in zip(lengths, densities):
            hamiltonian = vacuum_mix / (2.0 * energy * 1.0e9)
            hamiltonian = hamiltonian.copy()
            hamiltonian[0, 0] += MATTER_POTENTIAL_EV_PER_RHO_YE * density * ye
            evolution = expm(-1j * hamiltonian * length * KM_TO_EV_INV) @ evolution
        output.append(np.abs(evolution) ** 2)
    return np.asarray(output)


def figure_source_fold(result: dict, constant: dict) -> None:
    energy = result["energy_gev"]
    flux = result["flux_numu_per_gev_m2_pot"]
    contribution = result["differential_contribution_per_gev_per_pot"]
    p_prem = result["p_mumu"]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
    axes[0].plot(energy, flux, color="#4169a1", lw=1.5, label=r"Unoscillated $\nu_\mu$ flux")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"Flux [m$^{-2}$ GeV$^{-1}$ POT$^{-1}$]")
    axes[0].set_title("DUNE TDR FHC far-detector spectrum and selected-event fold")
    axes[0].grid(alpha=0.25, which="both")
    axp = axes[0].twinx()
    axp.plot(energy, p_prem, color="#c24b32", lw=1.2, label=r"$P(\nu_\mu\to\nu_\mu)$, PREM")
    axp.plot(energy, constant["p_mumu"], color="#c24b32", lw=1.0, ls="--",
             label=r"$P_{\mu\mu}$, constant 2.848 g cm$^{-3}$")
    axp.set_ylim(0, 1.05)
    axp.set_ylabel("Survival probability")
    l1, n1 = axes[0].get_legend_handles_labels()
    l2, n2 = axp.get_legend_handles_labels()
    axes[0].legend(l1 + l2, n1 + n2, loc="upper right")
    axes[1].plot(energy, contribution, color="#27815b", lw=1.5)
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"Neutrino energy $E_\nu$ [GeV]")
    axes[1].set_ylabel("Selected events / POT / GeV")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].set_xlim(0, 15)
    fig.tight_layout()
    save_figure(fig, "01_dune_source_propagation_fold")


def figure_earth_profile(benchmark: DuneBenchmark, result: dict) -> None:
    length = benchmark.baseline_km
    impact = np.sqrt(EARTH_RADIUS_KM**2 - (length / 2.0)**2)
    x = np.linspace(0.0, length, 1000)
    r = np.sqrt(impact**2 + (x - length / 2.0)**2)
    rho = prem_density_g_cm3(r)
    energy = result["energy_gev"]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.5))
    axes[0].plot(x, rho, color="#6e5945", lw=1.7)
    axes[0].axhline(benchmark.density_g_cm3, color="#b43c31", ls="--", lw=1.1,
                    label="DUNE GLoBES constant-density reference")
    axes[0].set_xlabel("Distance along chord [km]")
    axes[0].set_ylabel(r"PREM density [g cm$^{-3}$]")
    axes[0].set_title(f"Earth chord ({length:.1f} km)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    const_prob = oscillation_probabilities(
        energy, length, density_model="constant",
        constant_density_g_cm3=benchmark.density_g_cm3
    )
    axes[1].plot(energy, result["p_mumu"], label="PREM, layered", lw=1.5)
    axes[1].plot(energy, const_prob[:, 1, 1], label="Constant 2.848 g cm$^{-3}$",
                 ls="--", lw=1.2)
    axes[1].set_xlim(0.5, 8)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel(r"Neutrino energy $E_\nu$ [GeV]")
    axes[1].set_ylabel(r"$P(\nu_\mu\rightarrow\nu_\mu)$")
    axes[1].set_title("Three-flavor matter evolution")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "02_prem_chord_and_oscillation")


def figure_ber_sweep(rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    styles = {
        ("OOK", 0.0): ("#bd3c32", "-"),
        ("repeat-5", 0.0): ("#315c9b", "-"),
        ("OOK", 0.01): ("#bd3c32", "--"),
        ("repeat-5", 0.01): ("#315c9b", "--"),
    }
    for scheme in ("OOK", "repeat-5"):
        for b in (0.0, 0.01):
            selected = [r for r in rows if r["scheme"] == scheme and r["background_mean"] == b]
            color, linestyle = styles[(scheme, b)]
            ax.plot([r["signal_mean"] for r in selected],
                    [r["ber"] for r in selected],
                    color=color, ls=linestyle, lw=1.7,
                    label=f"{scheme}, b={b:g}/slot")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Expected signal events in an ON slot")
    ax.set_ylabel("Bit-error probability")
    ax.set_title("Few-event OOK channel: ideal MAP detection and soft repetition")
    ax.set_ylim(1e-5, 0.6)
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.tight_layout()
    save_figure(fig, "03_few_event_ber_and_repetition")


def figure_packet_success(packet_rows: list[dict]) -> None:
    fig, (ax, zoom) = plt.subplots(1, 2, figsize=(10.2, 4.8),
                                   gridspec_kw={"width_ratios": [1.55, 1.0]})
    colors = {(1, 0.0): "#bd3c32", (5, 0.0): "#315c9b",
              (1, 0.01): "#e3932d", (5, 0.01): "#2c897a"}
    for repeats in (1, 5):
        for b in (0.0, 0.01):
            chosen = [r for r in packet_rows
                      if r["repetitions"] == repeats and r["background_mean"] == b]
            chosen.sort(key=lambda r: r["signal_mean"])
            x = np.array([r["signal_mean"] for r in chosen])
            y = np.array([r["packet_success"] for r in chosen])
            lo = np.array([r["ci_low"] for r in chosen])
            hi = np.array([r["ci_high"] for r in chosen])
            y_display = np.where(y == 0.0, 0.5 / chosen[0]["packets"], y)
            label = f"repetition={repeats}, b={b:g}/slot"
            errors = [np.maximum(0.0, y_display - lo),
                      np.maximum(0.0, hi - y_display)]
            ax.errorbar(x, y_display, yerr=errors, marker="o", capsize=2.5,
                        color=colors[(repeats, b)], label=label, lw=1.2, ms=4)
            low = y_display <= 0.02
            zoom.errorbar(x[low], y_display[low],
                          yerr=[errors[0][low], errors[1][low]],
                          marker="o", capsize=2.5,
                          color=colors[(repeats, b)], label=label, lw=1.2, ms=4)
    ax.set_xscale("log")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Expected signal events in an ON slot")
    ax.set_ylabel("Correctly accepted 40-bit payload fraction")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7)
    ax.set_title("Full range")
    zoom.set_xscale("log")
    zoom.set_yscale("log")
    zoom.set_ylim(1e-5, 0.03)
    zoom.set_xlabel("Signal events per ON slot")
    zoom.set_ylabel("Packet success fraction (log scale)")
    zoom.set_title("Low-success region; zeros plotted at 0.5/N")
    zoom.grid(alpha=0.25, which="both")
    zoom.annotate("1/12,000 = 8.3e-5\n95% Wilson CI: 1.5e-5–4.7e-4",
                  xy=(1.0, 1.0/12000), xytext=(0.15, 0.004),
                  arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=7)
    fig.tight_layout()
    save_figure(fig, "04_finite_packet_success")


def figure_ppm_tradeoff(rows: list[dict]) -> None:
    """Plot error and raw rate against average signal events per slot.

    Equiprobable OOK has p_on=1/2. M-PPM sends one ON slot in M. The x-axis
    therefore uses the average signal event budget per channel slot, not the
    conditional mean in an ON slot. Error metrics are explicitly distinguished.
    """
    fig, (err, rate) = plt.subplots(1, 2, figsize=(10.0, 4.5),
                                   gridspec_kw={"width_ratios": [1.45, 1.0]})
    colors = {"OOK": "#bd3c32", "repeat-5": "#315c9b",
              "PPM-4": "#2b8c68", "PPM-8": "#8b63a9"}
    for scheme in ("OOK", "repeat-5", "PPM-4", "PPM-8"):
        for background, linestyle in ((0.0, "-"), (0.01, "--")):
            selected = sorted(
                [r for r in rows if r["scheme"] == scheme and
                 r["background_mean"] == background],
                key=lambda r: r["average_signal_events_per_slot"]
            )
            x = np.array([r["average_signal_events_per_slot"] for r in selected])
            y = np.array([r["ber"] if scheme in ("OOK", "repeat-5")
                          else r["ppm_symbol_error"] for r in selected])
            err.plot(x, y, color=colors[scheme], ls=linestyle, lw=1.5,
                     label=f"{scheme}, b={background:g}/slot")
    err.set_xscale("log")
    err.set_yscale("log")
    err.set_xlabel("Average signal events per channel slot")
    err.set_ylabel("Error probability (OOK BER; PPM symbol error)")
    err.set_title("Error at a matched average event budget")
    err.grid(alpha=0.25, which="both")
    err.legend(fontsize=7, ncol=2)
    rate.plot([0.006, 0.7], [1.0, 1.0], color=colors["OOK"], lw=1.6,
              label="OOK: 1 bit/slot")
    rate.plot([0.006, 0.7], [0.2, 0.2], color=colors["repeat-5"], lw=1.6,
              label="repeat-5: 0.2 bit/slot")
    rate.plot([0.006, 0.7], [np.log2(4)/4]*2, color=colors["PPM-4"], lw=1.6,
              label="PPM-4: 0.5 bit/slot")
    rate.plot([0.006, 0.7], [np.log2(8)/8]*2, color=colors["PPM-8"], lw=1.6,
              label="PPM-8: 0.375 bit/slot")
    rate.set_xscale("log")
    rate.set_xlim(0.006, 0.7)
    rate.set_ylim(0, 1.08)
    rate.set_xlabel("Average signal events per channel slot")
    rate.set_ylabel("Raw information rate [bits/slot]")
    rate.set_title("Uncoded raw rate")
    rate.grid(alpha=0.25, which="both")
    rate.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    save_figure(fig, "06_ppm_energy_rate_tradeoff")


def figure_quadrature_convergence(rows: list[dict]) -> None:
    order = np.array([r["quadrature_order"] for r in rows])
    relative = np.array([r["relative_difference_from_order_32"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(order, np.maximum(relative, 1e-16), marker="o", color="#315c9b")
    ax.set_yscale("log")
    ax.set_xlabel("Gauss–Legendre nodes per 250 MeV flux bin")
    ax.set_ylabel("Relative difference in folded event rate")
    ax.set_title("Energy-bin quadrature convergence (order 32 reference)")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    save_figure(fig, "07_quadrature_convergence")


def figure_latency(rate_rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    signal = np.array([r["signal_mean_per_on_slot"] for r in rate_rows])
    duration_h = np.array([r["slot_duration_hours"] for r in rate_rows])
    ax.plot(signal, duration_h, marker="o", color="#315c9b", lw=1.7)
    for x, y in zip(signal, duration_h):
        ax.annotate(f"{y:.2g} h", (x, y), textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Target mean selected events in an ON slot")
    ax.set_ylabel("ON-slot duration at nominal DUNE exposure [hours]")
    ax.set_title("Time required to accumulate few-event symbols in the 40 kt benchmark")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    save_figure(fig, "05_dune_slot_duration")


def main() -> None:
    benchmark = DuneBenchmark()
    flux = read_flux()
    xsec = read_xsec()
    result = selected_spectrum_per_pot(
        flux=flux, cross_sections=xsec,
        baseline_km=benchmark.baseline_km,
        fiducial_mass_kt=benchmark.fiducial_mass_kt,
        density_model="prem",
        max_step_km=benchmark.max_profile_step_km,
    )
    # Numerical convergence of the folded rate against a higher-order
    # quadrature reference. This tests bin integration separately from PREM
    # segment-size convergence.
    quadrature_rows = []
    quadrature_results = {}
    for order in (4, 6, 8, 12, 16, 24, 32):
        folded = selected_spectrum_per_pot(
            flux=flux, cross_sections=xsec,
            baseline_km=benchmark.baseline_km,
            fiducial_mass_kt=benchmark.fiducial_mass_kt,
            density_model="prem", max_step_km=benchmark.max_profile_step_km,
            quadrature_order=order,
        )
        quadrature_results[order] = folded["events_per_pot"]
    high_order_rate = quadrature_results[32]
    for order, rate in quadrature_results.items():
        quadrature_rows.append({
            "quadrature_order": order,
            "selected_events_per_pot": rate,
            "relative_difference_from_order_32": abs(rate / high_order_rate - 1.0),
            "reference_order": 32,
        })
    write_csv(DATA_OUT / "quadrature_convergence.csv", quadrature_rows)

    step_rows = []
    for step_km in (10.0, 5.0, 2.0, 1.0):
        folded = selected_spectrum_per_pot(
            flux=flux, cross_sections=xsec,
            baseline_km=benchmark.baseline_km,
            fiducial_mass_kt=benchmark.fiducial_mass_kt,
            density_model="prem", max_step_km=step_km,
            quadrature_order=12,
        )
        step_rows.append({"max_profile_step_km": step_km,
                          "selected_events_per_pot": folded["events_per_pot"]})
    step_reference = step_rows[-1]["selected_events_per_pot"]
    for row in step_rows:
        row["relative_difference_from_1km"] = abs(
            row["selected_events_per_pot"] / step_reference - 1.0
        )
    write_csv(DATA_OUT / "profile_step_convergence.csv", step_rows)

    # Independent numerical-method cross-check: compare the production
    # Hermitian-eigensolver propagation against SciPy's general matrix
    # exponential for constant PREM segments.
    crosscheck_rows = []
    check_energies = np.array([0.5, 1.0, 2.0, 2.5, 3.0, 5.0, 8.0])
    for density_model in ("constant", "prem"):
        production = oscillation_probabilities(
            check_energies, benchmark.baseline_km,
            density_model=density_model,
            constant_density_g_cm3=benchmark.density_g_cm3,
            max_step_km=benchmark.max_profile_step_km,
        )
        independent = expm_probability_crosscheck(
            check_energies, benchmark.baseline_km, density_model,
            constant_density_g_cm3=benchmark.density_g_cm3,
            ye=benchmark.ye_crust,
        )
        for i, energy in enumerate(check_energies):
            crosscheck_rows.append({
                "density_model": density_model,
                "energy_gev": energy,
                "p_mumu_production_eigh": production[i, 1, 1],
                "p_mumu_scipy_expm": independent[i, 1, 1],
                "max_abs_probability_difference": np.max(
                    np.abs(production[i] - independent[i])
                ),
            })
    write_csv(DATA_OUT / "oscillation_expm_crosscheck.csv", crosscheck_rows)
    constant = selected_spectrum_per_pot(
        flux=flux, cross_sections=xsec,
        baseline_km=benchmark.baseline_km,
        fiducial_mass_kt=benchmark.fiducial_mass_kt,
        density_model="constant",
        max_step_km=benchmark.max_profile_step_km,
    )

    spectrum_rows = []
    for i, energy in enumerate(result["energy_gev"]):
        spectrum_rows.append({
            "energy_gev": energy,
            "unoscillated_numu_flux_per_gev_m2_pot": result["flux_numu_per_gev_m2_pot"][i],
            "p_mumu_prem": result["p_mumu"][i],
            "p_mumu_constant_2p848": constant["p_mumu"][i],
            "cc_xsec_1e38_cm2": result["cc_xsec_1e38_cm2"][i],
            "selection_efficiency_proxy": result["efficiency_proxy"][i],
            "selected_events_per_gev_per_pot": result["selected_events_per_gev_pot"][i],
        })
    write_csv(DATA_OUT / "dune_fhc_selected_spectrum.csv", spectrum_rows)

    selected_per_pot = result["events_per_pot"]
    selected_per_year = selected_per_pot * benchmark.pot_per_year
    selected_per_second = selected_per_year / SECONDS_PER_YEAR
    pot_rate = beam_pot_per_second(benchmark.pot_per_year)
    rate_rows = []
    for signal_mean in (0.1, 0.3, 1.0, 3.0):
        duration = slot_duration_for_signal_mean(signal_mean, selected_per_second)
        rate_rows.append({
            "signal_mean_per_on_slot": signal_mean,
            "slot_duration_seconds": duration,
            "slot_duration_hours": duration / 3600.0,
            "pot_per_on_slot": pot_rate * duration,
            "selected_events_per_year_40kt": selected_per_year,
            "selected_events_per_second_40kt": selected_per_second,
            "status": "source-derived rate folded from DUNE TDR flux; effective continuous exposure; not a spill-level link simulation",
        })
    write_csv(DATA_OUT / "dune_rate_and_slot_scale.csv", rate_rows)

    sweep_rows = []
    signal_grid = np.geomspace(0.01, 5.0, 50)
    for background in (0.0, 0.01):
        for signal_mean in signal_grid:
            _, _, ber = ook_error_probabilities(signal_mean, background)
            _, _, ber5 = ook_error_probabilities(5 * signal_mean, 5 * background)
            cap, p_capacity = binary_ook_capacity(signal_mean, background, p_max=0.5)
            row_base = {
                "signal_mean": signal_mean,
                "background_mean": background,
                "average_signal_events_per_slot": 0.5 * signal_mean,
                "binary_ook_capacity_bits_per_slot_p_on_le_0p5": cap,
                "capacity_optimal_p_on": p_capacity,
                "ook_ber": ber,
                "repetition5_ber_per_payload_bit": ber5,
            }
            sweep_rows.append({
                **row_base, "scheme": "OOK", "ber": ber,
                "raw_bits_per_slot": 1.0,
                "average_signal_events_per_raw_bit": 0.5 * signal_mean,
                "error_metric": "bit-error probability",
            })
            sweep_rows.append({
                **row_base, "scheme": "repeat-5", "ber": ber5,
                "raw_bits_per_slot": 0.2,
                "average_signal_events_per_raw_bit": 2.5 * signal_mean,
                "error_metric": "bit-error probability",
            })
            for order in (4, 8):
                sweep_rows.append({
                    **row_base,
                    "scheme": f"PPM-{order}",
                    "ppm_symbol_error": ppm_symbol_error(signal_mean, background, order),
                    "ppm_raw_bits_per_slot": np.log2(order) / order,
                    "average_signal_events_per_slot": signal_mean / order,
                    "average_signal_events_per_raw_bit": signal_mean / np.log2(order),
                    "raw_bits_per_slot": np.log2(order) / order,
                    "error_metric": "symbol-error probability",
                })
    write_csv(DATA_OUT / "few_event_coding_sweep.csv", sweep_rows)

    equal_budget_rows = []
    for q in (0.1, 0.25, 0.5):
        for background in (0.0, 0.01):
            for scheme in ("OOK", "repeat-5", "PPM-4", "PPM-8"):
                if scheme == "OOK":
                    signal_on = 2.0 * q
                    error = ook_error_probabilities(signal_on, background)[2]
                    rate = 1.0
                    energy_bit = q
                    metric = "bit-error probability"
                elif scheme == "repeat-5":
                    signal_on = 2.0 * q
                    error = repetition_ber(signal_on, background, repetitions=5)
                    rate = 0.2
                    energy_bit = 5.0 * q
                    metric = "bit-error probability"
                else:
                    order = int(scheme.split("-")[1])
                    signal_on = order * q
                    error = ppm_symbol_error(signal_on, background, order)
                    rate = np.log2(order) / order
                    energy_bit = order * q / np.log2(order)
                    metric = "symbol-error probability"
                equal_budget_rows.append({
                    "average_signal_events_per_slot": q,
                    "background_events_per_slot_sensitivity": background,
                    "scheme": scheme,
                    "on_slot_signal_mean": signal_on,
                    "average_signal_events_per_raw_bit": energy_bit,
                    "proton_beam_energy_j_per_raw_bit":
                        energy_bit / selected_per_pot *
                        proton_beam_energy_j_per_pot(benchmark.proton_energy_gev),
                    "raw_bits_per_slot": rate,
                    "error_metric": metric,
                    "error_probability": error,
                    "input_convention": "equiprobable OOK; one uniformly placed ON slot per PPM symbol; repetition factor five",
                })
    write_csv(DATA_OUT / "equal_event_budget_comparison.csv", equal_budget_rows)

    packet_rows = []
    packet_points = (0.1, 0.3, 1.0, 3.0)
    trials = 12_000
    seed_index = 0
    for background in (0.0, 0.01):
        for signal_mean in packet_points:
            for repeats in (1, 5):
                result_packet = simulate_crc_packets(
                    signal_mean, background, repetitions=repeats,
                    packets=trials, seed=20260925 + seed_index
                )
                seed_index += 1
                p = result_packet["packet_success"]
                successes = int(round(p * trials))
                low, high = wilson_interval(successes, trials)
                slot_duration = signal_mean / selected_per_second
                packet_duration = 48 * repeats * slot_duration
                avg_pot_packet = 48 * repeats * 0.5 * pot_rate * slot_duration
                beam_j_packet = avg_pot_packet * (
                    benchmark.proton_energy_gev * 1.602176634e-10
                )
                delivered_bits = 40 * p
                packet_rows.append({
                    **result_packet,
                    "successful_packets": successes,
                    "signal_mean": signal_mean,
                    "background_mean": background,
                    "ci_low": low,
                    "ci_high": high,
                    "slot_duration_hours": slot_duration / 3600.0,
                    "mean_packet_latency_hours": packet_duration / 3600.0,
                    "delivered_payload_bits_per_second_no_retransmission":
                        delivered_bits / packet_duration,
                    "proton_beam_energy_j_per_successful_payload_bit":
                        beam_j_packet / delivered_bits if delivered_bits > 0 else float("inf"),
                    "facility_electrical_energy": "not estimated",
                })
    write_csv(DATA_OUT / "finite_crc_packet_simulation.csv", packet_rows)

    lengths, densities = prem_path_segments(
        benchmark.baseline_km, benchmark.max_profile_step_km
    )
    average_prem_density = float(np.average(densities, weights=lengths))
    summary = {
        "source": "DUNE TDR GLoBES FHC far-detector flux (G4LBNF v3r5p4, 2017 optimized engineered beam)",
        "flux_reference_distance_km": benchmark.flux_reference_km,
        "oscillation_baseline_km": benchmark.baseline_km,
        "flux_reference_vs_oscillation_baseline_difference_km":
            benchmark.flux_reference_km - benchmark.baseline_km,
        "earth_density_model": "PREM radial profile, midpoint-segmented at max 5 km",
        "mean_prem_density_g_cm3": average_prem_density,
        "constant_density_crosscheck_g_cm3": benchmark.density_g_cm3,
        "fiducial_mass_kt": benchmark.fiducial_mass_kt,
        "beam_power_mw": benchmark.beam_power_mw,
        "proton_energy_gev": benchmark.proton_energy_gev,
        "nominal_pot_per_year": benchmark.pot_per_year,
        "selected_numu_cc_events_per_pot_40kt": selected_per_pot,
        "selected_numu_cc_events_per_year_40kt": selected_per_year,
        "selected_numu_cc_events_per_second_40kt_calendar_average": selected_per_second,
        "calendar_average_pot_per_second": pot_rate,
        "selected_numu_cc_events_per_year_constant_density":
            constant["events_per_pot"] * benchmark.pot_per_year,
        "quadrature_order_32_reference_events_per_pot": high_order_rate,
        "quadrature_order_12_relative_difference_from_order_32":
            abs(result["events_per_pot"] / high_order_rate - 1.0),
        "max_abs_probability_difference_eigh_vs_scipy_expm": max(
            row["max_abs_probability_difference"] for row in crosscheck_rows
        ),
        "assumptions": [
            "The FD flux file includes G4LBNF beamline and geometric flux prediction but no oscillation probability.",
            "The provided flux remains at the G4LBNF far-detector histogram plane (documented as 1297 km downstream of Horn 1); oscillation evolution uses the separately configured 1284.9 km baseline. No inverse-square rescaling is applied because the flux is a location-specific beamline simulation, not a point-source fluence law. The coordinate mapping should be checked against the original geometry before a higher-precision prediction.",
            "Only nu_mu CC signal is counted; other flavors, NC, detector-specific cosmic backgrounds, migration smearing, and systematic correlations are omitted.",
            "DUNE post-selection efficiency is applied directly at true energy as a first-order proxy, not as the full reconstructed-energy response.",
            "Nominal annual POT is divided by calendar-year seconds; there is no explicit spill, live-time, or modulation schedule. The nominal POT/year is the published annual exposure convention and already reflects the accelerator scenario's projected annual exposure.",
            "Background means of 0 and 0.01 are fixed per-slot channel sensitivity scenarios, not measured predictions or fixed physical background rates. Since the slot duration varies by signal mean, these cases do not represent a constant background rate in time.",
            "Energy per delivered bit is proton-beam energy only; wall-plug or accelerator electrical energy is not estimated.",
        ],
    }
    (DATA_OUT / "simulation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    figure_source_fold(result, constant)
    figure_earth_profile(benchmark, result)
    ber_rows = []
    for background in (0.0, 0.01):
        for signal_mean in signal_grid:
            _, _, ber = ook_error_probabilities(signal_mean, background)
            _, _, ber5 = ook_error_probabilities(5*signal_mean, 5*background)
            ber_rows.extend([
                {"signal_mean": signal_mean, "background_mean": background,
                 "scheme": "OOK", "ber": ber},
                {"signal_mean": signal_mean, "background_mean": background,
                 "scheme": "repeat-5", "ber": ber5},
            ])
    figure_ber_sweep(ber_rows)
    figure_packet_success(packet_rows)
    figure_latency(rate_rows)
    figure_ppm_tradeoff(sweep_rows)
    figure_quadrature_convergence(quadrature_rows)

    print(f"Selected nu_mu CC events per POT in 40 kt: {selected_per_pot:.6e}")
    print(f"Selected events per DUNE exposure year: {selected_per_year:.3f}")
    print(f"Calendar-average selected event rate: {selected_per_second:.6e} s^-1")
    print(f"PREM path mean density: {average_prem_density:.4f} g/cm^3")
    print(f"Outputs written beneath {DATA_OUT} and {FIGURES}")


if __name__ == "__main__":
    main()
