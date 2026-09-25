"""DUNE-TDR-based through-Earth neutrino counting-channel model.

The source spectrum is the public DUNE TDR GLoBES far-detector flux input.
This module folds its unoscillated nu_mu spectrum through three-flavor matter
evolution, the associated GENIE CC cross-section table, and the published
FHC nu_mu-disappearance selection-efficiency vector. It is a reproducible
first-order channel model, not a detector simulation: migration smearing,
communication-specific backgrounds, accelerator spill timing, and actual
beam modulation are not represented.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import poisson


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data" / "inputs"
FLUX_FILE = INPUTS / "dune_tdr_fhc_fd_flux.txt"
XSEC_FILE = INPUTS / "dune_tdr_cc_cross_sections.dat"
EFF_FILE = INPUTS / "dune_tdr_numu_fhc_selection_efficiency.txt"
GLOBES_FILE = INPUTS / "dune_tdr_globes_config.glb"

EARTH_RADIUS_KM = 6371.0
KM_TO_EV_INV = 5.067730716e9
MATTER_POTENTIAL_EV_PER_RHO_YE = 7.63247e-14
SECONDS_PER_YEAR = 365.25 * 24 * 3600
AVOGADRO = 6.02214076e23

# Central values used in the public DUNE TDR GLoBES configuration (NuFIT 4.0,
# normal ordering). The phase is fixed to zero for the benchmark.
OSCILLATION = {
    "theta12": 0.5903,
    "theta13": 0.1500,
    "theta23": 0.8660,
    "delta_cp": 0.0,
    "dm21_ev2": 7.39e-5,
    "dm32_ev2": 2.451e-3,
}
OSCILLATION["dm31_ev2"] = OSCILLATION["dm32_ev2"] + OSCILLATION["dm21_ev2"]


@dataclass(frozen=True)
class DuneBenchmark:
    baseline_km: float = 1284.9
    flux_reference_km: float = 1297.0
    density_g_cm3: float = 2.848
    fiducial_mass_kt: float = 40.0
    proton_energy_gev: float = 120.0
    beam_power_mw: float = 1.2
    pot_per_year: float = 1.1e21
    ye_crust: float = 0.50
    max_profile_step_km: float = 5.0


def read_flux(path: Path = FLUX_FILE) -> np.ndarray:
    """Read GLoBES FD flux rows: E, nue, numu, nutau, nuebar, numubar, nutaubar.

    Flux columns have units neutrinos / GeV / m^2 / POT. The source file is a
    flux at the DUNE far-detector reference plane; it already contains the
    beamline and geometric flux simulation. It is not oscillated.
    """
    table = np.loadtxt(path)
    if table.ndim != 2 or table.shape[1] < 7:
        raise ValueError("Expected seven columns in the DUNE GLoBES flux file.")
    if not np.all(np.diff(table[:, 0]) > 0):
        raise ValueError("Flux energy bins must be strictly increasing.")
    return table[:, :7]


def read_xsec(path: Path = XSEC_FILE) -> np.ndarray:
    table = np.loadtxt(path)
    if table.ndim != 2 or table.shape[1] < 7:
        raise ValueError("Expected log10(E) and six flavor cross-section columns.")
    return table


def parse_reco_efficiency(path: Path = EFF_FILE,
                          globes_path: Path = GLOBES_FILE) -> tuple[np.ndarray, np.ndarray]:
    """Return DUNE reconstructed-bin centers and post-selection efficiencies."""
    text = globes_path.read_text(encoding="utf-8")
    match = re.search(r"\$binsize\s*=\s*\{([^}]*)\}", text, flags=re.S)
    if not match:
        raise ValueError("Could not find the GLoBES reconstructed bin widths.")
    widths = np.array([float(v) for v in match.group(1).split(",") if v.strip()])
    emin_match = re.search(r"\$emin\s*=\s*([0-9.eE+-]+)", text)
    emin = float(emin_match.group(1)) if emin_match else 0.0
    edges = np.concatenate(([emin], emin + np.cumsum(widths)))
    centers = 0.5 * (edges[:-1] + edges[1:])
    eff_text = path.read_text(encoding="utf-8")
    eff_match = re.search(r"\{([^}]*)\}", eff_text, flags=re.S)
    if not eff_match:
        raise ValueError("Could not parse the GLoBES efficiency vector.")
    efficiency = np.array([float(v) for v in eff_match.group(1).split(",") if v.strip()])
    if len(centers) != len(efficiency):
        raise ValueError(
            f"GLoBES binning ({len(centers)}) and efficiency ({len(efficiency)}) differ."
        )
    return centers, efficiency


def pmns_matrix(theta12: float, theta13: float,
                theta23: float, delta_cp: float) -> np.ndarray:
    """Standard three-flavor PMNS matrix in the PDG parameterization."""
    s12, c12 = np.sin(theta12), np.cos(theta12)
    s13, c13 = np.sin(theta13), np.cos(theta13)
    s23, c23 = np.sin(theta23), np.cos(theta23)
    ep = np.exp(1j * delta_cp)
    return np.array([
        [c12 * c13, s12 * c13, s13 / ep],
        [-s12 * c23 - c12 * s23 * s13 * ep,
         c12 * c23 - s12 * s23 * s13 * ep, s23 * c13],
        [s12 * s23 - c12 * c23 * s13 * ep,
         -c12 * s23 - s12 * c23 * s13 * ep, c23 * c13],
    ], dtype=complex)


def prem_density_g_cm3(radius_km: np.ndarray | float) -> np.ndarray:
    """PREM radial density, with r in km; polynomial coefficients use r/R."""
    r = np.asarray(radius_km, dtype=float)
    x = np.clip(r / EARTH_RADIUS_KM, 0.0, 1.0)
    rho = np.empty_like(x)
    # Dziewonski-Anderson PREM radial density polynomials and shallow layers.
    m = x < 1221.5 / EARTH_RADIUS_KM
    rho[m] = 13.0885 - 8.8381 * x[m] ** 2
    m = (x >= 1221.5 / EARTH_RADIUS_KM) & (x < 3480.0 / EARTH_RADIUS_KM)
    rho[m] = 12.5815 - 1.2638*x[m] - 3.6426*x[m]**2 - 5.5281*x[m]**3
    m = (x >= 3480.0 / EARTH_RADIUS_KM) & (x < 5701.0 / EARTH_RADIUS_KM)
    rho[m] = 7.9565 - 6.4761*x[m] + 5.5283*x[m]**2 - 3.0807*x[m]**3
    m = (x >= 5701.0 / EARTH_RADIUS_KM) & (x < 5771.0 / EARTH_RADIUS_KM)
    rho[m] = 5.3197 - 1.4836*x[m]
    m = (x >= 5771.0 / EARTH_RADIUS_KM) & (x < 5971.0 / EARTH_RADIUS_KM)
    rho[m] = 11.2494 - 8.0298*x[m]
    m = (x >= 5971.0 / EARTH_RADIUS_KM) & (x < 6151.0 / EARTH_RADIUS_KM)
    rho[m] = 7.1089 - 3.8045*x[m]
    m = (x >= 6151.0 / EARTH_RADIUS_KM) & (x < 6346.6 / EARTH_RADIUS_KM)
    rho[m] = 2.6910 + 0.6924*x[m]
    m = (x >= 6346.6 / EARTH_RADIUS_KM) & (x < 6356.0 / EARTH_RADIUS_KM)
    rho[m] = 2.9
    m = x >= 6356.0 / EARTH_RADIUS_KM
    rho[m] = 2.6
    return rho


_PREM_BOUNDARIES_KM = np.array(
    [1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6371.0]
)


def prem_path_segments(baseline_km: float, max_step_km: float = 5.0
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Return segment lengths and midpoint PREM densities along a straight chord."""
    if not 0.0 < baseline_km <= 2.0 * EARTH_RADIUS_KM:
        raise ValueError("Baseline must be positive and no longer than Earth's diameter.")
    half = baseline_km / 2.0
    impact = np.sqrt(max(EARTH_RADIUS_KM**2 - half**2, 0.0))
    edges = [0.0, baseline_km]
    for radius in _PREM_BOUNDARIES_KM[:-1]:
        if radius > impact:
            offset = np.sqrt(radius**2 - impact**2)
            edges.extend([half - offset, half + offset])
    edges = np.unique(np.clip(edges, 0.0, baseline_km))
    refined = [edges[0]]
    for lo, hi in zip(edges[:-1], edges[1:]):
        n = max(1, int(np.ceil((hi - lo) / max_step_km)))
        refined.extend(np.linspace(lo, hi, n + 1)[1:])
    refined = np.asarray(refined)
    mids = 0.5 * (refined[:-1] + refined[1:])
    radii = np.sqrt(impact**2 + (mids - half)**2)
    return np.diff(refined), prem_density_g_cm3(radii)


def oscillation_probabilities(energies_gev: np.ndarray,
                              baseline_km: float,
                              density_model: str = "prem",
                              constant_density_g_cm3: float = 2.848,
                              ye: float = 0.50,
                              max_step_km: float = 5.0,
                              parameters: dict[str, float] | None = None,
                              antineutrino: bool = False) -> np.ndarray:
    """Compute P[out_flavor, in_flavor] by piecewise-constant Hamiltonian evolution."""
    pars = OSCILLATION if parameters is None else parameters
    U = pmns_matrix(pars["theta12"], pars["theta13"], pars["theta23"], pars["delta_cp"])
    if antineutrino:
        U = U.conj()
    masses = np.diag([0.0, pars["dm21_ev2"], pars["dm31_ev2"]]).astype(complex)
    energies_ev = np.asarray(energies_gev, dtype=float) * 1.0e9
    if np.any(energies_ev <= 0.0):
        raise ValueError("Neutrino energies must be positive.")
    vacuum = (U @ masses @ U.conj().T)[None, :, :] / (2.0 * energies_ev[:, None, None])
    if density_model == "prem":
        lengths, densities = prem_path_segments(baseline_km, max_step_km)
    elif density_model == "constant":
        lengths = np.array([baseline_km])
        densities = np.array([constant_density_g_cm3])
    else:
        raise ValueError("density_model must be 'prem' or 'constant'.")
    evolution = np.broadcast_to(np.eye(3, dtype=complex),
                                (len(energies_ev), 3, 3)).copy()
    sign = -1.0 if antineutrino else 1.0
    for length, density in zip(lengths, densities):
        potential = sign * MATTER_POTENTIAL_EV_PER_RHO_YE * density * ye
        hamiltonian = vacuum.copy()
        hamiltonian[:, 0, 0] += potential
        eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
        phases = np.exp(-1j * eigenvalues * length * KM_TO_EV_INV)
        step = np.einsum("eik,ek,ejk->eij",
                         eigenvectors, phases, eigenvectors.conj(), optimize=True)
        evolution = np.einsum("eij,ejk->eik", step, evolution, optimize=True)
    return np.abs(evolution) ** 2


def cc_cross_section_1e38_cm2(energies_gev: np.ndarray,
                               table: np.ndarray | None = None,
                               flavor_column: int = 1) -> np.ndarray:
    """Interpolate GLoBES GENIE sigma/E table; return sigma in 10^-38 cm^2."""
    if table is None:
        table = read_xsec()
    energy = np.asarray(energies_gev, dtype=float)
    y = np.interp(np.log10(np.maximum(energy, 1e-6)),
                  table[:, 0], table[:, 1 + flavor_column])
    return np.maximum(y * energy, 0.0)


def selected_spectrum_per_pot(flux: np.ndarray | None = None,
                              cross_sections: np.ndarray | None = None,
                              baseline_km: float = 1284.9,
                              fiducial_mass_kt: float = 40.0,
                              density_model: str = "prem",
                              max_step_km: float = 5.0
                              ) -> dict[str, np.ndarray | float]:
    """Fold DUNE FHC nu_mu flux through P(mu->mu), CC xsec, and DUNE efficiency.

    Efficiency is indexed by the GLoBES reconstructed-energy bins and applied
    directly at true energy here; migration smearing is omitted. This is a
    transparent efficiency proxy, not the full DUNE far-detector response.
    """
    if flux is None:
        flux = read_flux()
    if cross_sections is None:
        cross_sections = read_xsec()
    energy = flux[:, 0]
    flux_numu = flux[:, 2]
    # The public flux is a histogram on 250 MeV bins. Oscillation probability
    # can vary rapidly inside those bins, so integrate each bin with Gaussian
    # quadrature instead of evaluating P at its center and aliasing the result.
    edges = np.empty(len(energy) + 1)
    edges[1:-1] = 0.5 * (energy[:-1] + energy[1:])
    edges[0] = max(0.0, energy[0] - 0.5 * (energy[1] - energy[0]))
    edges[-1] = energy[-1] + 0.5 * (energy[-1] - energy[-2])
    nodes, weights = np.polynomial.legendre.leggauss(12)
    halfwidth = 0.5 * np.diff(edges)
    quad_energy = (energy[:, None] + halfwidth[:, None] * nodes[None, :]).ravel()
    quad_weights = np.broadcast_to(halfwidth[:, None] * weights[None, :],
                                   (len(energy), len(nodes))).ravel()
    probs = oscillation_probabilities(
        quad_energy, baseline_km, density_model=density_model, max_step_km=max_step_km
    )
    p_mumu_quad = probs[:, 1, 1].reshape(len(energy), len(nodes))
    center_probs = oscillation_probabilities(
        energy, baseline_km, density_model=density_model, max_step_km=max_step_km
    )
    p_mumu_center = center_probs[:, 1, 1]
    p_mumu = np.sum(p_mumu_quad * quad_weights.reshape(len(energy), len(nodes)), axis=1) / np.diff(edges)
    reco_centers, efficiency_values = parse_reco_efficiency()
    efficiency = np.interp(energy, reco_centers, efficiency_values, left=0.0, right=0.0)
    quad_efficiency = np.interp(quad_energy, reco_centers, efficiency_values,
                                left=0.0, right=0.0).reshape(len(energy), len(nodes))
    sigma = cc_cross_section_1e38_cm2(energy, cross_sections, flavor_column=1)
    quad_sigma = cc_cross_section_1e38_cm2(
        quad_energy, cross_sections, flavor_column=1
    ).reshape(len(energy), len(nodes))
    nucleons = fiducial_mass_kt * 1.0e9 * AVOGADRO
    # flux is m^-2 GeV^-1 POT^-1; cross section table is 10^-38 cm^2.
    # Multiply by 1e-4 to convert cm^2 to m^2 and 1e-38 for table units.
    selected_per_gev_pot = (
        flux_numu * p_mumu * sigma * 1.0e-42 * nucleons * efficiency
    )
    binned_yield_per_pot = flux_numu * 1.0e-42 * nucleons * np.sum(
        p_mumu_quad * quad_sigma * quad_efficiency *
        (quad_weights.reshape(len(energy), len(nodes))), axis=1
    )
    yield_per_pot = float(np.sum(binned_yield_per_pot))
    return {
        "energy_gev": energy,
        "flux_numu_per_gev_m2_pot": flux_numu,
        "p_mumu": p_mumu,
        "p_mumu_bin_center": p_mumu_center,
        "efficiency_proxy": efficiency,
        "cc_xsec_1e38_cm2": sigma,
        "selected_events_per_gev_pot": selected_per_gev_pot,
        "selected_events_per_bin_per_pot": binned_yield_per_pot,
        "bin_width_gev": np.diff(edges),
        "differential_contribution_per_gev_per_pot": binned_yield_per_pot / np.diff(edges),
        "events_per_pot": yield_per_pot,
    }


def map_threshold(signal_mean: float, background_mean: float) -> int:
    """MAP threshold for equiprobable binary OOK symbols and Poisson counts."""
    if signal_mean < 0 or background_mean < 0:
        raise ValueError("Poisson means cannot be negative.")
    if background_mean == 0.0:
        return 1
    if signal_mean == 0.0:
        return np.iinfo(np.int32).max
    return max(0, int(np.ceil(signal_mean / np.log1p(signal_mean / background_mean))))


def ook_error_probabilities(signal_mean: float, background_mean: float = 0.0,
                            threshold: int | None = None) -> tuple[float, float, float]:
    """Return P(false alarm), P(miss), BER for equiprobable OOK."""
    if threshold is None:
        threshold = map_threshold(signal_mean, background_mean)
    if threshold <= 0:
        false_alarm = 1.0
    else:
        false_alarm = float(poisson.sf(threshold - 1, background_mean))
    miss = float(poisson.cdf(threshold - 1, background_mean + signal_mean))
    return false_alarm, miss, 0.5 * (false_alarm + miss)


def repetition_ber(signal_mean: float, background_mean: float = 0.0,
                   repetitions: int = 5) -> float:
    """MAP BER for repetition coding with soft combining of Poisson counts."""
    f, m, ber = ook_error_probabilities(
        repetitions * signal_mean, repetitions * background_mean
    )
    return ber


def binary_ook_mutual_information(signal_mean: float, background_mean: float = 0.0,
                                  p_on: float = 0.5) -> float:
    """I(X;N) in bits/slot for OOK with a fixed Bernoulli input probability."""
    if not 0.0 <= p_on <= 1.0:
        raise ValueError("p_on must lie in [0, 1].")
    if p_on in (0.0, 1.0):
        return 0.0
    maximum = int(poisson.ppf(1 - 1e-13, background_mean + signal_mean) + 2)
    counts = np.arange(maximum + 1)
    p0 = poisson.pmf(counts, background_mean)
    p1 = poisson.pmf(counts, background_mean + signal_mean)
    py = (1.0 - p_on) * p0 + p_on * p1
    result = 0.0
    for px, pxy in ((1.0 - p_on, p0), (p_on, p1)):
        mask = (pxy > 0.0) & (py > 0.0)
        result += px * float(np.sum(pxy[mask] * np.log2(pxy[mask] / py[mask])))
    return result


def binary_ook_capacity(signal_mean: float, background_mean: float = 0.0,
                        p_max: float = 1.0) -> tuple[float, float]:
    """Maximize binary-input OOK mutual information for 0 <= p_on <= p_max."""
    fit = minimize_scalar(
        lambda p: -binary_ook_mutual_information(signal_mean, background_mean, p),
        bounds=(1e-8, p_max), method="bounded",
        options={"xatol": 1e-8},
    )
    p = float(np.clip(fit.x, 0.0, p_max))
    return -float(fit.fun), p


def ppm_symbol_error(signal_mean: float, background_mean: float = 0.0,
                     order: int = 4) -> float:
    """Exact max-count PPM symbol-error probability; ties are broken uniformly."""
    if order < 2:
        raise ValueError("PPM order must be >= 2.")
    max_count = int(poisson.ppf(1 - 1e-13, background_mean + signal_mean) + 2)
    k = np.arange(max_count + 1)
    p_signal = poisson.pmf(k, background_mean + signal_mean)
    p_below = poisson.cdf(k - 1, background_mean)
    p_equal = poisson.pmf(k, background_mean)
    p_correct = np.zeros_like(k, dtype=float)
    for ties in range(order):
        p_correct += (
            math.comb(order - 1, ties)
            * p_equal**ties
            * p_below**(order - 1 - ties)
            / (ties + 1)
        )
    success = float(np.sum(p_signal * p_correct))
    return float(np.clip(1.0 - success, 0.0, 1.0))


def crc8_atm(data: bytes) -> int:
    """CRC-8/ATM, poly=0x07, init=0, refin/refout=false, xorout=0."""
    crc = 0
    for byte in data:
        crc ^= byte
        for _ in range(8):
            crc = ((crc << 1) ^ 0x07) & 0xFF if crc & 0x80 else (crc << 1) & 0xFF
    return crc


def _bits_from_bytes(raw: np.ndarray) -> np.ndarray:
    return np.unpackbits(raw, axis=1, bitorder="big")


def simulate_crc_packets(signal_mean: float, background_mean: float = 0.0,
                         repetitions: int = 1, packets: int = 20_000,
                         seed: int = 20260925) -> dict[str, float]:
    """Monte Carlo random 40-bit payload + CRC-8 packet through repeated OOK.

    Slot synchronization is supplied externally. The CRC detects most decoded
    errors, but this routine does not simulate acquisition or acknowledgements.
    """
    rng = np.random.default_rng(seed)
    payload = rng.integers(0, 256, size=(packets, 5), dtype=np.uint8)
    crc = np.fromiter((crc8_atm(bytes(row)) for row in payload),
                      dtype=np.uint8, count=packets)
    frame = np.concatenate([payload, crc[:, None]], axis=1)
    tx_bits = _bits_from_bytes(frame).astype(bool)
    means = background_mean + signal_mean * tx_bits[..., None]
    counts = rng.poisson(means, size=means.shape[:2] + (repetitions,)).sum(axis=2)
    threshold = map_threshold(repetitions * signal_mean,
                              repetitions * background_mean)
    rx_bits = counts >= threshold
    rx_bytes = np.packbits(rx_bits, axis=1, bitorder="big")
    rx_payload = rx_bytes[:, :5]
    rx_crc = rx_bytes[:, 5]
    valid_crc = np.fromiter(
        (crc8_atm(bytes(row)) == check for row, check in zip(rx_payload, rx_crc)),
        dtype=bool, count=packets
    )
    correct = np.all(rx_payload == payload, axis=1) & valid_crc
    undetected = valid_crc & np.any(rx_payload != payload, axis=1)
    return {
        "packets": float(packets),
        "packet_success": float(np.mean(correct)),
        "crc_acceptance": float(np.mean(valid_crc)),
        "undetected_wrong_payload": float(np.mean(undetected)),
        "payload_bits": 40.0,
        "frame_bits": 48.0,
        "repetitions": float(repetitions),
    }


def beam_pot_per_second(pot_per_year: float) -> float:
    return pot_per_year / SECONDS_PER_YEAR


def proton_beam_energy_j_per_pot(proton_energy_gev: float) -> float:
    return proton_energy_gev * 1.602176634e-10


def slot_duration_for_signal_mean(signal_mean: float, selected_events_per_second: float) -> float:
    if selected_events_per_second <= 0:
        return float("inf")
    return signal_mean / selected_events_per_second


def format_float(value: float) -> str:
    return f"{value:.10g}"
