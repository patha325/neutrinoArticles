"""Numerical invariants and limiting-case tests for the Paper II model."""

import unittest

import numpy as np

from neutrino_channel import (
    DuneBenchmark,
    binary_ook_capacity,
    crc8_atm,
    ook_error_probabilities,
    oscillation_probabilities,
    parse_reco_efficiency,
    pmns_matrix,
    ppm_symbol_error,
    prem_density_g_cm3,
    prem_path_segments,
    read_flux,
    read_xsec,
    selected_spectrum_per_pot,
    simulate_crc_packets,
)


class NeutrinoChannelTests(unittest.TestCase):
    def test_prem_profile_and_chord_are_physical(self):
        density = prem_density_g_cm3(np.array([0.0, 3480.0, 6371.0]))
        self.assertTrue(np.all(density > 0.0))
        lengths, segment_density = prem_path_segments(1284.9, 5.0)
        self.assertAlmostEqual(float(lengths.sum()), 1284.9, places=8)
        self.assertLessEqual(float(lengths.max()), 5.0 + 1e-10)
        self.assertTrue(np.all(segment_density > 0.0))

    def test_pmns_and_evolution_are_unitary(self):
        U = pmns_matrix(0.59, 0.15, 0.86, 0.4)
        np.testing.assert_allclose(U.conj().T @ U, np.eye(3), atol=1e-12)
        probs = oscillation_probabilities(
            np.array([0.8, 2.5, 5.0]), 1284.9, density_model="constant"
        )
        np.testing.assert_allclose(probs.sum(axis=1), np.ones((3, 3)), atol=1e-10)
        self.assertTrue(np.all((probs >= 0.0) & (probs <= 1.0)))

    def test_ook_background_free_limit_and_capacity(self):
        _, miss, ber = ook_error_probabilities(1.0, 0.0)
        self.assertAlmostEqual(miss, np.exp(-1.0), places=12)
        self.assertAlmostEqual(ber, 0.5 * np.exp(-1.0), places=12)
        capacity, p_on = binary_ook_capacity(1.0, 0.0, p_max=0.5)
        self.assertGreater(capacity, 0.0)
        self.assertGreaterEqual(p_on, 0.0)
        self.assertLessEqual(p_on, 0.5)

    def test_ppm_zero_signal_and_crc_vector(self):
        self.assertAlmostEqual(ppm_symbol_error(0.0, 0.0, 4), 0.75, places=12)
        self.assertEqual(crc8_atm(b"123456789"), 0xF4)

    def test_public_inputs_and_fold(self):
        flux, xsec = read_flux(), read_xsec()
        self.assertGreater(float(flux[:, 2].sum()), 0.0)
        self.assertEqual(len(parse_reco_efficiency()[0]), 80)
        one_tonne = selected_spectrum_per_pot(
            flux, xsec, baseline_km=DuneBenchmark().baseline_km,
            fiducial_mass_kt=0.001, density_model="constant"
        )["events_per_pot"]
        ten_tonne = selected_spectrum_per_pot(
            flux, xsec, baseline_km=DuneBenchmark().baseline_km,
            fiducial_mass_kt=0.01, density_model="constant"
        )["events_per_pot"]
        self.assertGreater(one_tonne, 0.0)
        self.assertAlmostEqual(ten_tonne / one_tonne, 10.0, places=10)

    def test_packet_simulation_is_reproducible(self):
        a = simulate_crc_packets(0.3, 0.01, repetitions=2, packets=200, seed=8)
        b = simulate_crc_packets(0.3, 0.01, repetitions=2, packets=200, seed=8)
        self.assertEqual(a, b)
        self.assertGreaterEqual(a["crc_acceptance"], a["packet_success"])
        self.assertLessEqual(a["undetected_wrong_payload"], a["crc_acceptance"])


if __name__ == "__main__":
    unittest.main()
