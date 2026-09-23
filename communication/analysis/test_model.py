import math
import unittest

from model import ber_no_background, monte_carlo_ber, neutrino_energy_per_selected_event


class ModelChecks(unittest.TestCase):
    def test_benchmark_zero_count_relation(self):
        self.assertAlmostEqual(ber_no_background(.81 * 5), math.exp(-4.05) / 2)
        self.assertLess(abs(monte_carlo_ber(.81 * 5) - ber_no_background(.81 * 5)), .002)

    def test_far_field_scaling(self):
        base = neutrino_energy_per_selected_event(1000, 10, .001)
        self.assertAlmostEqual(neutrino_energy_per_selected_event(5000, 10, .001) / base, 25)
        self.assertAlmostEqual(neutrino_energy_per_selected_event(1000, 40, .001) / base, .25)
        self.assertAlmostEqual(neutrino_energy_per_selected_event(1000, 10, .01) / base, 100)

    def test_one_percent_target(self):
        self.assertAlmostEqual(ber_no_background(-math.log(.02)), .01)


if __name__ == '__main__':
    unittest.main()
