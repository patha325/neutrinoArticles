import math
import unittest

from empirical_message import (
    BITS_PER_PACKET, GEV_J, PROTONS_PER_ON_PULSE, PROTON_GEV,
    REGULAR_SPACING_S, SUPERCYCLE_S, crc8, packet_latency_stats,
    scenario, send_packet, slot_time,
)


class MessageCheck(unittest.TestCase):
    def test_published_crc8_reference_vector(self):
        self.assertEqual(crc8(b"123456789"), 0xF4)

    def test_no_missed_events_gives_exact_valid_message(self):
        import random
        self.assertEqual(
            send_packet(b"hello", 5, random.Random(7), zero_probability=0),
            (True, True, 5 * sum(byte.bit_count() for byte in b"hello" + bytes([crc8(b"hello")]))),
        )

    def test_all_ones_missed_is_undetected_by_crc(self):
        import random
        exact, accepted, _ = send_packet(
            b"hello", 3, random.Random(7), zero_probability=1
        )
        # A one-sided channel can turn an entire codeword into zero bytes,
        # which pass a zero-initial-state CRC; detection is not guaranteed.
        self.assertFalse(exact)
        self.assertTrue(accepted)

    def test_supercycle_and_power_bookkeeping(self):
        row = scenario(5, trials=1000)
        self.assertAlmostEqual(slot_time(25), SUPERCYCLE_S)
        self.assertAlmostEqual(slot_time(24), 24 * REGULAR_SPACING_S)
        self.assertAlmostEqual(row[7], 586.8866, places=3)
        self.assertLess(row[8], row[7])
        self.assertGreater(row[9], row[7])
        self.assertAlmostEqual(row[3], math.exp(-(.81 * 5)) / 2, places=3)
        # One on pulse sends 2.25e13 protons at 120 GeV: ~433 kJ.
        self.assertAlmostEqual(PROTONS_PER_ON_PULSE * PROTON_GEV * GEV_J / 1e3, 432.588, places=2)

    def test_packet_reliability_improves_with_repetition(self):
        rows = [scenario(n, trials=4000) for n in (3, 5, 9)]
        self.assertLess(rows[0][4], rows[1][4])
        self.assertLess(rows[1][4], rows[2][4])

    def test_latency_accounts_for_long_supercycle_gap(self):
        average, minimum, maximum = packet_latency_stats(1)
        self.assertGreater(average, BITS_PER_PACKET * REGULAR_SPACING_S)
        self.assertGreater(maximum - minimum, 6.0)


if __name__ == "__main__":
    unittest.main()
