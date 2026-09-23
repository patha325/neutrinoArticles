import unittest
from throughput import raw_symbols_per_s,payload_bits_per_s

class ThroughputAccounting(unittest.TestCase):
    def test_normalization(self):
        self.assertAlmostEqual(raw_symbols_per_s(5000,10,100),4.098208355)
        self.assertAlmostEqual(payload_bits_per_s(5000,10,100),1.639283342,places=6)
    def test_chord_and_mass_scaling(self):
        self.assertAlmostEqual(raw_symbols_per_s(1000,10,1)/raw_symbols_per_s(5000,10,1),25)
        self.assertAlmostEqual(raw_symbols_per_s(5000,40,1)/raw_symbols_per_s(5000,10,1),4)
    def test_overhead_accounting(self):
        self.assertAlmostEqual(payload_bits_per_s(5000,10,100,.5,.8)/raw_symbols_per_s(5000,10,100),.4)
if __name__=='__main__':unittest.main()
