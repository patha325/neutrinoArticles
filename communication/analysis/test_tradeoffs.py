import unittest
from tradeoffs import required_mw, detector_kt, latency_window_ms

class TradeoffChecks(unittest.TestCase):
    def test_inverse_power_mass_relation(self):
        self.assertAlmostEqual(detector_kt(5000,required_mw(5000,10)),10)
        self.assertAlmostEqual(detector_kt(5000,100),2.4400906770596196)

    def test_shorter_symbol_costs_more_power(self):
        self.assertAlmostEqual(required_mw(5000,10,.001)/required_mw(5000,10,1),1000)

    def test_propagation_window(self):
        self.assertAlmostEqual(latency_window_ms(5000),8.476590041790178)

if __name__=='__main__':unittest.main()
