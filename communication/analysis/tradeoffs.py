"""Derived receiver-power and idealized latency examples from model.py.
Values are conditional on its fixed-energy, fixed-divergence assumptions.
"""
import csv
import math
from pathlib import Path
from model import neutrino_energy_per_selected_event, TARGET_BER

OUT=Path(__file__).resolve().parent.parent/'data'
C=299_792_458.0
EARTH_R=6_371_000.0
FIBER_GROUP_INDEX=1.4677  # illustrative Corning SMF-28e+ 1550 nm product value

def required_mw(chord_km, detector_kt=10, symbol_s=1.0):
    return -math.log(2*TARGET_BER)*neutrino_energy_per_selected_event(chord_km,detector_kt)/symbol_s/1e6

def detector_kt(chord_km, neutrino_mw, symbol_s=1.0):
    return required_mw(chord_km,10,symbol_s)*10/neutrino_mw

def latency_window_ms(chord_km=5000):
    chord=chord_km*1000
    arc=2*EARTH_R*math.asin(chord/(2*EARTH_R))
    return (FIBER_GROUP_INDEX*arc-chord)/C*1000

def main():
    OUT.mkdir(exist_ok=True)
    with (OUT/'power_detector_tradeoff.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['chord_km','beam_half_angle_mrad','neutrino_energy_GeV','selected_efficiency','target_raw_BER','on_symbol_duration_s','neutrino_carried_on_symbol_MW','required_detector_kt'])
        for l in (1000,5000,12000):
            for p in (1,10,100):
                w.writerow([l,1,3,.5,TARGET_BER,1,p,round(detector_kt(l,p),6)])
    with (OUT/'latency_example.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['chord_km','detector_kt','neutrino_carried_on_symbol_MW','idealized_required_on_symbol_s','neutrino_propagation_ms','ideal_shortest_surface_fiber_ms','ideal_fiber_advantage_ms'])
        l=5000
        arc=2*EARTH_R*math.asin(l*1000/(2*EARTH_R))
        for p in (10,100,1000,3000,10000,30000):
            w.writerow([l,10,p,round(required_mw(l,10)/p,9),round(l*1000/C*1000,6),round(FIBER_GROUP_INDEX*arc/C*1000,6),round(latency_window_ms(l),6)])
    print('5000 km latency window (ms):',latency_window_ms())
    print('5000 km 10 kt neutrino MW for 1 ms on-symbol:',required_mw(5000,10,.001))

if __name__=='__main__':main()
