"""Illustrative accounting for raw OOK slots, FEC, framing and source compression.
FEC performance is NOT simulated: code rate and frame fraction are input assumptions.
"""
import csv
from pathlib import Path
from tradeoffs import required_mw

OUT=Path(__file__).resolve().parent.parent/'data'
CODE_RATE=.5
FRAME_FRACTION=.8

def raw_slots_per_s(chord_km,detector_kt,neutrino_on_mw):
    return neutrino_on_mw/required_mw(chord_km,detector_kt,1)

def payload_bits_per_s(chord_km,detector_kt,neutrino_on_mw,code_rate=CODE_RATE,frame_fraction=FRAME_FRACTION):
    return raw_slots_per_s(chord_km,detector_kt,neutrino_on_mw)*code_rate*frame_fraction

def main():
    with (OUT/'illustrative_throughput.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['chord_km','detector_kt','peak_neutrino_carried_on_slot_MW','raw_slots_per_s_at_1pct_uncoded_BER','assumed_FEC_rate','assumed_frame_payload_fraction','nominal_payload_bits_per_s_before_resends','illustrative_4_to_1_original_equivalent_bits_per_s'])
        for l in (1000,5000,12000):
            for m in (1,10,40):
                for p in (1,10,100):
                    raw=raw_slots_per_s(l,m,p);payload=payload_bits_per_s(l,m,p)
                    w.writerow([l,m,p,round(raw,8),CODE_RATE,FRAME_FRACTION,round(payload,8),round(payload*4,8)])
    print('5000 km, 10 kt, 100 MW peak on-slot power:',raw_slots_per_s(5000,10,100),'raw slots/s;',payload_bits_per_s(5000,10,100),'nominal payload bit/s')

if __name__=='__main__':main()
