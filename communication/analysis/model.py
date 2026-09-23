"""Transparent, deliberately optimistic neutrino-link lower-bound model.
Run: python model.py. No third-party dependencies.
"""
import csv
import math
import random
from pathlib import Path

OUT = Path(__file__).parent
GEV_J = 1.602176634e-10
AMU_G = 1.66053906660e-24
E_GEV = 3.0
SIGMA_CM2 = 0.67e-38 * E_GEV  # illustrative neutrino CC per nucleon; energy-specific approximation
MASS_KT = 10.0
EFF = 0.5
THETA_RAD = 1e-3
TARGET_BER = .01
BACKGROUND = 0.0
F = -math.log(2*TARGET_BER)  # zero-background, threshold one, equally probable bits

def neutrino_energy_per_selected_event(l_km, mass_kt=MASS_KT, theta=THETA_RAD, eff=EFF):
    area_cm2 = math.pi * (theta*l_km*1e5)**2
    p = (mass_kt*1e9/AMU_G)*SIGMA_CM2/area_cm2*eff
    return E_GEV*GEV_J/p

def ber_no_background(signal):
    return .5*math.exp(-signal)

def monte_carlo_ber(signal, n=200_000, seed=20260923):
    rng=random.Random(seed)
    wrong=0
    for _ in range(n):
        if rng.getrandbits(1) and rng.random() < math.exp(-signal):
            wrong+=1
    return wrong/n

def main():
    # Baseline: Stancil et al., arXiv:1203.2847v2, on-pulse selected count 0.81.
    benchmark=[]
    for repeats in (1,5,9):
        lam=.81*repeats
        benchmark.append((repeats,lam,ber_no_background(lam),monte_carlo_ber(lam)))
    with (OUT/'benchmark.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['pooled_pulses','mean_signal','analytic_uncoded_BER','simulated_uncoded_BER']);w.writerows(benchmark)
    with (OUT/'sensitivity.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['chord_km','detector_kt','beam_divergence_mrad','selected_efficiency','minimum_neutrino_beam_MW_per_1bps','target_uncoded_BER'])
        for length in (1000,5000,12000):
            for mass in (1,10,40):
                for theta in (.1,1,10):
                    joules=F*neutrino_energy_per_selected_event(length,mass,theta*1e-3)
                    w.writerow([length,mass,theta,EFF,round(joules/1e6,8),TARGET_BER])
    for r in benchmark: print('benchmark',r)
    for l in (1000,5000,12000):
        print('lower bound',l, F*neutrino_energy_per_selected_event(l)/1e6,'MW neutrino beam at 1 bit/s')

if __name__=='__main__':main()
