"""Geometric extrapolation of the manuscript's idealized neutrino count model.

Run from repository root: python communication/analysis/solar_system_scaling.py
The output is neutrino-carried, full-duty-equivalent on-symbol power; no
source conversion, attenuation, oscillations, backgrounds, or coding are used.
"""
import csv
from pathlib import Path

from model import F, neutrino_energy_per_selected_event

C_KM_S = 299792.458
AU_KM = 149597870.7
BASELINES_KM = (
    ("terrestrial_reference", 5000.0),
    ("mean_earth_moon", 384400.0),
    ("one_au_illustrative", AU_KM),
    ("two_au_illustrative", 2 * AU_KM),
)


def rows():
    for label, distance in BASELINES_KM:
        yield (
            label,
            distance,
            F * neutrino_energy_per_selected_event(distance) / 1e6,
            distance / C_KM_S,
        )


def main():
    output = Path(__file__).resolve().parent.parent / "data" / "solar_system_scaling.csv"
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scenario", "baseline_km", "neutrino_on_symbol_MW_for_1_raw_symbol_s", "vacuum_flight_s"))
        writer.writerows(rows())
    print(output)


if __name__ == "__main__":
    main()
