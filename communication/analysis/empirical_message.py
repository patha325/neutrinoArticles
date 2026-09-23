"""Finite-message illustration calibrated to the published NuMI–MINERvA run.

This is a NEW, externally synchronized, repeated-OOK/CRC-8 protocol, not a
reproduction of Stancil et al.'s convolutional code or measured packet data.
The signal mean is their selected-event count per proton pulse, including
events produced by neutrino interactions in upstream rock. Proton energy to
target is estimated; facility electricity and emitted-neutrino energy are not.

Run: python communication/analysis/empirical_message.py
"""
import csv
import math
import random
from pathlib import Path

SELECTED_PER_ON_PULSE = 0.81
PROTONS_PER_ON_PULSE = 2.25e13
PROTON_GEV = 120.0
GEV_J = 1.602176634e-10
SLOTS_PER_SUPERCYCLE = 25
SUPERCYCLE_S = 61.267
SLOT_S = SUPERCYCLE_S / SLOTS_PER_SUPERCYCLE
PAYLOAD_BYTES = 5
CRC_BYTES = 1
BITS_PER_PACKET = 8 * (PAYLOAD_BYTES + CRC_BYTES)


def crc8(data):
    """CRC-8/ATM: polynomial 0x07, init 0, no reflection, xorout 0."""
    remainder = 0
    for byte in data:
        remainder ^= byte
        for _ in range(8):
            remainder = (remainder << 1) ^ (0x07 if remainder & 0x80 else 0)
            remainder &= 0xff
    return remainder


def send_packet(payload, repeats, rng, zero_probability=None):
    """Return (exact_message, crc_acceptance, on_pulses_sent)."""
    if len(payload) != PAYLOAD_BYTES or repeats < 1:
        raise ValueError("expected five payload bytes and >=1 repeats")
    if zero_probability is None:
        zero_probability = math.exp(-SELECTED_PER_ON_PULSE * repeats)
    tx = payload + bytes([crc8(payload)])
    received = bytearray(tx)
    on_bits = 0
    for index, byte in enumerate(tx):
        for shift in range(8):
            mask = 1 << shift
            if byte & mask:
                on_bits += 1
                # Pool independent Poisson slots: P(no event in n on slots).
                if rng.random() < zero_probability:
                    received[index] &= ~mask
    rx_payload = bytes(received[:-1])
    correct = rx_payload == payload
    accepted = crc8(rx_payload) == received[-1]
    return correct, accepted, on_bits * repeats


def scenario(repeats, trials=50000, seed=20260923):
    rng = random.Random(seed + repeats)
    correct = accepted = undetected = on_pulses = 0
    for _ in range(trials):
        payload = rng.getrandbits(PAYLOAD_BYTES * 8).to_bytes(PAYLOAD_BYTES, "big")
        exact, check, pulses = send_packet(payload, repeats, rng)
        correct += exact
        accepted += check
        undetected += check and not exact
        on_pulses += pulses
    elapsed = BITS_PER_PACKET * repeats * SLOT_S
    proton_j_per_on_pulse = PROTONS_PER_ON_PULSE * PROTON_GEV * GEV_J
    return (
        repeats,
        trials,
        .5 * math.exp(-SELECTED_PER_ON_PULSE * repeats),
        correct / trials,
        accepted / trials,
        undetected / trials,
        elapsed,
        PAYLOAD_BYTES * 8 * (accepted - undetected) / trials / elapsed,
        on_pulses / trials * proton_j_per_on_pulse,
    )


HEADER = (
    "on_pulses_per_bit", "simulated_packets", "ideal_uncoded_BER",
    "exact_payload_fraction", "CRC_accepted_fraction",
    "undetected_wrong_payload_fraction", "seconds_per_attempt",
    "CRC_accepted_correct_payload_bits_per_second", "mean_proton_beam_joules_per_attempt",
)


def main():
    output = Path(__file__).resolve().parent.parent / "data" / "empirical_messages.csv"
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(HEADER)
        writer.writerows(scenario(n) for n in (3, 5, 9))
    print(output)


if __name__ == "__main__":
    main()
