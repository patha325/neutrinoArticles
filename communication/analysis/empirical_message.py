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

EVENTS_USED_FOR_MEAN = 1402
RECORDS_USED_FOR_MEAN = 3454
ON_PULSES_ESTIMATED = RECORDS_USED_FOR_MEAN / 2
SELECTED_PER_ON_PULSE = 2 * EVENTS_USED_FOR_MEAN / RECORDS_USED_FOR_MEAN
PROTONS_PER_ON_PULSE = 2.25e13
PROTON_GEV = 120.0
GEV_J = 1.602176634e-10
SLOTS_PER_SUPERCYCLE = 25
SUPERCYCLE_S = 61.267
REGULAR_SPACING_S = 2.2
POST_TRAIN_GAP_S = 6.267
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


def slot_time(slot_index):
    """NuMI pulse timing: 25 slots at 2.2 s, then the 6.267 s train gap."""
    cycle, phase = divmod(slot_index, SLOTS_PER_SUPERCYCLE)
    return cycle * SUPERCYCLE_S + phase * REGULAR_SPACING_S


def packet_latency_stats(repeats):
    """Latency stats for issue time uniform in a supercycle; includes phase wait."""
    slots = BITS_PER_PACKET * repeats
    gaps = [
        slot_time(i + 1) - slot_time(i)
        for i in range(SLOTS_PER_SUPERCYCLE)
    ]
    weighted_mean = 0.0
    minimum = math.inf
    maximum = 0.0
    for phase in range(SLOTS_PER_SUPERCYCLE):
        previous_gap = gaps[(phase - 1) % SLOTS_PER_SUPERCYCLE]
        probability = previous_gap / SUPERCYCLE_S
        wait_mean = previous_gap / 2
        transmission = slot_time(phase + slots - 1) - slot_time(phase)
        weighted_mean += probability * (wait_mean + transmission)
        minimum = min(minimum, transmission)
        maximum = max(maximum, previous_gap + transmission)
    return weighted_mean, minimum, maximum


def simulate_packets(repeats, trials, seed, signal_mean):
    rng = random.Random(seed + repeats)
    correct = accepted = undetected = on_pulses = 0
    zero_probability = math.exp(-signal_mean * repeats)
    for _ in range(trials):
        payload = rng.getrandbits(PAYLOAD_BYTES * 8).to_bytes(PAYLOAD_BYTES, "big")
        exact, check, pulses = send_packet(payload, repeats, rng, zero_probability)
        correct += exact
        accepted += check
        undetected += check and not exact
        on_pulses += pulses
    return correct / trials, accepted / trials, undetected / trials, on_pulses / trials


def scenario(repeats, trials=50000, seed=20260923, signal_mean=SELECTED_PER_ON_PULSE):
    correct, accepted, undetected, on_pulses = simulate_packets(
        repeats, trials, seed, signal_mean
    )
    latency, latency_min, latency_max = packet_latency_stats(repeats)
    proton_j_per_on_pulse = PROTONS_PER_ON_PULSE * PROTON_GEV * GEV_J
    return (
        repeats,
        trials,
        signal_mean,
        .5 * math.exp(-signal_mean * repeats),
        correct,
        accepted,
        undetected,
        latency,
        latency_min,
        latency_max,
        PAYLOAD_BYTES * 8 * (accepted - undetected) / latency,
        on_pulses * proton_j_per_on_pulse,
    )


HEADER = (
    "on_pulses_per_bit", "simulated_packets", "selected_events_per_on_pulse",
    "ideal_uncoded_BER",
    "exact_payload_fraction", "CRC_accepted_fraction",
    "undetected_wrong_payload_fraction", "mean_latency_s_uniform_issue_phase",
    "minimum_latency_s", "maximum_latency_s",
    "CRC_accepted_correct_payload_bits_per_second", "mean_proton_beam_joules_per_attempt",
    "input_mean_95pct_low", "input_mean_95pct_high",
    "correct_packet_fraction_at_input_low", "correct_packet_fraction_at_input_high",
    "accepted_correct_bps_at_input_low", "accepted_correct_bps_at_input_high",
)


def main():
    output = Path(__file__).resolve().parent.parent / "data" / "empirical_messages.csv"
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(HEADER)
        # Stancil et al. estimate lambda = 2*1402/3454, assuming half the
        # recorded bits are on. Propagate its approximate Poisson counting
        # interval separately from packet Monte Carlo sampling uncertainty.
        mean = SELECTED_PER_ON_PULSE
        se = math.sqrt(EVENTS_USED_FOR_MEAN) / ON_PULSES_ESTIMATED
        low, high = mean - 1.96 * se, mean + 1.96 * se
        for n in (3, 5, 9):
            baseline = scenario(n, signal_mean=mean)
            lower = scenario(n, trials=50000, seed=20261003, signal_mean=low)
            upper = scenario(n, trials=50000, seed=20261013, signal_mean=high)
            writer.writerow(baseline + (low, high, lower[5] - lower[6], upper[5] - upper[6], lower[10], upper[10]))
    print(output)


if __name__ == "__main__":
    main()
