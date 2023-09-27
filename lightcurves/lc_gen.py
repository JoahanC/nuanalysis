import numpy as np
import matplotlib.pyplot as plt


def gen_background(time_domain, amp):
    bkg = amp / 15 * np.ones(len(time_domain))
    bkg_spread = amp / 9 * np.random.rand(len(time_domain))
    bkg_error_bars = amp / 5 / 4 * np.ones(len(time_domain))
    bkg_data = bkg + bkg_spread

    return bkg_data, bkg_error_bars


def gen_pulse(time_domain, amp, peak, sigma, true_gamma):

    # Transient peak
    peak_space = np.copy(time_domain)
    peak_mask = np.copy(time_domain)
    peak_mask[peak_mask < peak - sigma] = 0
    peak_mask[peak_mask > peak + sigma] = 0

    peak_data = amp / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((peak - peak_space + sigma) / 0.2 * sigma) ** 2) / 1.5
    peak_spread = np.multiply(peak_data / 2 / sigma, np.random.rand(len(time_domain)))
    peak_data += peak_spread
    peak_data[peak_data > amp] = 0

    # Adding transient tail

    decay_mask = np.copy(time_domain)
    decay_mask[decay_mask < peak + sigma] = 0
    decay_mask[decay_mask >= peak + sigma] = 1

    decay_data = np.multiply(amp * (time_domain - peak) ** (-true_gamma), decay_mask)
    where_are_NaNs = np.isnan(decay_data)
    decay_data[where_are_NaNs] = 0
    decay_data[decay_data > amp] = 0
    decay_spread = np.copy(decay_data)
    decay_spread = np.multiply(decay_spread, np.random.rand(len(time_domain)))
    decay_data += decay_spread

    return decay_data + peak_data


def generate_gaps(time_domain):
    interval_gap = np.max(time_domain) / 10
    weights = 100 * np.random.uniform(-1, 1, 10)
    intervals = list(range(10) * interval_gap + weights)
    
    blank_ranges = []
    for interval in intervals:
        interval_len = np.random.rand(1) * 30
        interval_min = interval - interval_len 
        interval_max = interval + interval_len
        if interval_min < 0:
            interval_min = 0
        if interval_max > np.max(time_domain):
            interval_max = 1000
        blank_ranges.append([interval_min, interval_max])
    
    gap_domain = np.ones(len(time_domain))
    for bound in blank_ranges:
        gaps = np.multiply(np.logical_and(time_domain > bound[0], time_domain < bound[1]), np.ones(len(time_domain)))
        gap_domain -= gaps
    gap_domain[gap_domain < 0] = 0
    return gap_domain


def gen_multicurve(time_domain, volatility, variability, n):

    random_amps = np.random.rand(n) * volatility
    random_peaks = np.random.rand(n) * np.max(time_domain)
    random_sigmas = np.random.rand(n) * variability
    random_gammas = 2 * np.random.rand(n) + 1

    data, bkg_error = gen_background(time_domain, 10)
    gap_domain = generate_gaps(time_domain)
    
    for idx, val in enumerate(random_amps):
        pulse_data = gen_pulse(time_domain, random_amps[idx], random_peaks[idx], random_sigmas[idx], random_gammas[idx])
        data += pulse_data
    
    return np.multiply(data, gap_domain), np.multiply(bkg_error, gap_domain)


