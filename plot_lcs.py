#!/usr/bin/env python3

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

def load_snana_photometry(file_path, target_filter):
    times, mags, magerrs = [], [], []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("OBS:"):
                parts = line.split()
                if len(parts) < 9:
                    continue
                mjd = float(parts[1])
                flt = parts[2]
                mag = float(parts[5])
                magerr = float(parts[6])
                if flt == target_filter and mag != -999:
                    times.append(mjd)
                    mags.append(mag)
                    magerrs.append(None if magerr < 0 or magerr > 90 else magerr)
    return np.array(times), np.array(mags), np.array(magerrs)

def sanitize_filter_name(filter_name):
    return filter_name.replace(":", "")

def plot_all_lcs(base_dir, target_filter, max_day_offset=200):
    phot_paths = sorted(glob.glob(os.path.join(base_dir, "snana_virgo_data", "SN*", "Photometry", "*.snana.dat")))
    plt.figure(figsize=(12, 6))

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    marker_cycle = ['o', 's', '^', 'D', '*', 'v', 'P', 'X']
    n_colors = len(color_cycle)
    n_markers = len(marker_cycle)

    for i, path in enumerate(phot_paths):
        sn_name = os.path.basename(path).split(".")[0]
        times, mags, magerrs = load_snana_photometry(path, target_filter)
        if len(times) == 0:
            continue

        min_mjd = np.min(times)
        day_offsets = times - min_mjd
        mask = day_offsets <= max_day_offset
        if not np.all(mask):
            print(f"Skipping {np.sum(~mask)} late-time points (> {max_day_offset} days) for {sn_name}")
        times, mags, magerrs = times[mask], mags[mask], magerrs[mask]

        if len(times) < 1:
            print(f"No usable points to plot for {sn_name}")
            continue

        min_mag_idx = np.argmin(mags)
        t0 = times[min_mag_idx]
        t_shifted = times - t0

        color = color_cycle[i % n_colors]
        marker = marker_cycle[(i // n_colors) % n_markers]

        if np.any([e is not None for e in magerrs]):
            magerrs_clean = np.array([e if e is not None else 0 for e in magerrs])
            plt.errorbar(t_shifted, mags, yerr=magerrs_clean, fmt=marker, ms=4, alpha=0.6,
                         label=sn_name, color=color, linestyle='none')
        else:
            plt.plot(t_shifted, mags, marker, ms=4, alpha=0.6, label=sn_name, color=color)

    plt.gca().invert_yaxis()
    plt.xlabel("Days since brightest point", fontsize=14)
    plt.ylabel(f"{target_filter} magnitude", fontsize=14)
    #plt.title(f"Aligned SN Light Curves in {target_filter}", fontsize=16)
    plt.legend(fontsize=9, loc='upper right', ncol=2)
    plt.grid(True)
    plt.tight_layout()

    out_dir = os.path.join(base_dir, "plots")
    os.makedirs(out_dir, exist_ok=True)
    safe_filter = sanitize_filter_name(target_filter)
    out_path = os.path.join(out_dir, f"all_sne_{safe_filter}.pdf")
    plt.savefig(out_path)
    print(f"Saved figure to {out_path}")

    # Don't show the plot interactively
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot aligned SN light curves in a specific filter.")
    parser.add_argument("--filter", type=str, required=True, help="Filter name, e.g., 'standard::b'")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    plot_all_lcs(base_dir, args.filter)