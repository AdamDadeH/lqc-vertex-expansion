"""Redo one line of the Bianchi 1 Spectrum/Spect scan: shoot the Theta recurrence for a few trial
eigenvalues e and plot |a_n| (cf. Recur.nb with p1=50, p2=100, where the log puts the transition near e=-5001.4).

Usage:  uv run examples/bianchi1_spectrum_scan.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lqc import bianchi1

p1, p2, T = 50, 100, 50000
fig, ax = plt.subplots(figsize=(8, 5))
for e in (-5000.0, -5001.0, -5001.4, -5002.0, -5005.0):
    psi = bianchi1.shoot(e, p1, p2, T, scale=1.0)
    n = np.arange(T + 1)
    ax.plot(n[1000:], np.abs(psi[1000:]) / np.abs(psi[1000]), lw=0.8, label=f"e={e}  ratio={bianchi1.envelope_ratio(psi):.2f}")
ax.set_xscale("log"); ax.set_xlabel("n"); ax.set_ylabel("|a_n| / |a_1000|"); ax.legend()
ax.set_title("Bianchi I Theta recurrence, p1=50 p2=100: damped vs divergent trial eigenvalues")
fig.tight_layout(); fig.savefig("bianchi1_spectrum_scan.png", dpi=130)
print("wrote bianchi1_spectrum_scan.png")
