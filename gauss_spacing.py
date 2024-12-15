"""
A script to plot and evaluate gaussian width and spacing.

See the module docs in `gaussian.py` for conclusions drawn based on this.
"""

import numpy as np
from matplotlib import pyplot as plt

NUM_X = int(1e4)
# Range to plot
X_RANGE = (0., 12.)

GAUSS_START = 2.
GAUSS_SPACING = 0.4

NUM_GAUS = 8

A = 1.
C_COEFF = 0.6
# C_COEFF = 0.55
C = GAUSS_SPACING * C_COEFF

# Used to offset the amplitude increase that arrises from combining gaussians.
# AMP_SCALER = 1. - C_COEFF**2
AMP_SCALER = 0.6649  # for coeff = 0.6
# AMP_SCALER = 0.7253 # for coeff = 0.55

def gauss(ctr: float, a: float, c: float, posit_eval: float) -> float:
    x = posit_eval - ctr
    return a * np.exp(-(x ** 2) / (2. * c ** 2))


def main():
    """Plot gaussians"""
    x = np.linspace(X_RANGE[0], X_RANGE[1], NUM_X)
    result = np.zeros(NUM_X)

    # Generate evenly spaced Gaussian centers
    gauss_centers = np.linspace(GAUSS_START, GAUSS_START + GAUSS_SPACING * (NUM_GAUS - 1), NUM_GAUS)

    # Sum up the gaussians
    for gauss_ctr in gauss_centers:
        result += gauss(gauss_ctr, A, C, x) * AMP_SCALER

    # Plot the result
    plt.plot(x, result)
    plt.title("Sum of Uniformly Spaced Gaussians")
    plt.xlabel("x")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.show()


main()