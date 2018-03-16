# For an autocorrelation sequence determines the lattice coefficients
# and prediction variances.

import numpy as np

# FIR filter of order 20:
autocorr_sequence = np.zeros(21)
autocorr_sequence[0] = 1
autocorr_sequence[1:21] = 0.5

def ShurAlgorithm(autocorr_sequence):
    FIRfilter_order = len(autocorr_sequence) - 1
    lattice_coefs = np.zeros(FIRfilter_order)
    prediction_vars = np.zeros(FIRfilter_order + 1)
    G = np.zeros((2, FIRfilter_order + 1))
    A = np.zeros((2, FIRfilter_order + 1))
    G[0], G[1] = autocorr_sequence, autocorr_sequence
    G[0,0] = 0

    for i in range(FIRfilter_order):
        lattice_coefs[i] = G[0, i+1] / G[1, i]
        prediction_vars[i] = G[1, i]
        A[0, :] = G[0, :]
        A[1, i+1:] = G[1, i:FIRfilter_order]
        G = np.matmul(np.array(([1, -lattice_coefs[i]], [-lattice_coefs[i], 1])), A)

    prediction_vars[FIRfilter_order] = G[1, FIRfilter_order]

    return lattice_coefs, prediction_vars


lattice_coefs, prediction_vars = ShurAlgorithm(autocorr_sequence)

print("These are the lattice coefficients:")
print(lattice_coefs)

print("--------------------------------------")

print("These are the prediction variances:")
print(prediction_vars)