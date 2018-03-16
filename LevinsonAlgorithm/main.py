'''
Author: Andac Demir
Date: Feb 14 2018

Description:
A function based on the Levinson algorithm, with the autocorrelation
sequence which is an ndarray of doubles as input, and with
ki(reflection coefficients), Ai and σ_i^2 (for 1 ≤ i ≤ M) as outputs.
M is the number of orders.
'''


import numpy as np


def Levinson(autocorr_sequence, M):
    # Initialize:
    r = np.array([autocorr_sequence[1]], dtype=complex)
    A = {}
    A[0] = np.array([1], dtype=complex)
    var = np.zeros(M+1, dtype=complex)
    var[0] = autocorr_sequence[0]
    k = np.zeros(M, dtype=complex)
    # Iterate
    for i in range(M):
        delta = np.dot(np.conjugate(r), np.flip(A[i], axis=0).T)
        k[i] = delta / var[i]
        A[i+1] = (np.append(A[i], 0) - k[i] * np.flip(
                                    np.append(np.conjugate(A[i]), 0), axis=0))
        var[i+1] = (1 - np.abs(k[i]) ** 2) * var[i]
        if i != M - 1:
            r = np.append(r, autocorr_sequence[i+2])
    return k, A, var


# Test with problem 8.10
autocorr_sequence = np.array([1.707, 0.5 + 0.5j, 0.707j], dtype=complex)
k, A, var = Levinson(autocorr_sequence, M=2)
print('Results for Problem 8.10:')
print('Ki (reflection coefficient values:)')
print(k)
print('A:')
print(A)
print('Variances')
print(var)
print('---------------------------------------------------------------------')

# Test with problem 8.11
autocorr_sequence = np.array([2.28, 0.768, 0.461], dtype=complex)
k, A, var = Levinson(autocorr_sequence, M=2)
print('Results for Problem 8.11:')
print('Ki (reflection coefficient values:)')
print(k)
print('A:')
print(A)
print('Variances')
print(var)
