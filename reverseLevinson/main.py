#A is a vector of doubles,
#first number corresponds to z^0, second z^-1, third z^-2 etc...
#You must enter 0 if a term does not exist.

import numpy as np
from fractions import Fraction

# For problem 1:
H = np.array([1, 2, Fraction(1,3)])
# For H1 in problem 2:
H1 = np.array([1, 0.4, -1.2, 2])
# For H2 in problem 3:
H2 = np.array([1, -0.6, 0.2, 0.5])


# Solves the lattice coefficients K recursively:
# Prints from K_n, K_n-1, ..., to K1.
def reverseLevinson(FIRfilter):
    filter_length = len(FIRfilter)

    if filter_length >= 2:
        k = FIRfilter[-1]

        if filter_length == 2 and FIRfilter[0] > 0:
            k = abs(k)

        print(k)
        FIRfilter = (FIRfilter - k * FIRfilter[::-1]) / (1 - k**2)
        FIRfilter = FIRfilter[:filter_length - 1]
        reverseLevinson(FIRfilter)


print("Lattice coefficients of problem 1")
reverseLevinson(H)

print("-----------------------------------")

print("Lattice coefficients of H1 in problem 2")
reverseLevinson(H1)

print("-----------------------------------")

print("Lattice coefficients of H2 in problem 2")
reverseLevinson(H2)