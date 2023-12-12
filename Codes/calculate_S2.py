# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:08:51 2023

@author: Walia Farzana
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
# Define the function for sech
def sech(x):
    return 2 / (np.exp(x) + np.exp(-x))

# Define the psi function
def psi(x, t, a, a1, x0, x1, b, b1):
    term1 = a * np.sqrt(2) * np.exp(1j * (b * x) / 2) * sech(a * (x - x0))
    term2 = a1 * np.sqrt(2) * np.exp(1j * (b1 * x) / 2) * sech(a1 * (x - x1))
    return term1 + term2

# Calculate x values
x_values = np.linspace(-100, 100, 10000)  # Define the range for x values

# Define the parameters
a = 1
a1 = 2  # a1=2*a
x0 = 400
x1 = 700
b = 3
b1 = -3  # b=-b1
# Calculate psi values for the x values
psi_values = psi(x_values, 0, 1, 2, 400, 700, 3, -3)

# Compute the gradient of psi with respect to x
gradient_psi = np.gradient(psi_values, x_values)

gradient_real_value=np.real(gradient_psi)
gradient_mean_value=np.abs(np.mean(gradient_real_value))**2
# Use the gradients as needed


# the function to get the absolute of psi (x,t=0)
def integrand(x, t, a, a1, x0, x1, b, b1):
    psi_value = psi(x, t, a, a1, x0, x1, b, b1)
    return (2*gradient_mean_value-0.5* np.abs(psi_value)**4)

## function to compute the definite integral from 0 to L
def calculate_S2(t, a, a1, x0, x1, b, b1, L):
    result, _ = quad(integrand, 0, L, args=(t, a, a1, x0, x1, b, b1))
    return result


# Calculate S0 for a specific time t and various L values
t = 0
L_values = [100, 500, 600, 900, 1000,5000,6000,7000,8000]  # Different values of L

S2_values = [] # empty list to append different value of S_0 integration for different values of L
for L in L_values:
    S2_t = calculate_S2(t, a, a1, x0, x1, b, b1, L)
    S2_values.append(S2_t)
    print(f"S0({t}) for L={L}: {S2_t}")
    
    
plt.figure(figsize=(8, 6))
plt.plot(L_values, S2_values, marker='o')
plt.xlabel('L values')
plt.ylabel('S2(t)')
plt.title('Numerical Solution for Different L Values')
plt.grid(True)
plt.show()