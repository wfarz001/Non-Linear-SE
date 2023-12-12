# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:11:36 2023

@author: Walia Farzana
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from sympy import symbols, exp, sech, diff, integrate, conjugate, sqrt

# Define the symbols
x, t, a, b, x0, a1, b1, L = symbols('x t a b x0 a1 b1 L')

#function for sech, as python does not have built-in sech function
def sech(x):
    return 2 / (np.exp(x) + np.exp(-x))

#function to compute the psi(x,t=0) as per the reference paper (Dr. Vahala's Paper,equation(24))
def psi(x, t, a, a1, x0, x1, b, b1):
    term1 = a*np.sqrt(2) * np.exp(1j * (b * x)/2) * sech(a * (x - x0))
    term2 = a1*np.sqrt(2) * np.exp(1j * (b1 * x)/2) * sech(a1 * (x - x1))
    return term1 + term2

# the function to get the absolute of psi (x,t=0)
def integrand(x, t, a, a1, x0, x1, b, b1):
    psi_value = psi(x, t, a, a1, x0, x1, b, b1)
    return np.abs(psi_value)**2

## function to compute the definite integral from 0 to L
def calculate_S0(t, a, a1, x0, x1, b, b1, L):
    result, _ = quad(integrand, 0, L, args=(t, a, a1, x0, x1, b, b1))
    return result

# Set your values for a, a1, x0, x1, b, b1, L
a = 1
a1 = 2 # a1=2*a
x0 = 400
x1 = 700
b = 3
b1 = -3 #b=-b1# Plotting the numerical solution for different L values

# L = 6000


# # Calculate S0 for a specific time t
# t = 0 # As per the initial condition equation
# S0_t = calculate_S0(t, a, a1, x0, x1, b, b1, L)
# print(f"S0({t}) = {S0_t}")


# Calculate S0 for a specific time t and various L values
#t = 0
L_values = [100, 500, 600, 900, 1000,5000,6000,7000,8000]  # Different values of L

S0_values = [] # empty list to append different value of S_0 integration for different values of L
for L in L_values:
    S0_t = calculate_S0(t, a, a1, x0, x1, b, b1, L)
    print(f"S0({t}) for L={L}: {S0_t}")
    # Compute the time derivative of S0(t)
    dS0_dt = diff(S0_t, t)
    print('The result of the Derivation of S0 with respect to time:',dS0_dt)
    
    
# # Compute the time derivative of S0(t)
# dS0_dt = diff(S0_t, t)

# # Print the result
# print('The result of the Derivation of S0 with respect to time')
# print(dS0_dt)