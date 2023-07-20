from scipy import optimize
import math
import numpy as np
import matplotlib.pyplot as plt

def Energy(diego):
    diego = phi1, phi2, M
    return #Energy Equation

#Energy Equation
def calculate_energy(H, M, Phi_1, Phi_2, J_AF, E_ANIS):
    cos_phi_1 = math.cos(Phi_1)
    cos_phi_2 = math.cos(Phi_2)
    energy = -0.5 * H * M * (cos_phi_1 + cos_phi_2) - J_AF * math.cos(Phi_2 - Phi_1) + E_ANIS
    return energy


minimize(0


# Cubic Symmetry [001] Growth [Anistropic Terms]

'''

Case 1: H || [100] in-plane

'''

##K1 = 0.1 
##sin_squared_2phi_1 = math.sin(Phi_1 * 2) ** 2
##sin_squared_2phi_2 = math.sin(Phi_2 * 2)) ** 2
##E_ANIS = (1/8) * K1 * (sin_squared_2phi_1 + sin_squared_2phi_2)
##
##'''
##
##Case 2: H || [110] in-plane
##
##'''
##K1 = 0.1
##cos_squared_2phi_1 = math.cos(Phi_1 * 2) ** 2
##cos_squared_2phi_2 = math.cos(Phi_2 * 2)) ** 2
##E_ANIS = (1/8) * K1 * (cos_squared_2phi_1 + cos_squared_2phi_2)
##
##'''
##
##Case 3: H || [100] in-plane
##
##'''
##K1 = 0.1
##sin_squared_2phi_1 = math.sin(Phi_1 * 2) ** 2
##sin_squared_2phi_2 = math.sin(Phi_2 * 2)) ** 2
##sin_squared_4phi_1 = math.sin(Phi_1) ** 4
##sin_squared_4phi_2 = math.sin(Phi_2)) ** 4
##E_ANIS = (1/8) * K1 (sin_squared_2phi_1 + sin_squared_2phi_2 + sin_squared_4phi_1 + sin_squared_4phi_2)
##
##'''
##
##Case 4: H || [011_]
##
##'''
##K1 = 0.1
##sin_squared_2phi_1 = math.sin(Phi_1 * 2) ** 2
##sin_squared_2phi_2 = math.sin(Phi_2 * 2)) ** 2
##sin_squared_4phi_1 = math.sin(Phi_1) ** 4
##sin_squared_4phi_2 = math.sin(Phi_2)) ** 4
##E_ANIS = (1/8) * K1 (sin_squared_2phi_1 + sin_squared_2phi_2 + sin_squared_4phi_1 + sin_squared_4phi_2)
##
##
##'''
##
##To do: Other cases..
##'''
##

#For loop, with different H fields. 




