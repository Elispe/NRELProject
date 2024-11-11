import numpy as np
from math import pi


#Power system parameters
tau_z = 10
M_DERS = 40
D_DERS = 1.5
beta = -0.1
Kg = 1/0.05
Tsm = 2
ohm_per_km = 0.3
omega_s = 2 * pi * 60
SBASE = 1000  # MW
VBASE = 380  # kV
ZBASE = (VBASE * 1e3) ** 2 / (SBASE * 1e6)