import numpy as np
import matplotlib.pyplot as plt
import math as m
import netCDF4 as nc
np.set_printoptions(suppress=True)

filename = "energy_spectrum.dat"
energy_spectrum = nc.Dataset(filename)
energy_spectrum = np.array(energy_spectrum["grid data"][:], dtype=np.float64)

ev_to_ry = 13.605693122
n_atoms = 256

energy_spectrum = energy_spectrum/(n_atoms/(ev_to_ry*1000))

plt.figure()
plt.eventplot(energy_spectrum, orientation='horizontal', colors='b', linewidths=0.1)
plt.xlabel("Energy (meV/atom)")
plt.show()
