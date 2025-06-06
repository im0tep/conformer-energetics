#!/usr/bin/env python3

import numpy as np
forces_mace = np.loadtxt("forces_mace_lpg.csv", delimiter=",", skiprows=1)  # eV/angstrom from mace atoms.get_forces() 
forces_psi4 = np.loadtxt("forces_psi4.csv", delimiter=",", skiprows=1)    # Hartree/Bohr from psi4 .gradient() 

forces_psi4_ev_ang= forces_psi4 / 51.42206747627547 #convert Hartree/Bohr into eV / angstrom 
rmse_mev = 1000 * np.sqrt(np.mean((forces_mace - forces_psi4_ev_ang)**2))
print ('Force RMSE is: ' + str (rmse_mev)  + 'in meV / angstrom') 
