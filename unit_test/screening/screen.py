#!/usr/bin/env python3

import re
import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from pynucastro.rates import Nucleus, Rate
from pynucastro.networks import Composition

network = sys.argv[1]
inputs = sys.argv[2]

with open(f"../../networks/{network}/reaclib_rates.F90", "r") as file:
    
    reaclib_rates = file.read()
    
start = reaclib_rates.index("subroutine net_screening_init()")
end = reaclib_rates.index("end subroutine net_screening_init")
reaclib_rates = reaclib_rates[start:end]

scr_pat = re.compile(r"call add_screening_factor\(zion\(j([a-z0-9]{1,6})\), " +
    r"aion\(j\1\), &\s+zion\(j([a-z0-9]{1,6})\), aion\(j\2\)\)")

screened = scr_pat.findall(reaclib_rates)
screened = [(Nucleus(n1), Nucleus(n2)) for n1, n2 in screened]

with open(f"../../networks/{network}/actual_network.F90", "r") as file:
    
    actual_network = file.read()
    
start = actual_network.index("! Nuclides")
end = actual_network.index("! Reactions")
actual_network = actual_network[start:end]

nuclide_pat = re.compile(r"j([a-z0-9]{1,6})")
nuclides = nuclide_pat.findall(actual_network)
nuclides = list(map(Nucleus, nuclides))

with open(inputs, "r") as file:
    
    inputs_screen = file.read()

comp = Composition(nuclides)
mfrac_pat = re.compile(r"massfractions\(([0-9]+)\)\s+=\s+([0-9.ed]+)")

for i, mfrac in mfrac_pat.findall(inputs_screen):
    
    i = int(i) - 1
    mfrac = float(mfrac)
    nuc = nuclides[i].raw
    comp.set_nuc(nuc, mfrac)
    
fscr_pat = re.compile(r"fscr\(\s+[0-9]+\s+\):\s+([0-9.E+-]+)")

rhos = np.logspace(6.9, 9.3, num=50)
Ts = np.logspace(6.9, 9.3, num=50)
rhos, Ts = np.meshgrid(rhos, Ts)
rhos = rhos.flatten()
Ts = Ts.flatten()

Gammas = np.zeros_like(rhos)
aj1978 = np.zeros_like(rhos)
chug2007 = np.zeros_like(rhos)

for j in tqdm(range(len(rhos))):
    
    rho, T = rhos[j], Ts[j]
    inputs_screen = re.sub(r"density\s+=\s+[0-9.ed]+", f"density = {rho}", inputs_screen)
    inputs_screen = re.sub(r"temperature\s+=\s+[0-9.ed]+", f"temperature = {T}", inputs_screen)

    with open("__inputs_screen__", "w") as file:
        
        file.write(inputs_screen)
        
    proc = subprocess.run("./main.Linux.gfortran.exe __inputs_screen__".split(),
            stdout=subprocess.PIPE)

    fscr_micro = fscr_pat.findall(str(proc.stdout))
    fscr_micro = list(map(float, fscr_micro))
    fscr_micro = np.array(fscr_micro, dtype=np.float64)

    Zbar = comp.zbar()
    Abar = comp.abar()
    fscr_pnuc = np.zeros_like(fscr_micro)

    for i in range(len(fscr_pnuc)):
        
        nuc1, nuc2 = screened[i]
        fscr_pnuc[i] = Rate.screen_chugunov2007(nuc1, nuc2, rho, T, Zbar, Abar)

    # Some physical constants in CGS units
    amu = 1.6605390666e-24 # g
    qe = 4.80320425e-10 # esu
    hbar = 1.05457266e-27 # erg * s
    k_B = 1.38064903e-16 # erg / K

    # Average mass and total number density
    mbar = Abar * amu
    ntot = rho / mbar

    # Compute electron sphere radius
    a_e = (0.75 / (np.pi * Zbar * ntot))**(1/3)

    # Coulomb coupling factor
    Gamma_e = qe**2 / (a_e * k_B * T)

    rat = fscr_micro / fscr_pnuc
    rat[rat < 1] **= -1
    rat -= 1
    err = rat.mean()
    
    Gammas[j] = Gamma_e
    aj1978[j] = fscr_micro.mean()
    chug2007[j] = fscr_pnuc.mean()
    j += 1
    
    with open(f"fscr_{rho:.1e}_{T:.1e}.dat", "w") as file:
        
        print("AJ", "Chug2007", file=file)
        i = 0
        
        for item in zip(fscr_micro, fscr_pnuc):
            
            print(i, *item, file=file)
            i += 1
            
    sys.exit(0)
    
    """
    plt.plot(np.arange(len(fscr_micro)), fscr_micro, label="AJ")
    plt.plot(np.arange(len(fscr_pnuc)), fscr_pnuc, label="Chug2007")
    plt.title(f"ρ = {rho:.1e}, T = {T:.1e}")
    plt.legend()
    if rat.max() > 5: plt.yscale("log")

    plt.gcf().set_size_inches(15, 10)
    plt.savefig(f"plots/fscr_{rho:.1e}_{T:.1e}.png")
    plt.gcf().clear()
    """

indx = Gammas.argsort()
Gammas = Gammas[indx]
aj1978 = aj1978[indx]
chug2007 = chug2007[indx]

with open("fscr.dat", "w") as file:
    
    print("Gamma_e", "AJ", "Chug2007", file=file)
    
    for item in zip(Gammas, aj1978, chug2007):
        
        print(*item, file=file)

plt.plot(Gammas, aj1978, label="AJ")
plt.plot(Gammas, chug2007, label="Chug2007")
plt.title(r"F$_{scr}$ vs. $\Gamma_e$")
plt.xlabel(r"$\Gamma_e$")
plt.ylabel(r"F$_{scr}$")
plt.legend()
plt.yscale("log")
plt.show()
    
os.remove("__inputs_screen__")
