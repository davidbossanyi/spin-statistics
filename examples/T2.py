from spinstats.models import Model2 as Model
from spinstats import SpinHamiltonian
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# set spin Hamiltonian parameters
sh = SpinHamiltonian()
sh.J = 0
sh.D = 6.45e-6
sh.E = -6.45e-7
sh.X = 6.45e-9
sh.rAB = (0, 0, 0)
sh.alpha, sh.beta, sh.gamma, sh.theta, sh.phi, sh.B = 0, 0, 0, 0, 0, 0
sh.calculate_everything()


# set up the model
m = Model()
# excited states
m.kSF = 100
m.k_SF = 100
m.kTS = 10
m.kTF = 10/2
m.kD = 0.1
m.kTTA = 1
# rates of decay
m.kS = 0.0625
m.kT = 0
m.initial_weighting = {'T': 1}
# gaplaw parameters
m.gaplaw_A = 4900
m.gaplaw_gamma = 0.845
m.gaplaw_Evib = 0.17
m.gaplaw_kT = 0.025

# set the energy levels of T1 (fixed) and T2 varied
ET1 = 1.14
ET2s = np.linspace(2.36, 2.5, 100)

# initialise a results DataFrame
results = pd.DataFrame(index=ET2s, columns=['perp, no RISC', 'perp, RISC', 'para, no RISC', 'oblique, no RISC', 'para, RISC'])

# loop through the T2 energies and calculate the eta value at each in 5 different cases
for ET2 in ET2s:

    # set internal conversion rates
    m.kIC1 = m.calculate_internal_conversion_rate(-ET1)
    m.kIC2 = m.calculate_internal_conversion_rate(ET2-2*ET1)
    m.kIC21 = m.calculate_internal_conversion_rate(ET1-ET2)

    # perp, no RISC
    m.kRISC = 0
    m.kSF = 0
                
    sh.alpha = np.pi/2
    sh.beta = 0
    sh.gamma = 0
    
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    
    results.loc[ET2, 'perp, no RISC'] = m.eta
    
    # perp, RISC
    m.kRISC = 5000
    m.kSF = 100
                
    sh.alpha = np.pi/2
    sh.beta = 0
    sh.gamma = 0
    
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    
    results.loc[ET2, 'perp, RISC'] = m.eta
    
    # para, no RISC
    m.kRISC = 0
    m.kSF = 0
                
    sh.alpha = 0
    sh.beta = 0
    sh.gamma = 0
    
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    
    results.loc[ET2, 'para, no RISC'] = m.eta
    
    # oblique, no RISC
    m.kRISC = 0
    m.kSF = 0
                
    sh.alpha = 6.045128780774355
    sh.beta = 2.014966076324582
    sh.gamma = 3.3797107890263254
    
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    
    results.loc[ET2, 'oblique, no RISC'] = m.eta
    
    # para, RISC
    m.kRISC = 5000
    m.kSF = 100
                
    sh.alpha = 0
    sh.beta = 0
    sh.gamma = 0
    
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    
    results.loc[ET2, 'para, RISC'] = m.eta


# subtract TT energy
results.index -= 2*ET1

# plot the results
fig, ax = plt.subplots()
results *= 100
results.plot(ax=ax)
ax.set_xlabel('TT-T2 energy gap (eV)')
ax.set_ylabel(r'$\eta$ (%)')
ax.legend(frameon=False)

# save the results
# results.to_csv('T2.csv', header=True, index=True)