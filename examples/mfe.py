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

# set the threshold for spin character to 5%
threshold = 0.05

Bs = np.linspace(0, 0.25, 250)
mfe = np.zeros_like(Bs)

Cs = np.zeros_like(Bs)
Ct = np.zeros_like(Bs)
Cq = np.zeros_like(Bs)


# set energy levels of T1 and T2 states
ET1 = 1.14
ET2 = 2.43

# set up the model
m = Model()
# energy gaps and IC rates
m.kIC1 = m.calculate_internal_conversion_rate(-ET1)
m.kIC2 = m.calculate_internal_conversion_rate(ET2-2*ET1)
m.kIC21 = m.calculate_internal_conversion_rate(ET1-ET2)
m.kRISC = 0
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

# loop through the magnetic fields and calculate the singlet population at each
for i, B in enumerate(Bs):

    sh.B = B
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.simulate()
    mfe[i] = m.S
    
    # get the number of TT states witn spin character above threshold
    Cs[i] = np.where(sh.cslsq > threshold)[0].shape[0]
    Ct[i] = np.where(sh.sum_ctlsq > threshold)[0].shape[0]
    Cq[i] = np.where(sh.sum_cqlsq > threshold)[0].shape[0]

# combine the results into a pandas DataFrame
results = pd.DataFrame(index=Bs*1000, columns=['MFE', 'CS', 'CT', 'CQ'], data=np.transpose(np.vstack((mfe, Cs, Ct, Cq))))
results = results.apply(lambda x: 100*(x-x[0])/x[0] if x.name in ['MFE'] else x)
results.index.name = 'B (mT)'

# plot the results
fig, ax = plt.subplots()
ax.axhline(0, color='0.5')
ax.plot(results.index, results['MFE'], 'k-')
ax.set_xlabel('Magnetic field strength (mT)')
ax.set_ylabel('MFE (%)')

fig, ax = plt.subplots()
ax.plot(Bs, Cs, label='> 5% S')
ax.plot(Bs, Ct, label='> 5% T')
ax.plot(Bs, Cq, label='> 5% Q')
ax.set_xlabel('Magnetic field strength (mT)')
ax.set_ylabel('# TT states')
ax.legend(frameon=False)

# save the results
#results.to_csv('mfe.csv', header=True, index=True)