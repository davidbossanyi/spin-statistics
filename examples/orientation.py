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

angles = np.linspace(0, np.pi, 250)
eta = np.zeros_like(angles)

Cs = np.zeros_like(angles)
Ct = np.zeros_like(angles)
Cq = np.zeros_like(angles)


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

# loop through the angles and calculate the eta value at each
for i, angle in enumerate(angles):

    sh.alpha = angle
    sh.calculate_everything()
    m.set_spinhamiltonian_parameters(sh)
    m.calculate_eta()
    eta[i] = m.eta
    
    # get the number of TT states witn spin character above threshold
    Cs[i] = np.where(sh.cslsq > threshold)[0].shape[0]
    Ct[i] = np.where(sh.sum_ctlsq > threshold)[0].shape[0]
    Cq[i] = np.where(sh.sum_cqlsq > threshold)[0].shape[0]

# combine the results into a pandas DataFrame
results = pd.DataFrame(index=angles, columns=['eta', 'CS', 'CT', 'CQ'], data=np.transpose(np.vstack((eta, Cs, Ct, Cq))))
results = results.apply(lambda x: 100*x if x.name in ['eta'] else x)
results.index.name = 'J (eV)'

# plot the results
fig, ax = plt.subplots()
ax.plot(results.index, results['eta'], 'k-')
ax.set_xlabel('Angle (radians)')
ax.set_ylabel(r'$\eta$ (%)')

fig, ax = plt.subplots()
ax.plot(angles, Cs, label='> 5% S')
ax.plot(angles, Ct, label='> 5% T')
ax.plot(angles, Cq, label='> 5% Q')
ax.set_xlabel('Angle (radians)')
ax.set_ylabel('# TT states')
ax.legend(frameon=False)

# save the results
#results.to_csv('angle.csv', header=True, index=True)