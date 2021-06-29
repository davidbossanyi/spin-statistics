Spin Hamiltonian
================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
Introduction
------------

The Hamiltonian is the sum of the Zeeman, zero-field, intertriplet dipole-dipole and intertriplet exchange terms.
	
For details, refer to `this paper <https://doi.org/10.1021/acs.jpcc.6b04934>`_.

The class uses:

- :math:`\hbar=1`
- Tesla for magnetic field strength
- electronvolts for energies
- radians for angles
- ZX'Z'' convention for Euler angles
- the molecular system of molecule A throughout, which takes the x-axis parallel to the long molecular axis, the y-axis parallel to the short molecular axis and the z-axis perpendicular to the plane of the molecular backbone

The zero-field basis is used for the construction of the Hamiltonian:

.. math::

   \{|\phi_i\rangle\}=\{|xx\rangle, |xy\rangle, |xz\rangle, |yx\rangle, |yy\rangle, |yz\rangle, |zx\rangle, |zy\rangle, |zz\rangle\}

The :class:`spinstats.SpinHamiltonian` can be used to compute the total spin Hamiltonian, its eigenvectors and eigenvalues and the overlap factors that are required for further simulations. See the examples for usage.

API
---

.. autoclass:: spinstats.SpinHamiltonian
	:members: