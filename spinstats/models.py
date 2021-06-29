import numpy as np


class _RateModel:

    def __init__(self):
        self._number_of_states = 2
        self.states = ['S', 'T']
        self.rates = []
        self.model_name = 'base'
        self._time_resolved = False
        self.G = 2.7e13
        self._allowed_initial_states = {'S', 'T'}
        self._initial_state_mapping = {'S': 0, 'T': -1}
        self.initial_weighting = {'S': 1}
        self.gaplaw_A = 4900
        self.gaplaw_gamma = 0.845
        self.gaplaw_Evib = 0.17
        self.gaplaw_kT = 0.025
        # couplings
        self.cs = np.ones(9)/9
        self.ct_x = np.ones(9)/9
        self.ct_y = np.ones(9)/9
        self.ct_z = np.ones(9)/9
        self.cq_a = np.ones(9)/9
        self.cq_b = np.ones(9)/9
        self.cq_x = np.ones(9)/9
        self.cq_y = np.ones(9)/9
        self.cq_z = np.ones(9)/9
        
    def _check_initial_weighting(self):
        for starting_state in self.initial_weighting.keys():
            if starting_state not in self._allowed_initial_states:
                raise ValueError('invalid state {0} in initial_weighting'.format(starting_state))
            if self.initial_weighting[starting_state] < 0:
                raise ValueError('weightings must be positive')
        return
            
    def _set_initial_condition(self):
        self._y0 = np.zeros(self._number_of_states)
        total_weights = np.sum(np.array(list(self.initial_weighting.values())))
        for key in self.initial_weighting.keys():
            idx = self._initial_state_mapping[key]
            weight = self.initial_weighting[key]/total_weights
            self._y0[idx] = weight*self.G
        return
    
    def calculate_internal_conversion_rate(self, deltaE):
        """
        Compute an internal conversion rate based on the energy gap law.

        Parameters
        ----------
        deltaE : float
            Energy gap between the initial and final state. Positive if the final state is higher in energy.

        Returns
        -------
        kIC : float
            Internal conversion rate in unts of per ns.

        """
        if deltaE >= 0:
            kIC = self.gaplaw_A*np.exp(-deltaE/self.gaplaw_kT)
        else:
            kIC = self.gaplaw_A*np.exp(deltaE*self.gaplaw_gamma/self.gaplaw_Evib)
        return kIC
    
    def set_spinhamiltonian_parameters(self, spinhamiltonian):
        """
        Import calculated overlaps from the spin Hamiltonian

        Parameters
        ----------
        spinhamiltonian : spinstats.SpinHamiltonian
            Instance of :class:`spinstats.SpinHamiltonian`, with parameters already computed.

        Returns
        -------
        None.

        """
        self.cs = spinhamiltonian.cslsq
        self.ct_x = spinhamiltonian.ctlsq_x
        self.ct_y = spinhamiltonian.ctlsq_y
        self.ct_z = spinhamiltonian.ctlsq_z
        self.cq_a = spinhamiltonian.cqlsq_a
        self.cq_b = spinhamiltonian.cqlsq_b
        self.cq_x = spinhamiltonian.cqlsq_x
        self.cq_y = spinhamiltonian.cqlsq_y
        self.cq_z = spinhamiltonian.cqlsq_z
        return
    
                
class Model2(_RateModel):
    r"""
    A class for performing calculations of the spin-statistical factor.
    
    
    Attributes
    ----------
    states : list of str
        The names of the excited state species.
    rates : list of str
        The names of the different rate constants in the model.
    model_name : str
        The name of the model.
    initial_weighting : dict
        Dictionary of (str, float) pairs. Key is the state name (str) and value is its initial weight (float).
    gaplaw_A : float
        Energy gap law prefactor in unts of per ns.
    gaplaw_gamma : float
        Material dependent gamma factor from the energy gap law.
    gaplaw_Evib : float
        Vibrational energy to use in the energy gap law in units of eV.
    gaplaw_kT : float
        Thermal energy for internal conversion calculations in units of eV.
    G : float
        The exciton generation rate for :attr:`initial_state`. Units of per volume per time.
    kD : float
        Rate of triplet-pair dissociation into free triplets. Units of per ns.
    kTTA : float
        (Arbitrary) effective linear rate of triplet-triplet annihilation . UNits of per ns.
    kIC1 : float
        Rate of internal conversion from :math:`^3(TT)` to :math:`T_1`. Units of per ns.
    kIC2 : float
        Rate of internal conversion from :math:`^3(TT)` to :math:`T_2`. Units of per ns.
    kIC21 : float
        Rate of internal conversion from :math:`T_2` to :math:`T_1`. Units of per ns.
    kRISC : float
        Rate of high-level intersystem crossing from :math:`T_2` to :math:`S_1`. Units of per ns.
    kTF : float
        Rate of triplet-pair fusion from weakly exchange-coupled to strongly exchange-coupled. Units of per ns.
    kTS : float
        Rate of triplet-pair separation from strongly exchange-coupled to weakly exchange-coupled. Units of per ns.
    kSF : float
        Rate of singlet fission from :math:`S_1` to :math:`^1(TT)`. Units of per ns.
    k_SF : float
        Rate of triplet-pair fusion from :math:`^1(TT)` to :math`S_1`. Units of per ns.
    kS : float
        Radiative rate of the :math:`S_1` state. Units of per ns.
    kT : float
        Nonradiative rate of the :math:`T_1` state. Should be set to zero for true calculations of spin statistics. Units of per ns.
    cs : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^1(TT)`.
    ct_m : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^3(TT)_m`, where m is x, y or z.
    cq_m : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^5(TT)_m`, where m is a, b, x, y or z.
    simulation_results : dict
        Produced by :meth:`simulate`. Keys are the excited-state names (str), values are the simulated populations (float).
    plqy : float
        Calculated photoluminescence quantum yield from :meth:`calculate_plqy`.
    ucqy : float
        Calculated upconversion quantum yield from :meth:`calculate_ucqy`.
    eta : float
        Calculated spin-statistical factor from :meth:`calculate_eta`.
    """
    
    def __init__(self):
        super().__init__()
        # metadata
        self.model_name = 'Model 2'
        self._number_of_states = 21
        self.states = ['S', 'T']
        self.rates = ['kD', 'kTTA', 'kIC1', 'kIC2', 'kIC21', 'kRISC', 'kT', 'kTF', 'kTS', 'kSF', 'k_SF', 'kS']
        # rates between excited states
        self.kD = 0.1
        self.kTTA = 1
        self.kIC1 = 30
        self.kIC2 = 30
        self.kIC21 = 15
        self.kRISC = 1000
        self.kTF = 10/6
        self.kTS = 10
        self.kSF = 100
        self.k_SF = 100
        # rates of decay
        self.kS = 0.0625
        self.kT = 0
        # initial stuff
        self._allowed_initial_states = {'S', 'T'}
        self._initial_state_mapping = {'S': -2, 'T': 0}
        self.initial_weighting = {'T': 1}
        
    def _generate_rate_equation_matrix(self):
        #
        self.ct = [self.ct_x, self.ct_y, self.ct_z]
        self.cq = [self.cq_a, self.cq_b, self.cq_x, self.cq_y, self.cq_z]
        #
        self._rem = np.zeros((self._number_of_states, self._number_of_states))
        # T1
        self._rem[0, 0] = -(2*self.kTTA+self.kT)
        self._rem[0, 1:10] = 2*self.kD
        self._rem[0, 15:18] = self.kIC1
        self._rem[0, 20] = self.kIC21
        # TT_l
        for i in range(9):
            self._rem[i+1, i+1] = -(self.kD+self.kTF*(self.cs[i]+self.ct_x[i]+self.ct_y[i]+self.ct_z[i]+self.cq_a[i]+self.cq_b[i]+self.cq_x[i]+self.cq_y[i]+self.cq_z[i]))
            self._rem[i+1, 0] = (1/9)*self.kTTA
            for m in range(5):
                self._rem[i+1, m+10] = self.kTS*self.cq[m][i]
            for m in range(3):
                self._rem[i+1, m+15] = self.kTS*self.ct[m][i]
            self._rem[i+1, 18] = self.kTS*self.cs[i]
        # 5TT
        for m in range(5):
            self._rem[m+10, m+10] = -self.kTS*np.sum(self.cq[m])
            for i in range(9):
                self._rem[m+10, i+1] = self.kTF*self.cq[m][i]
        # 3TT
        for m in range(3):
            self._rem[m+15, m+15] = -(self.kTS*np.sum(self.ct[m])+self.kIC1+self.kIC2)
            for i in range(9):
                self._rem[m+15, i+1] = self.kTF*self.ct[m][i]
        # 1TT
        self._rem[18, 18] = -(self.kTS*np.sum(self.cs)+self.k_SF)
        for i in range(9):
            self._rem[18, i+1] = self.kTF*self.cs[i]
        self._rem[18, 19] = self.kSF
        # S1
        self._rem[19, 19] = -(self.kS+self.kSF)
        self._rem[19, 18] = self.k_SF
        self._rem[19, 20] = self.kRISC
        # T2
        self._rem[20, 20] = -(self.kIC21+self.kRISC)
        self._rem[20, 15:18] = self.kIC2
        return
    
    def _generate_generation_rate_matrix(self):
        self._grm = np.zeros((self._number_of_states, 1))
        self._grm[-2, 0] = -self._GS
        self._grm[0, 0] = -self._GT
        return
            
    def simulate(self):
        """
        Calculate the excited state populations.

        Returns
        -------
        None.

        """
        self._set_generation_rates()
        self._generate_generation_rate_matrix()
        self._generate_rate_equation_matrix()
        self.sol = np.linalg.solve(self._rem, self._grm)
        self._unpack_simulation()
        self._wrap_simulation_results()
        return
    
    def calculate_plqy(self):
        """
        Calculate the photoluminescence quantum yield.

        Returns
        -------
        None.

        """
        self.initial_weighting = {'S': 1}
        self.simulate()
        self.plqy = self.kS*self.S/self.G
        return
        
    def calculate_ucqy(self):
        """
        Calculate the upconversion quantum yield.

        Returns
        -------
        None.

        """
        self.initial_weighting = {'T': 1}
        self.simulate()
        self.ucqy = 2*self.kS*self.S/self.G
        return
    
    def calculate_eta(self):
        """
        Calculate the spin-statistical factor.

        Returns
        -------
        None.

        """
        self.calculate_plqy()
        self.calculate_ucqy()
        self.eta = self.ucqy/self.plqy
        return
    
    def _unpack_simulation(self):
        self.S = self.sol[-2, 0]
        self.T = self.sol[0, 0]
        return
    
    def _set_generation_rates(self):
        self._check_initial_weighting()
        self._set_initial_condition()
        self._GS = self._y0[self._initial_state_mapping['S']]
        self._GT = self._y0[self._initial_state_mapping['T']]
        return
    
    def _wrap_simulation_results(self):
        self.simulation_results = dict(zip(self.states, [self.S, self.T]))
        return
    
class Model1(_RateModel):
    r"""
    A class for performing calculations of the spin-statistical factor.
    
    
    Attributes
    ----------
    states : list of str
        The names of the excited state species.
    rates : list of str
        The names of the different rate constants in the model.
    model_name : str
        The name of the model.
    initial_weighting : dict
        Dictionary of (str, float) pairs. Key is the state name (str) and value is its initial weight (float).
    gaplaw_A : float
        Energy gap law prefactor in unts of per ns.
    gaplaw_gamma : float
        Material dependent gamma factor from the energy gap law.
    gaplaw_Evib : float
        Vibrational energy to use in the energy gap law in units of eV.
    gaplaw_kT : float
        Thermal energy for internal conversion calculations in units of eV.
    G : float
        The exciton generation rate for :attr:`initial_state`. Units of per volume per time.
    kD : float
        Rate of triplet-pair dissociation into free triplets. Units of per ns.
    kTTA : float
        (Arbitrary) effective linear rate of triplet-triplet annihilation . UNits of per ns.
    kIC1 : float
        Rate of internal conversion from :math:`^3(TT)` to :math:`T_1`. Units of per ns.
    kIC2 : float
        Rate of internal conversion from :math:`^3(TT)` to :math:`T_2`. Units of per ns.
    kIC21 : float
        Rate of internal conversion from :math:`T_2` to :math:`T_1`. Units of per ns.
    kRISC : float
        Rate of high-level intersystem crossing from :math:`T_2` to :math:`S_1`. Units of per ns.
    kTF : float
        Rate of triplet-pair fusion to form :math:`S_1`. Units of per ns.
    kSF : float
        Rate of singlet fission from :math:`S_1` to triplet-pair states. Units of per ns.
    kS : float
        Radiative rate of the :math:`S_1` state. Units of per ns.
    kT : float
        Nonradiative rate of the :math:`T_1` state. Should be set to zero for true calculations of spin statistics. Units of per ns.
    cs : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^1(TT)`.
    ct_m : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^3(TT)_m`, where m is x, y or z.
    cq_m : numpy.ndarray
        1D array containing the overlap factors between the 9 :math:`(T..T)` states and :math:`^5(TT)_m`, where m is a, b, x, y or z.
    simulation_results : dict
        Produced by :meth:`simulate`. Keys are the excited-state names (str), values are the simulated populations (float).
    plqy : float
        Calculated photoluminescence quantum yield from :meth:`calculate_plqy`.
    ucqy : float
        Calculated upconversion quantum yield from :meth:`calculate_ucqy`.
    eta : float
        Calculated spin-statistical factor from :meth:`calculate_eta`.
    """
    
    def __init__(self):
        super().__init__()
        # metadata
        self.model_name = 'Model 1'
        self._number_of_states = 12
        self.states = ['S', 'T']
        self.rates = ['kD', 'kTTA', 'kIC1', 'kIC2', 'kIC21', 'kRISC', 'kT', 'kTF', 'kSF', 'kS']
        # rates between excited states
        self.kD = 0.1
        self.kTTA = 1
        self.kIC1 = 30
        self.kIC2 = 30
        self.kIC21 = 15
        self.kRISC = 1000
        self.kTF = 100
        self.kSF = 100
        # rates of decay
        self.kS = 0.0625
        self.kT = 0
        # initial stuff
        self._allowed_initial_states = {'S', 'T'}
        self._initial_state_mapping = {'S': -2, 'T': 0}
        self.initial_weighting = {'T': 1}
        
    def _generate_rate_equation_matrix(self):
        #
        self.ct = [self.ct_x, self.ct_y, self.ct_z]
        self.cq = [self.cq_a, self.cq_b, self.cq_x, self.cq_y, self.cq_z]
        #
        self._rem = np.zeros((self._number_of_states, self._number_of_states))
        # T1
        self._rem[0, 0] = -(2*self.kTTA+self.kT)
        self._rem[0, 11] = self.kIC21
        for i in range(9):
            self._rem[0, i+1] = 2*self.kD+self.kIC1*(self.ct_x[i]+self.ct_y[i]+self.ct_z[i])
        # TT_l
        for i in range(9):
            self._rem[i+1, i+1] = -(self.kD+self.kTF*self.cs[i]+(self.kIC1+self.kIC2)*(self.ct_x[i]+self.ct_y[i]+self.ct_z[i]))
            self._rem[i+1, 0] = (1/9)*self.kTTA
            self._rem[i+1, 10] = self.kSF*self.cs[i]
        # S1
        self._rem[10, 10] = -(self.kS+self.kSF*np.sum(self.cs))
        self._rem[10, 11] = self.kRISC
        for i in range(9):
            self._rem[10, i+1] = self.kTF*self.cs[i]
        # T2
        self._rem[11, 11] = -(self.kIC21+self.kRISC)
        for i in range(9):
            self._rem[11, i+1] = self.kIC2*(self.ct_x[i]+self.ct_y[i]+self.ct_z[i])
        return
    
    def _generate_generation_rate_matrix(self):
        self._grm = np.zeros((self._number_of_states, 1))
        self._grm[-2, 0] = -self._GS
        self._grm[0, 0] = -self._GT
        return
            
    def simulate(self):
        """
        Calculate the excited state populations.

        Returns
        -------
        None.

        """
        self._set_generation_rates()
        self._generate_generation_rate_matrix()
        self._generate_rate_equation_matrix()
        self.sol = np.linalg.solve(self._rem, self._grm)
        self._unpack_simulation()
        self._wrap_simulation_results()
        return
    
    def calculate_plqy(self):
        """
        Calculate the photoluminscence quantum yield.

        Returns
        -------
        None.

        """
        self.initial_weighting = {'S': 1}
        self.simulate()
        self.plqy = self.kS*self.S/self.G
        return
        
    def calculate_ucqy(self):
        """
        Calculate the upconversion quantum yield.

        Returns
        -------
        None.

        """
        self.initial_weighting = {'T': 1}
        self.simulate()
        self.ucqy = 2*self.kS*self.S/self.G
        return
    
    def calculate_eta(self):
        """
        Calculate the spin-statistical factor.

        Returns
        -------
        None.

        """
        self.calculate_plqy()
        self.calculate_ucqy()
        self.eta = self.ucqy/self.plqy
        return
    
    def _unpack_simulation(self):
        self.S = self.sol[-2, 0]
        self.T = self.sol[0, 0]
        return
    
    def _set_generation_rates(self):
        self._check_initial_weighting()
        self._set_initial_condition()
        self._GS = self._y0[self._initial_state_mapping['S']]
        self._GT = self._y0[self._initial_state_mapping['T']]
        return
    
    def _wrap_simulation_results(self):
        self.simulation_results = dict(zip(self.states, [self.S, self.T]))
        return
