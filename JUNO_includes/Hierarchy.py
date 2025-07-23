import numpy as np
import pandas as pd

class hierarchy:
    """
    This class is used to contain the PDG values for the oscillation parameters 
    but is also used to make computate other properties of the two different hierarchies (such as the PMNS matrix)
    and will later be used to create "custom" hierarchies with parameters shifted from their PDG values to 
    study the influence of the parameters values on observable antineutrino phenomena.
    """

    def __init__(self, hiera, X):
        """
        Call the right constructor, because this class has of two "constructors" that are called 
	depending on the context and their arguments. 

	Unlike in C++, python does not allow to surcharge methods using different arugments, so it is necessary to have a "real"
	constructor hierarchy.__init__(hiera,X) that receives hiera (0 or 1 depending on if the hierarchy is inverted or not) and X
	which can be either a string (path to hierarchies.csv) or a list.
	Depending on the type of X, __init__ will call either one of the two constructors.
        """
        if len(X)==1:
            self.init_csv(X[0],hiera)
        elif len(X[1]) == 3:
            self.init_manu(X[0],X[1],X[2],X[3])



    def init_csv(self, file, hiera):
        """
        Constructor using pdg2024 values, stored in hierarchies.csv
        Takes as argument the path of the CSV
        The hierarchy object has the same attributes than if it was created with the manual constructor
        """
        dh = pd.read_csv(file)

        self.inverted = bool(dh["inverted"][hiera])                 # Normal or inverted
        self.delta_radians = np.radians(dh["delta_deg"][hiera]*180) # CP phase (conversion because it is stored in degrees)

        # sin & cos
        self.s_12 = dh["s12"][hiera]
        self.s_13 = dh["s13"][hiera]
        self.s_23 = dh["s23"][hiera]
        self.c_12 = 1 - self.s_12
        self.c_13 = 1 - self.s_13
        self.c_23 = 1 - self.s_23

        # squared mass differences
        self.m_21 = dh["m21"][hiera]
        if not self.inverted:       # Normal hierarchy
            self.m_32_NH = dh["m32_NH"][hiera]
            self.m_31_NH = self.m_21 + self.m_32_NH
        else:                       # Inverted hierarchy
            self.m_23_IH = dh["m23_IH"][hiera]
            self.m_13_IH = self.m_23_IH - self.m_21

        self.PMNS_matrix() # constructs the PMNS matrix

        

    def init_manu(self, inverted, Dmasses2, sins, delta_radians):
        """
        Constructor using custom values, when shifting the parameter value by some amount (in ResolutionEstimation.Res_estim)
        Takes as argument the characteristics of the hierachy object to be created:
        - inverted: 0 or 1 if the hierachy is normal or inverted 
        - Dmasses2: list of the 3 defining masses (m12, m23 and m31_NH or m13_IH)
        - sins: list of the 3 sines (s12,s23,s13)
        - delta_radians: the value in radians of the CP phase
        """
        if len(Dmasses2)!=3 or len(sins)!=3:
            print("Error")
            return(-1)

        self.inverted = inverted   # Normal or inverted
        self.delta_radians = delta_radians

        # sin & cos
        [self.s_12, self.s_13, self.s_23] = sins
        self.c_12 = 1 - self.s_12
        self.c_13 = 1 - self.s_13
        self.c_23 = 1 - self.s_23

        # squared mass differences
        self.m_21 = Dmasses2[0]
        if not self.inverted:      # normal hierarchy
            self.m_32_NH = Dmasses2[1]
            self.m_31_NH = Dmasses2[2]
        else:                       # Inverted hierarchy
            self.m_23_IH = Dmasses2[1]
            self.m_13_IH = Dmasses2[2]

        self.PMNS_matrix() # constructs the PMNS matrix



    def PMNS_matrix(self):
        """
        Generates the PMNS matrix modelling the neutrinos' mixing and oscillations 
        """
        self.PMNS = np.array([
                [np.sqrt(self.c_12*self.c_13), np.sqrt(self.s_12*self.c_13), np.sqrt(self.s_13)*np.exp(-1j*self.delta_radians)],
                [-np.sqrt(self.s_12*self.c_23) - np.sqrt(self.c_12*self.s_13*self.s_23)*np.exp(1j*self.delta_radians), np.sqrt(self.c_12*self.c_23) - np.sqrt(self.s_12*self.s_13*self.s_23)*np.exp(1j*self.delta_radians), np.sqrt(self.c_13*self.s_23)],
                [np.sqrt(self.s_12*self.s_23) - np.sqrt(self.c_12*self.s_13*self.c_23)*np.exp(1j*self.delta_radians), -np.sqrt(self.c_12*self.s_23) - np.sqrt(self.s_12*self.s_13*self.c_23)*np.exp(1j*self.delta_radians), np.sqrt(self.c_13*self.c_23)]
            ])


    def oscillation_anti(self, energy, baseline, initial_state, final_state):
        """
        Computes the oscillation probability from some initial antineutrino flavor to a final one, 
        for a given energy, baseline and hierarchy.
        This is a general method allowing to consider all flavor pairs. 

        We will consider only electron antineutrinos for the study of Juno 
        but doing so allows to obtain the oscillation patterns for any combination of neutrinos. 
        """
        U = np.conj(self.PMNS)

        if not self.inverted:   # NH
            self.mass = [[0, -self.m_21, -self.m_31_NH], [self.m_21, 0, -self.m_32_NH], [self.m_31_NH, self.m_32_NH, 0]]
        else:                   # IH
            self.mass = [[0, -self.m_21, self.m_13_IH], [self.m_21, 0, self.m_23_IH], [-self.m_13_IH, -self.m_23_IH, 0]]

        self.P_list = [] # List of oscillation probabilities 

        # Denominator of the probability (normalisation factor)
        den = (U @ np.conj(U.T))[initial_state][initial_state]*(U @ np.conj(U.T))[final_state][final_state] 
        
        for E in energy:  
            E = E*10**(-3)   # Converts E from MeV to GeV for a properly dimensioned computation
            P = 0
            for i in range(0, 3): 
                P += np.abs(U[initial_state][i])**2*np.abs(U[final_state][i])**2
                for j in range(i+1, 3):
                    P += 2*np.real(U[initial_state][i]*U[final_state][j]*np.conj(U[initial_state][j])*np.conj(U[final_state][i]))*np.cos(self.mass[i][j]*(baseline/E)*2.54) \
                        - 2*np.imag(U[initial_state][i]*U[final_state][j]*np.conj(U[initial_state][j])*np.conj(U[final_state][i]))*np.sin(self.mass[i][j]*(baseline/E)*2.54)
            self.P_list.append(P/den)
        return self.P_list





