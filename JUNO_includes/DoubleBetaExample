import numpy as np
from matplotlib import pyplot as plt

"""
This class is used to provide an example of application of the hierarchy 
parameter constraints and the difference that hierarchies make.

Studying neutrinoless double beta decay allows to probe the nature of neutrinos (Dirac or Majorana fermions).

The parameter of interest is the effective Majorana mass which is a combination of the neutrino
mass eigenstates and the neutrino mixing matrix terms.

The effective mass can be expressed as a function of the lightest neutrino mass and relates to the inverse
of the processe's lifetime.

This code generates the accessible parameter space depending on the hierarchy, providing an example of 
consequence from Juno's hierarchy measurements.
"""


class DoubleBetaExample:
    def __init__(self):
        self.mini_mass = np.linspace(1e-4, 1, num = 2048*6)


    def mass_bb(self, alpha, beta) :
        """
        The PMNS matrix is supplemented by two Majorana phases alpha and beta.
        This function explicitly computes the effective mass for a pair of phases.
        """
        m21 = 7.53e-5
        m32_NH = 2.455e-3
        m31_IH = 2.454e-3
        s_12 = 3.07e-1
        s_13NH = 2.19e-2
        s_13IH = 2.19e-2
        delta = 1.19

        mbb_NH = []
        mbb_IH = []
        
        for mass in self.mini_mass :
            # This loop parametrizes the absolute neutrino masses and computes the effective mass
            m1_NH = mass
            m2_NH = np.sqrt(mass**2 + m21)
            m3_NH = np.sqrt(mass**2 + m32_NH + m21/2)
            
            mbb_auxNH = np.abs((1-s_12)*(1-s_13NH)*np.exp(2*1j*alpha)*m1_NH 
                                + (1-s_13NH)*s_12*np.exp(2*1j*beta)*m2_NH
                                + s_13NH*m3_NH*np.exp(-2*1j*delta))
            
            m1_IH = np.sqrt(mass**2 + m31_IH - m21/2)
            m2_IH = np.sqrt(mass**2 + m31_IH + m21/2)
            m3_IH = mass
        
            mbb_auxIH = np.abs((1-s_12)*(1-s_13IH)*np.exp(2*1j*alpha)*m1_IH 
                                + (1-s_13IH)*s_12*np.exp(2*1j*beta)*m2_IH
                                + s_13IH*m3_IH*np.exp(-2*1j*delta))
        
            mbb_NH.append(mbb_auxNH)
            mbb_IH.append(mbb_auxIH)

        return mbb_NH, mbb_IH


    def DoubleBetaExample_graph(self, GUI=False):
        """
        Plots the graph for several phases pairs 
        """
        fig, ax = plt.subplots()
        alpha = np.linspace(0, np.pi/2, num = 15) 
        beta = np.linspace(0, np.pi/2, num = 15)

        first_loop = True
        for a in alpha :
            for b in beta :
                # loops over the phases
                NH_00, IH_00 = self.mass_bb(a, b)
                ax.plot(self.mini_mass, NH_00, color = 'red', alpha = 0.3, linewidth=1.5, label="NH" if first_loop else "")
                ax.plot(self.mini_mass, IH_00, color = 'blue', alpha = 0.3, linewidth=1.5, label="IH" if first_loop else "")
                first_loop = False

        ax.set_xlim(10**(-4), 1)
        ax.set_ylim(10**(-4), 1)
        ax.set_ylabel("$|m_{ββ}|$ [eV]", fontsize=12)
        ax.set_xlabel("$m_{min}$ [eV]", fontsize=12)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title("Invariant mass of double beta", fontsize=14)
        ax.legend(fontsize=10)
        if not GUI: 
            plt.show()
        else:
            return fig
