import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.stats import norm


class FluxEstimation:
    """
    This class models the antineutrino flux to be expected in the detector, coming from the nuclear plants
    There are several methods, each taking into account specific aspects.

    S_iso:                      characteristic flux for a given isotope
    IBD_cross_section:          contribution from inverse beta decay
    flux_tot:                   ideal flux (without detector effects)
    compute_flux_integrals:     ideal number of antineutrinos events per day
    energy_resolution:          energy resolution of each photomultiplier tube
    detector_response:          maps antineutrinos energies to visible detector energies
    simulate_visible_spectrum:  applies all effects to get a realistic antineutrino flux
    """
    
    def __init__(self, NH, IH, Cores, Isotopes, time_step_days, nb_points=4096):
        """
        Loads the data needed into the class.
        The antineutrino energies are at least equal to 1.8 MeV due to the beta decay threshold.  
        """
        self.Energies = np.linspace(1.805, 10, nb_points) # Emin= 1.8MeV so we only plot for E>1.8MeV
        self.NH = NH
        self.IH = IH
        self.Cores = Cores
        self.Isotopes = Isotopes
        self.time_step_days = time_step_days

        # Detector efficiency and number of protons per target
        self.epsilon = 0.822
        self.N_p = 1.44e33

        # Juno uses two sizes of photomultiplier tubes (Large and Small).
        # Each one with its own energy resolution parametrized as in [8]
        self.Resol_params = {"LPMT":np.array([0.0261, 0.0082, 0.0123]),
                             "SPMT":np.array([0.1536, 0.0082, 0.0677])
                             }
        self.E_vis = np.linspace(1.5, 11, nb_points) # Energy visible in the detector 


    # Isotopes contributions

    def S_iso(self, isotope):
        """
        Each reactor core emits antineutrinos through the beta decay of fission products 
        from the elements given in isotopes.csv. 
        
        The characteristic spectra (depending on the released energy fraction, etc...) follows the 
        parametrization of [3],[4]

        This function returns the characteristic antineutrino spectra of a given isotope.
        """
        Nb_nu_params = self.Isotopes[isotope].FluxesPerFission 

        res = 0
        for param_poly in range(0, 6) :
            res += Nb_nu_params[param_poly]*self.Energies**param_poly
        return np.exp(res)


    def plot_Nb_nu_pfission(self,GUI=False):
        """
        Complementary method to S_iso that plots the number of neutrinos per MeV
         per fission for all considered isotopes 
        """
        fig = plt.figure()
        ax=plt.gca()
        ax.set_xlabel(r'$E_{\nu}$ [MeV]', fontsize=12)
        ax.set_ylabel(r'$N_{\nu}$ / fission', fontsize=12)
        ax.set_title("Number of antineutrinos per fission for different isotopes", fontsize=14)
        for isotope in self.Isotopes.keys() :
            ax.plot(self.Energies, self.S_iso(isotope),linewidth= 1.5)
        ax.legend(self.Isotopes.keys(), fontsize=10)
        
        if not GUI:
            plt.show(block=True)
        else:
            return(fig)

    # Modelling ideal antineutrino flux  

    def IBD_cross_section(self):  # IBD = Inverse beta decay
        """
        Returns the Inverse Beta Decay cross section (approximation from [5])
        """
        m_e = 0.511   # Electron mass in MeV
        E_e = self.Energies - 939.5654 + 938.2721   # Electron's energy = E(antineutrino) - m(neutron) + m(proton)
        p_e = np.sqrt(E_e**2 - m_e**2)  # Electron's momentum
        return 0.0952*E_e*p_e*1e-52  # cross section converted in km²


    def flux_tot(self):
        """
        Computes the total ideal flux expected in the detector (without detector effects)
        assuming all reactor contributions, 
        for a given data taking time
            - without antineutrino oscillations
            - with both hierarchies oscillations
    
        As the expected spectrum results from the detected IBD events, the reactor flux has to be convoluted
            with the IBD cross section, the number of free protons in the detection target, 
            the electron antineutrino survival probability from the PMNS matrix [6] 
            and the detector efficiency (82.2%)
        """
        self.S_tot = 0          # Antineutrino flux without oscillations
        self.Oscil_tot_NH = 0   # Antineutrino flux with normal ordering oscillations
        self.Oscil_tot_IH = 0   # Antineutrino flux with inverted ordering oscillations

        for core in self.Cores.keys(): # Takes into account all cores
            S_aux = 0 
            den = 0
            aux = 0
            for iso in self.Isotopes.keys():
                den += self.Isotopes[iso].FuelFrac*self.Isotopes[iso].Ereleased
                aux += self.Isotopes[iso].FuelFrac*self.S_iso(iso)

            # Individual core flux
            S_aux = self.Cores[core].PowerMeVd/(4*np.pi*self.Cores[core].Baseline**2)*aux/den*self.epsilon*self.N_p*self.IBD_cross_section()*self.time_step_days

            self.S_tot += S_aux   # Total flux

            if self.NH != None:
                self.Oscil_tot_NH += S_aux*self.NH.oscillation_anti(self.Energies, self.Cores[core].Baseline, 0, 0)
            if self.IH != None:
                self.Oscil_tot_IH += S_aux*self.IH.oscillation_anti(self.Energies, self.Cores[core].Baseline, 0, 0)


    def compute_flux_integrals(self):
        """
        Compute the number of antineutrino to be expected per day
        by integrating the flux using Simpson's method.
        """
        self.integral_no_osc = np.real(simpson(self.S_tot, self.Energies)/self.time_step_days)
        self.integral_NH = np.real(simpson(self.Oscil_tot_NH, self.Energies)/self.time_step_days)
        self.integral_IH = np.real(simpson(self.Oscil_tot_IH, self.Energies)/self.time_step_days)
        

    def plot_flux(self, needplot=True, GUI=False):
        """
        Plots the graph of the ideal flux (with/without oscillations) 
        and give the number of antineutrinos expected per day

        needplot argument is True when we want to plot the result, otherwise we don't want to lose time
        GUI argument indicates that the method was called from the GUI and we want to return the figures to show them 
            in the GUI instead of plt.showing them
        """
        self.flux_tot()

        if needplot:
            self.compute_flux_integrals()
            fig = plt.figure()
            ax = plt.gca()

            ax.plot(self.Energies, self.S_tot, label='No oscillation (all cores)', linewidth=1.5, linestyle="dashed", color="black", alpha=0.7)

            #oscillations
            ax.plot(self.Energies, self.Oscil_tot_NH, label='NH',linewidth=1.5, color='red', alpha=0.7)
            ax.plot(self.Energies, self.Oscil_tot_IH, label='IH',linewidth=1.5, color='blue', alpha=0.7)

            ax.set_xlabel(r'$E_{\nu}$ [MeV]', fontsize=12)
            ax.set_ylabel('Number of events / MeV', fontsize=12)
            ax.text(0.009, 0.98, f'{np.round(self.time_step_days/365.25,2)} years of data', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
            ax.set_title('Comparison of neutrino fluxes with/without oscillation', fontsize=14)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(3,3))

            # Adds the integrals to the plot
            text_x = self.Energies[-1]*0.995
            text_y = max(np.max(self.S_tot), np.max(self.Oscil_tot_NH), np.max(self.Oscil_tot_IH))*0.65
            ax.text(text_x, text_y, f"No Oscillation : {self.integral_no_osc:.4} / day \nNH : {self.integral_NH:.3} / day \nIH : {self.integral_IH:.3} / day", fontsize=10, ha='right', va='bottom', bbox=dict(facecolor="#f0f0f0" if needplot and GUI else 'white', alpha=0.8), multialignment='left')
        
            ax.legend(fontsize=10)
            
            if not GUI:
                plt.show(block=True)
            else:
                return fig,self.integral_no_osc,self.integral_NH,self.integral_IH


    # Modelling of detector effects 

    def energy_resolution(self, PMT_key):
        """
        Takes a photomultiplier (PM) into argument.
        The visible energy resolution is defined by a gaussian.
        It takes 3 contributions:
            - 1 driven by the Poisson statistic of the number of detected photoelectrons
            - 1 dominated by the detector's spatial non uniformity 
            - 1 dominated by dark noise
        """
        sigma_Eres = np.sqrt((self.Resol_params[PMT_key][0]/np.sqrt(self.E_vis))**2 + self.Resol_params[PMT_key][1]**2 + (self.Resol_params[PMT_key][2]/self.E_vis)**2)*self.E_vis
        return sigma_Eres

    def detector_response(self, PMT_key):
        """
        Apply the detector response function, which includes energy nonlinearity and resolution.
            Parameters:
            - self.Energies: True neutrino energies
            - self.E_vis: Measured visible energies
            Returns:
            - Response matrix R_map(self.Energies, self.E_vis) which associates each antineutrino energy 
              with a probability to be detected as a visible energy
        """
        R_map = np.zeros((len(self.Energies), len(self.E_vis)))
        for i, e in enumerate(self.Energies):
            sigma_Eres = self.energy_resolution(PMT_key)
            R_map[i,:] = norm.pdf(self.E_vis, loc = e, scale = sigma_Eres)
        return R_map   # The element R_map[i][j] is the probability density that a neutrino with energy E_ν[i] is measured as E_vis[j]


    def simulate_visible_spectrum(self, flux, PMT_key):
        """
        Simulate the flux and response of the detector.
        Returns the expected visible antineutrino spectrum observed at Juno.
        We only take into account energy resolution effects for the two PMT and neglect energy non linearities. 
        """
        R_map = self.detector_response(PMT_key)
        S_vis = (flux @ R_map) * (self.Energies[1] - self.Energies[0])  # Integrate flux with response
        return S_vis   # Visible flux


    def plot_visible_spectrum(self, needplot=True, GUI=False):
        """
        Plots the visible flux with all detector effects, with and without PMTs.
        Displays the real number of antineutrinos per day.

        The effect of the PMTs tends to attenuate fine structures in the oscillation spectrum.
        We only compute small PMTs if we want to compare on the plot but future functions will only
        consider the effects of large PMTs (_L), because small PMTs receive a smaller portion of the light.

        needplot argument is True when we want to plot the result, otherwise we don't want to lose time
        GUI argument indicates that the method was called from the GUI and we want to return the figures to show 
            them in the GUI instead of plt.showing them
        """
        if self.NH != None:
            self.S_vis_NH_L = self.simulate_visible_spectrum(self.Oscil_tot_NH, 'LPMT')
        if self.IH != None:
            self.S_vis_IH_L = self.simulate_visible_spectrum(self.Oscil_tot_IH, 'LPMT')

        if needplot:
            self.S_vis_NH_S = self.simulate_visible_spectrum(self.Oscil_tot_NH, 'SPMT')
            self.S_vis_IH_S = self.simulate_visible_spectrum(self.Oscil_tot_IH, 'SPMT')

            # compute the integrals for the real number of detected antineutrinos
            i_NH_L = np.real(simpson(self.S_vis_NH_L, self.E_vis))/self.time_step_days
            i_IH_L = np.real(simpson(self.S_vis_IH_L, self.E_vis))/self.time_step_days
            i_NH_S = np.real(simpson(self.S_vis_NH_S, self.E_vis))/self.time_step_days
            i_IH_S = np.real(simpson(self.S_vis_IH_L, self.E_vis))/self.time_step_days

            fig = plt.figure()
            ax = plt.gca()
            # 0.02 factor comes from the fact that we regroup by 20 keV for convenience.
            ax.plot(self.Energies, self.Oscil_tot_NH*0.02, label=r'Ideal NH', linewidth=1.5, linestyle="dashed", color='red', alpha=0.3)
            ax.plot(self.Energies, self.Oscil_tot_IH*0.02, label=r'Ideal IH', linewidth=1.5, linestyle="dashed", color='blue', alpha=0.3)

            ax.plot(self.E_vis, self.S_vis_NH_L*0.02, label=r'NH LPMT', linewidth=1.5, color='red', alpha = 0.7)
            ax.plot(self.E_vis, self.S_vis_IH_L*0.02, label=r'IH LPMT', linewidth=1.5, color='blue', alpha = 0.7)
            ax.plot(self.E_vis, self.S_vis_NH_S*0.02, label=r'NH SPMT', linewidth=1.5, color='black', alpha = 0.2)
            ax.plot(self.E_vis, self.S_vis_IH_S*0.02, label=r'IH SPMT', linewidth=1.5, color='black', alpha = 0.3)

            ax.set_xlabel('$E_{vis}$ [MeV]', fontsize=12)
            ax.set_ylabel('Events / 20 keV', fontsize=12)
            ax.set_title('Comparison between ideal and visible energy \n spectrum after detector response', fontsize=14)
            ax.text(0.009, 0.98, f'{np.round(self.time_step_days/365.25,2)} years of data', transform=plt.gca().transAxes, fontsize=10, va='top', ha='left')
            ax.legend(fontsize=10)

            # Adds the integrals to the plot
            text_x = self.E_vis[-1]*0.995 
            ax.text(text_x, 200, f"NH LPMT : {i_NH_L:.3} / day \nNH SPMT : {i_NH_S:.3} / day \nIH LPMT : {i_IH_L:.3} / day \nIH SPMT : {i_IH_S:.3} / day",
                    fontsize=10, ha='right', va='bottom', bbox=dict(facecolor="#f0f0f0" if needplot and GUI else 'white', alpha=0.8), multialignment='left')

            if not GUI:
                plt.show()
            else:
                return fig,i_NH_L,i_IH_L,i_NH_S,i_IH_S