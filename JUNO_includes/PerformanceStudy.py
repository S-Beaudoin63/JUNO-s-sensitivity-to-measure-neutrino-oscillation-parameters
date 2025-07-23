import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

from JUNO_includes.FluxEstimation import FluxEstimation
from JUNO_includes.ResolutionEstimation import ResolutionEstimation
import JUNO_includes.Utils as JI

def nb_points_coarse(x):
    """
    Empirical way to get the number of points for the coarse estimation
    """
    return 3 if x<=3 else int(20*(1-np.exp(-x*np.log(2)/20)))


class PerformanceStudy:
    """
    This class automates the realization of different types of studies:
        - full_resolution_computation: estimates expected flux with detector effects,
            then realizes a coarse resolution estimation 
            into a precise resolution estimation,
            all of this for one parameter and a given time.

        - loop_time: runs full_resolution_computation for several data-taking times 
            to get the time evolution of the relative precision of one parameter.
        
        - loop_param: runs loop_time for all parameters of the given hierarchy
            to get the time evolution of the relative precision for all parameters.
        
        - multiparam_study: automates the call to Res_estim_multi
            by computing first the expected flux with detector effects. 
        
        - analyze_hierarchy_distinction: estimates the minimum data-taking time required for JUNO to be able to
            distinguish between both hierarchies with at least 3 sigmas.
    """


    def __init__(self, NH, IH, Cores, Isotopes):  
        # For easy access to PDG values
        self.data = [NH,IH] 
        self.Cores = Cores
        self.Isotopes = Isotopes


    def full_resolution_computation(self, param, hiera, time, nb_points=20,
                                 needplot=True, GUI=False):
        """
        Estimates expected flux with detector effects, 
	    then realizes a coarse resolution estimation,
        into a precise resolution estimation, 
	    all of this for one parameter and a given time.
	
	    Regarding nb_points, it's the number of points used in the precise estimation.
	    The one used for the coarse estimation is deduced from it using f(x) = 20(1-exp(-nb_points*ln(2)/20)
	    so that, if nb_points = 20: coarse study of 10 points into precise study of 20 points. 
        """
        if JI.incompability_params(param,hiera): # check if parameter is compatible with hierarchy
            return
        
        # Expected flux
        FluxStudy = FluxEstimation(self.data[hiera] if hiera==0 else None, self.data[hiera] if hiera==1 else None, self.Cores, self.Isotopes,time)  
        FluxStudy.plot_flux(needplot=False)
        FluxStudy.plot_visible_spectrum(needplot=False) 

        Visible_flux_expected = (FluxStudy.S_vis_NH_L if hiera==0 
                                 else FluxStudy.S_vis_IH_L)  # Selects flux for correct hierarchy

        # Created instance of ResolutionEstimation 
        Full_ResolutionEstimation = ResolutionEstimation(time,
                                    np.real(Visible_flux_expected).copy(),
                                    self.data,
                                    self.Cores,
                                    self.Isotopes)

        # Coarse estimation 
        Epsilon_coarse = Full_ResolutionEstimation.Res_estim(param, 
                                            hiera, # 0 or 1 for NH or IH 
                                            nb_points_coarse(nb_points),  
                                            coarse=True,
                                            needplot=needplot,
                                            GUI=GUI
                                            )

        # Precise estimation 
        res = Full_ResolutionEstimation.Res_estim(param,
                            hiera,  # NH
                            nb_points,#20, # nb points (more precise estimation)
                            coarse=False,
                            ep_mm = Epsilon_coarse[0],
                            ep_mp = Epsilon_coarse[1],
                            ep_pm = Epsilon_coarse[3],
                            ep_pp = Epsilon_coarse[4],
                            needplot=needplot,
                            GUI=GUI,
                            )
        
        if not GUI:
            return res
        else:
            return (res[0], res[1], res[2], res[3], res[4], res[5], Epsilon_coarse[-1], res[-1])

    
    def Graph_time(self, param, hiera, TIME, Precisions, nb_points,
                   GUI=False):    
        """
        Plots the results (Precisions) from the loop_time method
        """
        official_times = [100, 6*365.25, 20*365.25] # 100 days, 6 years, 20 years
        trueValue, articleValues, ep_pdg = JI.trueValue(self.data, param, hiera)
        
        lower_precisions = np.abs(np.array(Precisions).T)[2]/trueValue*100 # lower bound of relative precision
        upper_precisions = (np.array(Precisions).T)[5]/trueValue*100       # upper bound of relative precision

        fig = plt.figure()
        ax = plt.gca()

        ax.scatter(TIME, lower_precisions, marker = '^', color = 'blue', alpha = 0.7, label = "Lower bound")
        ax.scatter(TIME, upper_precisions, marker = 'v', color = 'red',  alpha = 0.7, label = "Upper bound")

        ax.plot(TIME, lower_precisions, color="blue", alpha=0.7, linewidth=1.5)
        ax.plot(TIME, upper_precisions, color="red",  alpha=0.7, linewidth=1.5)

        if articleValues!=None:
            article_precisions = np.array(list(articleValues.values()))/trueValue*100
            ax.scatter(official_times, article_precisions, color = 'black', alpha = 0.7, marker = 'x', label = "[2] estimation")
            ax.plot(official_times,article_precisions,color="black", linewidth=1.5, linestyle='dotted',alpha=0.5)

        ax.axhline(y = 1, color = 'gray', linestyle = 'dotted') 
        ax.set_xscale('log')
        ax.set_yscale('log')
        txthiera = "NH" if hiera==0 else "IH"
        ax.set_title(f"Time evolution of {param} relative precision for the {txthiera}", fontsize=14)
        ax.set_xlabel('JUNO data-taking times [days]', fontsize=12)
        ax.set_ylabel('Relative precision [%]', fontsize=12)
        ax.legend(fontsize=10)

        if not GUI:
            plt.show()
        else:
            return fig


    def loop_time(self, param, hiera, TIME, nb_points,
                  needplot=True, GUI=False):
        """
        Runs full_resolution_computation for several data-taking times 
        to get the time evolution of the relative precision of one parameter.
        """
        if JI.incompability_params(param,hiera): # check if parameter is compatible with hierarchy
            return


        print(f"\nBegining loop for parameter {param}")
        Precisions = []
        for t in TIME :
            print(f"Currently computing for {t} days")
            res = self.full_resolution_computation(param, hiera, t, nb_points,
                                                   needplot=False)
            Precisions.append(res) 

        print(Precisions)
        if needplot or GUI:
            fig = self.Graph_time(param, hiera, TIME, Precisions, nb_points, GUI=GUI)
        if not needplot:
            plt.close("all")

        if not GUI:
            return Precisions
        else:
            return fig


    def Graph_params(self, Lparam, hiera, TIME, Precisions, nb_points,
                     GUI=False):
        """
        Plots the results (Precisions) from the loop_params method
        """    
        official_times = [100, 6*365.25, 20*365.25] # 100 days, 6 years, 20 years
        colors = {"s12":"red" , "s13":"blue" , "m21":"green" , "m31_NH":"black" , "m13_IH":"black"}

        fig = plt.figure()
        ax = plt.gca()
       
        for param in Lparam:
            trueValue, articleValues, ep_pdg = JI.trueValue(self.data, param, hiera)

            Precs2 = np.array(Precisions[param]).T
            lower_precisions = np.abs(Precs2)[2]/trueValue*100  # lower bound 
            upper_precisions = (Precs2)[5]/trueValue*100        # upper bound
            
            if articleValues!=None:                                                                                         
                article_precisions = np.array(list(articleValues.values()))/trueValue*100                                                 
                ax.plot(official_times, article_precisions, color=colors[param], marker="x", linewidth=1.5, linestyle='dotted', alpha=0.7) 
            
            ax.plot(TIME, np.mean([lower_precisions,upper_precisions], axis=0), color=colors[param], linewidth=1.5, alpha=0.7, label=param, marker=".")

        ax.axhline(y = 1, color = 'gray', linestyle = 'dotted') 
        ax.set_xscale('log')
        ax.set_yscale('log')
        txthiera = "NH" if hiera==0 else "IH"
        ax.set_title(f"Time evolution of relative precisions for {txthiera}", fontsize=14)     
        ax.set_xlabel('JUNO data-taking time [days]', fontsize=12)
        ax.set_ylabel('Relative precision [%]', fontsize=12)
        ax.legend(fontsize=10)

        if not GUI:
            plt.show()
        else:
            return fig


    def loop_params(self, Lparam, hiera, TIME, nb_points,
                    needplot=True, GUI=False):
        """
        Runs loop_time for all parameters of the given hierarchy
        to get the time evolution of the relative precision for all parameters.
        """

        Lparam = Lparam[hiera]
        Precisions = {}
        for param in Lparam:
            res = self.loop_time(param, hiera, TIME, nb_points, needplot=False)
            Precisions[param] = res
        
        print("\n",Precisions,"\n")
        if needplot or GUI:
            fig = self.Graph_params(Lparam, hiera, TIME, Precisions, nb_points, GUI=GUI)
        if not needplot:
            plt.close("all")

        if not GUI:
            return Precisions
        else:
            return fig


    def multiparam_study(self, param1, param2, hiera, time,
                         nbpoints=10, needplot=True, GUI=False):
        """
        Automates the call to Res_estim_multi
        by computing first the expected flux with detector effects. 
        """
        if JI.incompability_params(param1,hiera) or JI.incompability_params(param2,hiera): # check if parameter is compatible with hierarchy
            return


        # Expected flux
        FluxStudy = FluxEstimation(self.data[0], self.data[1], self.Cores, self.Isotopes, time)
        FluxStudy.plot_flux(needplot=False)
        FluxStudy.plot_visible_spectrum(needplot=False)

        Visible_flux_expected = (FluxStudy.S_vis_NH_L if hiera==0 
                                        else FluxStudy.S_vis_IH_L)

        # Create instance of ResolutionEstimation
        Resolution_multi = ResolutionEstimation(time,
                                                np.real(Visible_flux_expected).copy(),
                                                self.data,
                                                self.Cores,
                                                self.Isotopes)

        # Does the Multiparametric resolution estimation
        res = Resolution_multi.Res_estim_multi(hiera,
                                               nbpoints,
                                               param1, param2,
                                               needplot=needplot, GUI=GUI)
        return res
    

    def analyze_hierarchy_distinction(self, time_min, time_max, nb_points, 
				      series=5, GUI=False):
        """
        Estimates the minimum data-taking time required for JUNO to be able to
        distinguish between both hierarchies with at least 3 sigmas.

        The method works as follows:
            - we assume NH to be True (although the same reasoning with IH yields the same results) 
              and we compute all the associated binned spectra. 
            - then we compute all IH spectra and add a gaussian noise in order to mimick statistical fluctuations
            - we compare NH and IH noisy through a Chi² method
            - we get the p-value and associated significance thresholds using Wilks' theorem  
            - we loop over time to obtain the significance evolution
            - we repeat the process 'series' times to average over all simulations and get a more stable result        
        """

        time = np.linspace(time_min, time_max, nb_points)
        
        # Creates empty array to be filled with significance time evolution 
        sigma_means = np.zeros_like(time)
        bins=510
        E_vis = np.linspace(1.5, 11, bins*10)

        for _ in range(series): # does series loops
            print(f"\nSeries n° {_+1} out of {series}")
            chi2_values = []
            for k,t in enumerate(time):
                print(f"{k+1} out of {len(time)}")
                
                # Computes fluxes for NH and IH
                FluxStudy = FluxEstimation(self.data[0], self.data[1], self.Cores, self.Isotopes, t, nb_points=bins*10)
                FluxStudy.plot_flux(needplot=False)
                FluxStudy.plot_visible_spectrum(needplot=False) # Detector effects 

                Flux_NH = np.real(FluxStudy.S_vis_NH_L).copy()
                Flux_IH = np.real(FluxStudy.S_vis_IH_L).copy()

                # Turns fluxes into histograms 
                histo_NH = np.histogram(E_vis, bins=bins, weights=Flux_NH/bins)[0]
                histo_IH = np.histogram(E_vis, bins=bins, weights=Flux_IH/bins)[0]

                # Adds gaussian noise to each bin of IH  
                noise = np.random.normal(0, np.sqrt(histo_IH), size = np.shape(histo_IH))
                histo_IH_noisy = histo_IH + noise

                # Computes Chi² 
                chi2 = np.sum((histo_IH_noisy - histo_NH) ** 2 / histo_NH)
                chi2_values.append(chi2)

            p_values = stats.chi2.sf(chi2_values, bins)
            sigma_values = stats.norm.isf(p_values)
            sigma_means += sigma_values
        sigma_means /= series

        # Finds the first times at which significance >= 3 or 4 (nb of sigmas)
        t_3sigma = [t for t, chi2 in zip(time, sigma_means) if chi2 >= 3]
        t_4sigma = [t for t, chi2 in zip(time, sigma_means) if chi2 >= 4]

        if len(t_3sigma)!=0:
            t_3sigma = t_3sigma[0]
        else:
            t_3sigma = -1
        if len(t_4sigma)!=0:
            t_4sigma = t_4sigma[0]
        else:
            t_4sigma = -1
        
        if t_3sigma < 100:
            t_3sigma = None
        if t_4sigma < 100:
            t_4sigma = None

        print(f"Data-taking time to distinguish hierarchies at 3 sigma: {t_3sigma} ; at 4 sigma: {t_4sigma}")

        fig1 = plt.figure()
        ax1 = plt.gca()
        ax1.plot(histo_IH_noisy, label="IH + gaussian noise", color = 'blue', alpha = 0.3, linewidth=1.5)
        ax1.plot(histo_IH, label = "IH", color = 'blue', alpha = 0.7, linewidth=1.5)
        ax1.plot(histo_NH, label="NH (expected)", color = 'red', alpha = 0.7, linewidth=1.5)
        ax1.set_title("Visualisation of NH and IH+noise spectra", fontsize=14)
        ax1.set_xlabel("$E_{vis}$ [MeV]", fontsize=12)
        ax1.set_ylabel("Visible flux",    fontsize=12)
        ax1.legend(fontsize=10)
        if not GUI:
            plt.show()

        fig2 = plt.figure()
        ax2 = plt.gca()
        ax2.plot(time, sigma_means, color = 'black', alpha = 0.7, linewidth=1.5, marker='.')
        ax2.axhline(3, linewidth=1.5, linestyle='dashed', color='yellow', label=r"$3\sigma$")
        ax2.axhline(4, linewidth=1.5, linestyle='dashed', color='orange', label=r"$4\sigma$")
        ax2.axhline(5, linewidth=1.5, linestyle='dashed', color='red', label=r"$5\sigma$")
        ax2.axvline(3*365.25, linewidth=1.5, linestyle='dashed', color='black', alpha=0.5)
        ax2.axvline(6*365.25, linewidth=1.5, linestyle='dashed', color='black', alpha=0.5)
        ax2.fill_between(time, 3, 4, color='yellow', alpha=0.1)
        ax2.fill_between(time, 4, 5, color='orange', alpha=0.1)
        ax2.fill_between(time, 5, 6, color='red', alpha=0.1)
        ax2.set_title("Time evolution of the significance of hierarchy distinction", fontsize=14)
        ax2.set_xlabel("Data-taking time (days)", fontsize=12)
        ax2.set_ylabel(r"Significance ($\sigma$)", fontsize=12)
        ax2.set_xlim(time_min, time_max)
        ax2.set_ylim(0, 6)
        ax2.set_xscale('log')
        ax2.legend(fontsize=10)

        if not GUI:
            plt.show()
            return t_3sigma, t_4sigma
        else:
            return  t_3sigma, t_4sigma, fig1, fig2

