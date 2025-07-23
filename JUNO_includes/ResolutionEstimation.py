import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import matplotlib.ticker as ticker

from JUNO_includes.Hierarchy import hierarchy
from JUNO_includes.FluxEstimation import FluxEstimation
import JUNO_includes.Utils as JI

def scientific_formatter(x, pos):
     """
     Used for the format of some matplotlib plots
     """
     return '{:.1e}'.format(x)

def double_gaussian(x, mu1, sigma1, mu2, sigma2):
     """
     Returns a double gaussian function for later convenience
     """
     return stats.norm.pdf(x, mu1, sigma1) + stats.norm.pdf(x, mu2, sigma2)


def adaptive_linspace(x_min, x_max, n_points, mu1, sigma1, mu2, sigma2):
     """
     Equivalent to the linspace function but following a double gaussian distribution.
     This is used to optimize the number of points used in the loop of function Res_estim down below,
     by focusing the points in the two regions of interest.
     """
     x = np.linspace(x_min, x_max, 1000)
     pdf = double_gaussian(x, mu1, sigma1, mu2, sigma2)
     cdf = np.cumsum(pdf)  
     cdf = (cdf - cdf.min()) / (cdf.max() - cdf.min())    # Normalize to [0, 1]
    
     uniform_samples = np.linspace(0, 1, n_points)        # Evenly spaced values in the CDF domain
     adaptive_points = np.interp(uniform_samples, cdf, x) # Interpolation to map uniform samples to x values
     return adaptive_points


def search_for_epsilon(EPSILON, DeltaCHI2):
     """
     This function is used in the Res_estim function further down.
     It is used to find a bounding of the points where the DeltaCHI2 is equal to 1.
     The DeltaCHI2 function is a parabola pointing upwards 
     so there are two points (epsilon- and epsilon+) to find. 
     """
     i=0
     while DeltaCHI2[i] > 1:
          i+=1
     ep_mm = EPSILON[i-1]      # lower bound for epsilon-
     imm = i-1
     ep_mp = EPSILON[i]        # upper bound for epsilon-
     imp=i

     while DeltaCHI2[i] < 1:
          i+=1
     ep_pm = EPSILON[i-1]     # lower bound for epsilon+
     ipm = i-1
     ep_pp = EPSILON[i]       # upper bound for epsilon+
     ipp=i
     return(ep_mm, imm, ep_mp, imp, ep_pm, ipm, ep_pp, ipp)  # returns the bounds and their ranks

 
def reconstruct_parabola(x1, y1, x2, y2, x3, y3):
     """
     A general function that reconstruct a parabola from 3 points
     """
     den = (x1-x2) * (x1-x3) * (x2-x3)
     A = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / den
     B = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / den
     C = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / den
     return A,B,C
 

class ResolutionEstimation:
     """
     This class is used to estimate JUNO's resolutions to the oscillation parameters of interest.
     The first method (Res_estim) assumes uncorrelated parameters:
          - each parameter is associated to a range so that we can make them vary by adding an epsilon value.
          - we compare the binned best known spectrum (defined with parameters from PDG 2024)
            to the spectrum with a slightly shifted parameter through a Chi² procedure. 
          - the plot_chi2_parabola method plots the graph associated to Res_estim 
            and compares results with PDG and article [2]
     
     The second method (Res_estim_multi) allows to verify the correlation factors between parameters:
          - we select a pair of parameters, the other two are fixed to their best known values 
          - we make the two parameters vary and draw a Chi² map allowing to represent the different significance regions
            according to Wilks' theorem
          - then we relate it to the correlation factor through the computation of the Hessian matrix
     """

     def __init__(self, time_step_days, S_expect, data, Cores, Isotopes):
          self.time_step_days = time_step_days    # Data-taking time
          self.E_vis = np.linspace(1.5, 11, 4096) # Visible energy range
          self.S_expected = S_expect              # Expected flux without any shift on the parameters (with LPMT)
          self.data = data                        # For easy access to PDG values
          self.Cores = Cores
          self.Isotopes = Isotopes




     def plot_chi2_parabola(self, param, hiera, param_shifted, DeltaCHI2, epsilon_bounds,
                            epsilon_indexes, epsilon_means, start, stop, GUI=False):
          """
          Plots the DeltaChi² = Chi² - min(Chi²) computed in Res_estimation,
          the paraboloas reconstructed using the uncertainties provided from PDG 2024
          and gives the expected outcome of [2]
          """
          
          fig,ax = plt.subplots()  
 
          # Retrieves the parameter value from PDG and 
          # the parameters uncertainties from [2] (given at 100 days, 6 years and 20 years) and from PDG
          trueValue, ep_article, ep_pdg = JI.trueValue(self.data, param, hiera) 

          # Reconstructs a DeltaChi² (parabolic) profile corresponding to current PDG uncertainties  
          [x1,y1] = [trueValue - ep_pdg , 1]
          [x2,y2] = [trueValue,0]
          [x3,y3] = [trueValue + ep_pdg, 1]
          # Calculates the unknowns of the equation y=ax^2+bx+c
          a, b, c = reconstruct_parabola(x1, y1, x2, y2, x3, y3)    

          # Plot the PDG profile, vertical lines indicates parameters shifted by one sigma 
          linspace_parabola = np.linspace(start,stop,100) + trueValue
          pdg_parabola = [a*((x)**2)+b*(x)+c for x in linspace_parabola]

          ax.plot([trueValue-ep_pdg,trueValue-ep_pdg], [0,10],  color="red", alpha=0.7, linewidth=1.5, linestyle='dotted')
          ax.plot([trueValue+ep_pdg,trueValue+ep_pdg], [0,10],  color="red", alpha=0.7, linewidth=1.5, linestyle='dotted')
          ax.plot(linspace_parabola,pdg_parabola,label = "PDG",  color="red", alpha=0.7, linewidth=1.5, linestyle='dashed')


          ax.plot(param_shifted, DeltaCHI2, label=f"JUNO {np.round(self.time_step_days/365.25,2)} yrs", color="blue", alpha=0.7, linewidth=1.5)

          # Wilks' theorem significance and corresponding DeltaChi²
          ax.axhline(1,color="black", linewidth=1.5, linestyle='dotted', alpha=0.7) # 3 sigma
          ax.axhline(4,color="black", linewidth=1.5, linestyle='dotted', alpha=0.7) # 4 sigma
          ax.axhline(9,color="black", linewidth=1.5, linestyle='dotted', alpha=0.7) # 5 sigma

          # Approximate location of the DeltaChi²=1 crossing points
          ax.scatter(np.array(epsilon_bounds)+trueValue,[DeltaCHI2[i] for i in epsilon_indexes],s=20,c="k")

          # Values from [2] represented as vertical lines
          if (param in ["m21","s12","s13","m31_NH"]) and (self.time_step_days in [100,6*365.25,20*365.25]):
               ep_article = ep_article[self.time_step_days]
               ax.plot([trueValue-ep_article, trueValue-ep_article], [0,10], color="green", alpha=0.7, linewidth=1.5, linestyle='dashed', label="[2] estimation")  
               ax.plot([trueValue+ep_article, trueValue+ep_article], [0,10], color="green", alpha=0.7, linewidth=1.5, linestyle='dashed')                   

          acq_time = np.round(self.time_step_days/365.25,2)
          acq_time = str(acq_time)+" years" if acq_time>1 else str(self.time_step_days)+" days"
          ax.set_title(r"$\Delta\chi^2$" + f" distribution for {param} after {acq_time}", fontsize=14)   
          ax.set_xlabel(param, fontsize=12)                                                              
          ax.set_ylabel(r"$\Delta\chi^2$", fontsize=12)
          ax.tick_params(axis='both', labelsize=10)
          ax.legend(fontsize=10)                                                                         
          ax.set_ylim(-1,10)

          ax.xaxis.set_major_formatter(ticker.FuncFormatter(scientific_formatter))

          # Adds the epsilon means to the plot
          ax.text(0.98, 0.02, r"$\epsilon_{-}$"+f" : {np.round(epsilon_means[0],8):.2e}\n"+r"$\epsilon_{+}$"+f" : {np.round(epsilon_means[1],8):.2e}",
                  transform=plt.gca().transAxes, ha='right', va='bottom', bbox=dict(facecolor="#f0f0f0" if GUI else 'white', alpha=0.8), multialignment='left', fontsize=10)    ############

          if not GUI:
               plt.show()
          else:
               return fig





     def Res_estim(self, param, hiera, nb_points,
                  coarse=False, ep_mm=None, ep_mp=None, ep_pm=None, ep_pp=None, needplot=True, GUI=False):
          """
          Estimates the relative precision for a non correlated parameter 'param' (s12, s13, m21, m31_NH or m13_IH) after a data-taking time self.time_step_days
          hiera = 0 or 1 is NH or IH
          bins = 510 is the number of energy bins assumed from previous Daya Bay information
          nb_points is the number of shifted parameter we evaluate

          This method aims at finding the one sigma uncertainty on a given parameter using a Chi² procedure 
          between an expected binned flux and a shifted parameter flux.  
          To do so it needs to explore different shifted parameter values.
        
	  It is called twice in (PerformanceStudy.full_resolution_computation):
               - a first time with coarse = True. It explores a small linspace to approximately locate the DeltaChi²=1 crossing points
               - then, ep_mm, ep_mp, ep_pm and ep_pp are filled with the coarse estimates 
                 and we reapply the method with a double gaussian addaptative linspace centered on the previously 
                 estimated ep points.
          This allows to use a small amount of points while assuring a good accuracy.
          """
          bins=510

          # The bounds of the linspace to explore depending on the parameter, empirically adjusted
          bounds = {"m21":    {"coarse":(-0.2e-5,0.3e-5),         "precise":(-0.2e-5,0.3e-5)},
                    "s12":    {"coarse":(-0.02,0.02),             "precise":(-0.027,0.027)},
                    "s13":    {"coarse":(-0.0219,0.03),           "precise":(-0.0219,0.03)},
                    "m31_NH": {"coarse":(-3.806e-5,3.806e-5),     "precise":(-3.806e-5,3.806e-5)},
                    "m13_IH": {"coarse":(-3.8944e-5,3.8944e-5),   "precise":(-3.8944e-5,3.8944e-5)} 
                    }
          
          # For easy to PDG values
          trueValue = JI.trueValue(self.data,param,hiera)[0]
          data = self.data[hiera]
   
          # Defines the coarse and precise linspaces
          start,stop = bounds[param]["coarse" if coarse else "precise"]    
          if coarse:
               EPSILON = np.linspace(start,stop,nb_points)  
          else:
               den_sigma = 1.5
               EPSILON = adaptive_linspace(start,stop,nb_points,np.mean([ep_mm,ep_mp]),abs(ep_mm-ep_mp)/den_sigma,np.mean([ep_pm,ep_pp]),abs(ep_pm-ep_pp)/den_sigma)

          if coarse:
               fig1 = plt.figure()
               ax1 = plt.gca() 

          # Discretized expected antineutrino flux
          histo_expected = np.histogram(self.E_vis,bins=bins,weights=self.S_expected.copy()/bins)[0]

          # Loop over the linspace
          chi2_stats = [] 
          for loop,ep in enumerate(EPSILON):
               print(f"{param} ; {self.time_step_days} : index epsilon {loop+1} out of {len(EPSILON)}")

               # Initializes the hierarchy parameters          
               if hiera==0:
                    m2 = data.m_32_NH
                    m3 = data.m_31_NH
               elif hiera==1:
                    m2 = data.m_23_IH
                    m3 = data.m_13_IH
               m1 = data.m_21
               s12 = data.s_12
               s13 = data.s_13

               # Adds epsilon to chosen parameter
               if param == "m21":
                    m1 += ep
               elif param == "s12":
                    s12 += ep 
               elif param == "s13":
                    s13 += ep
               elif (param=="m31_NH" and hiera==0) or (param=="m13_IH" and hiera==1):
                    m3 += ep
               else:
                    print(f"Parameter {param} incompatible with hierarchy {hiera}")

               # Creates hierarchy object with one parameter shifted
               h = hierarchy(hiera,[data.inverted,
                                    [m1, m2, m3],
                                    [s12, s13, data.s_23],
                                    data.delta_radians
                                    ])
          
               # Create instances of FluxEstimation with the hierarchy depending on its ordering
               FluxStudy_obs = FluxEstimation(h if hiera == 0 else None, None if hiera == 0 else h, self.Cores, self.Isotopes, self.time_step_days)

               # Computes the observed flux with a parameter shifted (see FluxEstimation.py)
               # option needplot=False means that we don't want to lose time ploting graphs
               FluxStudy_obs.plot_flux(needplot=False) # ideal flux
               FluxStudy_obs.plot_visible_spectrum(needplot=False) # applies detector effects

               # Retrieves the flux for large PMTs for the correct hierarchy
               S_obs = np.real(FluxStudy_obs.S_vis_NH_L if hiera==0
                               else FluxStudy_obs.S_vis_IH_L).copy()  

               if (needplot or GUI) and coarse:
                    ax1.plot(self.E_vis,S_obs,linewidth=0.5)

               # Discretized observed antineutrino flux
               histo_obs = np.histogram(self.E_vis,bins=bins,weights=S_obs.copy()/bins)[0]
             
               # Compute Chi² for the shifted parameter
               chi2_stats.append( np.sum(np.power(histo_obs-histo_expected,2)/histo_expected) )
       
          chi2_stats = np.array(chi2_stats)

          # Plots the flux (expected and observed) 
          if (needplot or GUI) and coarse:
               txtcoarse = "(coarse)" if coarse else ""
               ax1.plot(self.E_vis, self.S_expected, color="black", linewidth=1.5, label="Expected") # expected curve
               ax1.legend(fontsize=10)
               ax1.set_xlabel(param, fontsize=12)
               ax1.set_ylabel("Visible flux", fontsize=12)
               ax1.tick_params(axis='both', labelsize=10)      
               ax1.set_title(f"Flux comparison for {param} {txtcoarse}", fontsize=14)
               if not GUI:
                    plt.show() # Flux graph 

          # List of shifted parameters
          param_shifted = EPSILON + trueValue
        
          minchi2 = np.min(chi2_stats)
          # Check in case there is a nan (not a number) in the Chi² list (otherwise breaks everything)
          if str(min(chi2_stats))=='nan': 
               for k in range(len(chi2_stats)):
                    if str(chi2_stats[k])=='nan':
                         chi2_stats[k] = np.inf
               minchi2 = min(chi2_stats)

          # Subtracts the min value to the Chi² array 
          # This is why the min of DeltaCHI2 is 0 (as Chi² is positive)
          DeltaCHI2 = chi2_stats - minchi2

          # Epsilon intervals such that DeltaCHI2(epsilon) crosses 1
          ep_mm, imm, ep_mp, imp, ep_pm, ipm, ep_pp, ipp = search_for_epsilon(EPSILON, DeltaCHI2)
          ep_m_mean = np.mean([ep_mm, ep_mp])
          ep_p_mean = np.mean([ep_pm, ep_pp])

          print(f"epsilon-lower : {ep_mm}\nepsilon-upper : {ep_mp}\nepsilon-mean : {ep_m_mean}")
          print(f"epsilon+lower : {ep_pm}\nepsilon+upper : {ep_pp}\nepsilon+mean : {ep_p_mean}\n")

          if (needplot or GUI) and not coarse:
               fig2 = self.plot_chi2_parabola(param,hiera,param_shifted,DeltaCHI2,[ep_mm,ep_mp,ep_pm,ep_pp],[imm,imp,ipm,ipp],[ep_m_mean,ep_p_mean],start,stop,GUI=GUI)
          if not needplot:
               plt.close("all")
          if not GUI:
               return ep_mm, ep_mp, ep_m_mean, ep_pm, ep_pp, ep_p_mean
               
          if coarse:
               return ep_mm, ep_mp, ep_m_mean, ep_pm, ep_pp, ep_p_mean, fig1
          else:
               return ep_mm, ep_mp, ep_m_mean, ep_pm, ep_pp, ep_p_mean, fig2




     def Res_estim_multi(self, hiera, nb_points, param1_name, param2_name,
                         needplot=True, GUI=False):
          """
          Shows a two-dimensional oscillation parameter space and confidence levels contours assuming 
          pairs of correlated parameters.

          Oscillation parameters vary simultaneously, allowing to minimize the Chi² function and project 
          it onto two-dimensional spaces.
          From the map obtained, the computation of the Hessian matrix allows to go back to the correlation factor between
          pairs of parameters. 
          This way, we can quantify how correlated those parameters are.

          Handleling correlated parameters requires modifying the Chi² method.
          We need to introduce a covariance matrix, filled with both correlated and uncorrelated uncertainties (statistics and systematics).
          Those uncertainties refer to the flux, detection and background parts.
               - correlated flux uncertainties: 2%
               - uncorrelated flux:             0.8%    
               - detection:                     1%
               - background:                    1%
          To these systematics, we need to add statistical uncertainties corresponding to the standard deviation of count bin. 
          """
          bins=510
          # Ranges of the explored parameter space empirically adjusted
          ranges = {"s12":[0.3,0.314],
                    "s13":[1e-2,3.55e-2],
                    "m21":[7.42e-5,7.64e-5],
                    "m31_NH":[2.505e-3,2.555e-3],
                    "m13_IH":[2.43e-3,2.48e-3]}
          param1_range = ranges[param1_name]
          param2_range = ranges[param2_name]

          data = self.data[hiera] # For easy access to PDG values

          
          fixed_params = {"s12": data.s_12, "s13": data.s_13, "m21": data.m_21, "m31_NH": data.m_31_NH} \
               if hiera==0 else {"s12": data.s_12, "s13": data.s_13, "m21": data.m_21, "m13_IH": data.m_13_IH}
          
          # Initializes the two parameters that will vary 
          fixed_params[param1_name] = None  
          fixed_params[param2_name] = None  
          param1_values = np.linspace(param1_range[0], param1_range[1], nb_points)
          param2_values = np.linspace(param2_range[0], param2_range[1], nb_points)

          # Discretized expected antineutrino flux
          histo_expected = np.histogram(self.E_vis, bins=bins, weights=self.S_expected.copy()/bins)[0]

          # Initializes the Chi² matrix
          chi_square_results = np.zeros((nb_points, nb_points))  

          # Loop over the two parameters
          for i, param1 in enumerate(param1_values):
               for j, param2 in enumerate(param2_values):
                    print(f"1: {i} out of {len(param1_values)}   ;   2: {j} out of {len(param2_values)}")
                    # Shift the two parameters 
                    fixed_params[param1_name] = param1
                    fixed_params[param2_name] = param2

                    # Creates hierarchy object with two parameters shifted
                    h = hierarchy(hiera,[data.inverted,
                                         [fixed_params["m21"],
                                          data.m_32_NH if hiera == 0 else data.m_23_IH,
                                          fixed_params["m31_NH"] if hiera == 0 else fixed_params["m13_IH"]],
                                          [fixed_params["s12"], fixed_params["s13"], data.s_23],
                                          data.delta_radians])
        

                    # Create instances of FluxEstimation with the hierarchy depending on its ordering
                    FluxStudy_obs = FluxEstimation(h if hiera == 0 else None, None if hiera == 0 else h, self.Cores, self.Isotopes, self.time_step_days)
                    
                    # Computes the observed flux with two parameter shifted (see FluxEstimation.py)
                    FluxStudy_obs.plot_flux(needplot=False)
                    FluxStudy_obs.plot_visible_spectrum(needplot=False)

                    # Retrieves the flux for large PMTs for the correct hierarchy
                    S_obs = np.real(FluxStudy_obs.S_vis_NH_L if hiera == 0 else FluxStudy_obs.S_vis_IH_L).copy()
                    
                    # Discretized observed antineutrino flux
                    histo_obs = np.histogram(self.E_vis, bins=bins, weights=S_obs.copy()/bins)[0]

                    # Correlation matrix filled with both correlated and uncorrelated uncertainties 
                    V = np.zeros((bins, bins))
                    sigma_flux_c, sigma_flux_nc, sigma_detect, sigma_bkg = 0.02, 0.008, 0.01, 0.01 
                 
                    for k in range(bins):
                         for l in range(bins):
                              V[k, l] += sigma_flux_c**2 * histo_obs[k] * histo_obs[l] + sigma_detect**2 * histo_obs[k] * histo_obs[l]
                         V[k, k] += histo_obs[k] + sigma_flux_nc**2 * histo_obs[k]**2 + sigma_bkg**2 * histo_obs[k]**2

                    # Compute Chi² for the shifted parameter
                    chi2 = (histo_expected - histo_obs).T @ np.linalg.pinv(V) @ (histo_expected - histo_obs)
                    chi_square_results[i, j] = chi2  
        
          minchi2 = np.min(chi_square_results)
          delta_chi2 = np.array(chi_square_results) - minchi2

          if param1_name != param2_name:  # avoid inverting the correlation matrix which is not inversible when we consider the same parameter
               # Part where we estimate the correlation factor, 
               # by computing the Hessian matrix corresponding to the DeltaChi² variations and inverting it to go back to the correlation matrix 
               dx = param1_values[1] - param1_values[0] 
               dy = param2_values[1] - param2_values[0] 

               # Second derivatives using finite differences
               d2chi_dx2 = (delta_chi2[2:, 1:-1] - 2*delta_chi2[1:-1, 1:-1] + delta_chi2[:-2, 1:-1])/dx**2
               d2chi_dy2 = (delta_chi2[1:-1, 2:] - 2*delta_chi2[1:-1, 1:-1] + delta_chi2[1:-1, :-2])/dy**2
               d2chi_dxdy = (delta_chi2[2:, 2:] - delta_chi2[2:, :-2] - delta_chi2[:-2, 2:] + delta_chi2[:-2, :-2])/(4*dx*dy)
          
               # Means to approximate Hessian elements 
               H11 = np.mean(d2chi_dx2)
               H22 = np.mean(d2chi_dy2)
               H12 = np.mean(d2chi_dxdy)

               H = np.array([[H11, H12], [H12, H22]])  # Hessian matrix
               C = np.linalg.inv(H)                    # Corelation matrix
          
               rho = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]) # Correlation 
          else:
               rho = 1
          print(f"Estimated correlation factor between {param1_name} and {param2_name}: {rho:.5f}")

          if needplot:
               fig, ax = plt.subplots()
               X, Y = np.meshgrid(param1_values, param2_values)
               sigma_levels = [2.30, 6.18, 11.83, 18.42, 28.74] # DeltaChi² values corresponding to a significance level for 2D parameter space
               c = ax.pcolormesh(X, Y, delta_chi2.T, cmap='coolwarm', shading='auto', vmin=0, vmax=30, alpha = 0.5) 
               contours = ax.contour(X, Y, delta_chi2.T, levels=sigma_levels, colors=['black', 'black', 'black', 'black', 'black'], linewidths=1.5, linestyles=['solid', 'solid', 'solid', 'solid', 'solid'])
               ax.clabel(contours, fmt={2.30: r'1$\sigma$', 6.18: r'2$\sigma$', 11.83: r'3$\sigma$', 18.42: r'4$\sigma$', 28.74: r'5$\sigma$'}, inline=True, fontsize=10)
               cbar = fig.colorbar(c, ax=ax)
               cbar.set_label(r'$\Delta\chi^2$',fontsize=10) 

               ax.scatter(JI.trueValue(self.data,param1_name,hiera)[0],
                          JI.trueValue(self.data,param2_name,hiera)[0],
                          s=70, color="red", marker='x', linewidth=1, label='Central PDG value')
               ax.set_xlabel(param1_name, fontsize=12)
               ax.set_ylabel(param2_name, fontsize=12)
               ax.tick_params(axis='both', labelsize=10)
               ax.ticklabel_format(axis='both', style='sci')
               ax.set_title(r'$\Delta\chi^2$ CL Contour Map', fontsize=14)

               ax.set_xlim(param1_range[0], param1_range[1])
               ax.set_ylim(param2_range[0], param2_range[1])
               ax.legend(fontsize=10)
        
               ax.set_box_aspect(1)
     
               if not GUI:
                    plt.show()
                    return chi_square_results
               else:
                    return fig, rho

