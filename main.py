import sys
if "./JUNO_includes" not in sys.path:
    sys.path.append("./JUNO_includes")
print("\n")
from JUNO_includes.Fill_cores import Dic_of_cores
from JUNO_includes.Fill_isotopes import Dic_of_isotopes
from JUNO_includes.Hierarchy import hierarchy
from JUNO_includes.FluxEstimation import FluxEstimation
from JUNO_includes.PerformanceStudy import PerformanceStudy 
from JUNO_includes.OscillationsViewer import neutrino_oscillations
from JUNO_includes.DoubleBetaExample import DoubleBetaExample

"""
main.py is made so you only have to uncomment the lines 
corresponding to the type of analysis you want to run.
All functions called are commented in their respective files.
All functions are accessible with pre-generated graphs in the GUI
to save time for just checking the graphs. 
But if you want to check the well-functionning of the data generation, 
main.py will generate new graphs. 
"""


# Creates instances of hierarchy class for easy access to PDG 2024 values.
NH = hierarchy(0,["./data/hierarchies.csv"]) # Normal hierarchy
IH = hierarchy(1,["./data/hierarchies.csv"])  # Inverted hierarchy

# Loads information about the cores : Name, Power [GW], Baseline [km]
Cores = Dic_of_cores()


# Loads information about the isotopes : Name, Energy released per fission (MeV), fission rate, Fraction of each fuel component
# The neutrino flux per energy per fission for an isotope can be approximated as the exponential of a polynomial
Isotopes = Dic_of_isotopes()

year = 365.25 # Number of days in a year
TIME =  [100, 6*365.25, 20*365.25] 

"""
1 - Visualisation of neutrino oscillations
Runs an interface in the terminal for the user to select an oscillation phenomenon to study,
as well as the initial and final flavour.
"""
osc_interface = neutrino_oscillations(NH, IH)
#osc_interface.plot_oscillations_interface()  



"""
2 - Antineutrino flux estimation
Estimates the antineutrino flux to be expected in JUNO at a given time.
Plots:
    - number of antineutrino per fission and energy for all considered isotopes
    - ideal flux without detectors effects
    - effective flux, taking detectors effects into account
"""
t = 6*year
FluxStudy = FluxEstimation(NH, IH, Cores, Isotopes, t)
#FluxStudy.plot_Nb_nu_pfission()  
#FluxStudy.plot_flux(needplot=True)

## Detector response
#FluxStudy.plot_visible_spectrum(needplot=True)




"""
3 - Resolution Estimation
Estimates the relative precision on one parameter at a given time.
"""
param = "s12" # "s12" ; "s13" ; "m21" ; "m31_NH" if hiera=0 ; "m13_IH" if hiera=1
hiera = 0 # 0 (NH) or 1 (IH) 
t = 6*year # 100 ; 6*year ; 20*year will be compared with the resolutions from [2] but any other value should work 
Full_ResolutionEstimation = PerformanceStudy(NH, IH, Cores, Isotopes)
#Full_ResolutionEstimation.full_resolution_computation(param, hiera, t, 20, needplot=True) 
 



"""
4 - Time evolution of the relative precision for one parameter
"""
param = "s12" # "s12" ; "s13" ; "m21" ; "m31_NH" if hiera=0 ; "m13_IH" if hiera=1
hiera = 0 # 0 (NH) or 1 (IH)
TimeLoop = PerformanceStudy(NH, IH, Cores, Isotopes)
#Precisions = TimeLoop.loop_time(param, hiera, TIME, 20, needplot=True)




"""
5 - Time evolution of the relative precision for all parameters
Loops previous method (4) for all parameters
"""
Lparam = {0:["s12" , "s13" , "m21", "m31_NH"] , 1:["s12" , "s13" , "m21", "m13_IH"]}
hiera = 0 # 0 (NH) or 1 (IH)
ParamLoop = PerformanceStudy(NH, IH, Cores, Isotopes)
# Precisions2 = ParamLoop.loop_params(Lparam, hiera, TIME, 20, needplot=True)




"""
6 - Correlated parameters analysis
Better handling of correlation and systematics.
Visualizes the Chi² function over 2D projections of the parameter space,
while printing estimates of the correlation coefficients.
"""
param1 = "s12" # "s12" ; "s13" ; "m21" ; "m31_NH" if hiera=0 ; "m13_IH" if hiera=1
param2 = "s13" # "s12" ; "s13" ; "m21" ; "m31_NH" if hiera=0 ; "m13_IH" if hiera=1
hiera = 0 # 0 (NH) or 1 (IH)
t = 6*year
Resolution_multi = PerformanceStudy(NH, IH, Cores, Isotopes)
#Resolution_multi.multiparam_study(param1, param2, hiera, t, nbpoints=10, needplot=True)




"""
7 - Distinction between hierarchies
Estimates the minimum data-taking time required for JUNO to be
able to distinguish between hierarchies (at 3 and 4 sigma).
/!\ This method uses random so we average on a certain number of runs.  
"""
start = 100 # Starting point of time interval (days)
end = 10*year # Ending point of time interval (days)
nb_points = 20 # Number of points in the time interval
series = 5 # Number of runs to average on 
diff_hieras = PerformanceStudy(NH, IH, Cores, Isotopes)
#diff_hieras.analyze_hierarchy_distinction(start, end, nb_points, series=series) 




"""
Bonus: illustrates a broader particle physics example of 
experiment affected by the outcome of JUNO,
by showing the accessible parameter space of the invariant 
mass of the neutrinoless double beta decay depending on the hierarchy.
"""
betaStudy = DoubleBetaExample()
betaStudy.DoubleBetaExample_graph()
