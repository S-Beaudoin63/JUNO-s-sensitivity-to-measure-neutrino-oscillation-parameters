# JUNO-s-sensitivity-to-measure-neutrino-oscillation-parameters
Master 2 computational physics project


## Introduction ##
The Jiangmen Underground Neutrino Observatory (JUNO) is a next-generation liquid scintillator detector
designed for precise neutrino oscillation measurement, located 52.2 km from several nuclear reactors in
China. Our results confirm that JUNO should be able to determine the oscillation parameters sin(θ12),
sin(θ13), ∆m12 and ∆m13 with world-leading precision under 0.5% as well as determining the neutrino mass
hierarchy with 3 - 4σ significance after at least 6 years of data-taking.



## Overview ##
This program is divided into several parts, each corresponding to a file in the folder JUNO_includes.
It takes data from three csv (Comma Separated Values) files in the data folder.

All functionnalities can be used from the main.py file which has been commented to make it intuitive.

To check the results, we also have made a Graphic User Interface (GUI) which can be run from the JunoApp.py
file or using the JunoApp.exe. The latter option is useful because it doesn't require to install libraries
used for the GUI. Both should be supported on Windows and Linux. The assets folder contains the files
required for the GUI to work. More on that at the end of this README.

There is a last folder, UML_diagram, containing UML diagrams to visualize classes and custom packages,
generated using pyreverse.



## Layout of the project  ##
All the codes related to the modeling and physical computations are stored in the JUNO_includes folder.
For the best-fit values of the parameters, we use the values from PDG 2024.

# 1 - Data loading #
The physical constants and the properties of the reactors and elements at play have been stored in three csv files in the data folder.
This data is loaded at the begining of the main (respectively in the constructor of the App) as they will be used pretty much everywhere.
The folder includes:

- cores.csv: contains the properties of the nuclear reactors close to the JUNO experiment (Power [GW], Baseline, PowerW/h).
    Note: not all cores are taken into account in the simulation as three of them (TS-C3, TS-C4, HZ) won't be operational
    during the experiment (see article for the parametrization of the reactor antineutrino flux [3], [4]).
    The csv is read in Fill_cores.py, creating an instance of the "core" class defined in the same file for all cores that
    need to be taken into account and storing them in a dictionnary by the function Dic_of_cores().

- isotopes.csv: contains the characteristics of the four isotopes present in the aforementionned cores (Energy released,
    Fission Rate, Fuel Fraction, the coefficients of the polynomial parametrizing the flux [3], [4]).
    The csv is read in Fill_isotopes.py, which works based on the same layout as Fill_cores.py: a function Dic_of_isotopes()
    dictionnary containing an isotope object for all considered isotopes.

- hierarchies.csv: contains PDG values for the antineutrino oscillation parameters depending on the hierarchy
    (s12=sin(θ12), s13=sin(θ13), s23=sin(θ23), m21=∆m2_12, m23=∆m2_23, m31_NH=∆m2_31 or m13_IH=∆m2_12, delta (CP phase).
    Note: we name differently the third squared masses differences depending on the hierarchy to always handle positive values.
    Indeed for the normal hierarchy, ∆m2_31 > 0 and ∆m2_13 <0 while for the inverted hierarchy, ∆m2_31 < 0 and ∆m2_13 > 0.



# 2 - Hierachies #
The hierarchy class beeing more complex than core and isotope, it deserves its own segment.
It is used to contain the PDG values for the oscillation parameters but also to compute quantities associated with
the two hierarchies (such as the PMNS matrix). It will later be used to create "custom" hierarchies with
parameters shifted from their PDG values to study the influence of parameters values on observable antineutrino phenomena.

The class consists of two "constructors" that are called depending on the context and their arguments.
Unlike in C++, python does not allow to surcharge methods using different arugments, so it is necessary to have a "real"
constructor hierarchy.__init__(hiera,X) that receives hiera (0 or 1 depending on if the hierarchy is inverted or not) and X
which can be either a string (path to hierarchies.csv) or a list.
Depending on the type of X, __init__ will call either:

- init_csv(string file, bool hiera):
    Uses PDG 2024 values, read in hierarchies.csv, to fill in a hierarchy object.
    Takes as argument the type of hierarchy and the path to the CSV (for compatibility between Linux and Windows)

- init_manu(bool inverted, list(float) Dmasses2, list(float) sins, float delta_radians):
        "Manual" constructor, uses custom values given as parameters, to fill in a hierarchy object
            - inverted: 0 or 1 if the hierachy is normal or inverted
            - Dmasses2: list of the 3 defining masses (m12, m23 and m31_NH or m13_IH)
            - sins: list of the 3 sines (s12,s23,s13)
            - delta_radians: the value in radians of the CP phase
          The hierarchy object has the same attributes than if it was created with the manual constructor

- PMNS_matrix():
    Generates the PMNS matrix modelling neutrinos mixing and oscillations (called in both constructors)
	/!\ Warnings "casting real/complex values" when discarding real/imaginary parts

- oscillation_anti(linspace energy, float baseline, flavor initial_state, flavor final_state): flavour is 0, 1 or 2 for electron/muon/tau
        Computes the oscillation probability from some initial antineutrino flavor to a final one,
        for a given energy, baseline and hierarchy.
        We consider only electron antineutrinos for the study of Juno



# 3 - Common variables #
Accross the different files, several variables will be redundant (they may also be found as attributes of various classes):
- hiera: int, 0 or 1 depending on if the hierarchy is inverted or not

- NH (hierarchy object): corresponds to the normal hierarchy filled with PDG 2024 values
- IH (hierarchy object): corresponds to the inverted hierarchy filled with PDG 2024 values
- data (dictionnary of hierarchy objects) = {0:NH, 1:IH}, used for easier access to PDG values

- Cores (dictionnary of core objects):          contains the properties of all used cores to the project
- Isotopes (dictionnary of isotope objects): contains the properties of the isotopes present in the cores of the nearby reactors

- year (float): 365.25 days/year

- param (string): names of the parameters accessible through JUNO's measurements ("s12", "s13", "m21", "m31_NH", "m13_IH")

- E, Energies (np.linspace): antineutrino energies, between 1.805 and 10 MeV
- E_vis (np.linspace):         antineutrino energies, visible to the detector, between 1.5 and 11 MeV

- Flux or S (np.array/list of float): antineutrino flux in function of the energy
- bins (int) = 510: the number of energy bins assumed from previous Daya Bay information

- nb_points (int): the number of times we shift a given parameter we evaluate when we estimate a resolution

- needplot (bool): indicates if we want to plot the results of a function, so that we don't lose time when we don't need them (default:True).
- GUI (bool): indicates that the method was called from the GUI and we want to return the figures instead of plt.showing them (default=False).



# 4 - Oscillations # (Proof of concept)
To validate the behavior of the hierarchy class methods, we had previously done a similar code for neutrinos, plot their
oscillation probabilities. This was a checkpoint to ensure that the code worked well as we could easily compare
the graphs obtained with litterature.

We kept this code as additional features, that allow to visualize the global properties of neutrino oscillations
as an introduction to the subject.
It can be used boht in the GUI and in the terminal:
osc_interface = neutrino_oscillations(NH,IH): creates an instance of the text interface in the terminal
osc_interface.plot_oscillations_interface():  launches latter interface to dialogue with the user



# 5 - Antineutrino flux estimation #
The first step of the project was to estimate the number of antineutrinos that would be detected
in JUNO, knowing the different nuclear reactor cores that emit them.
Each reactor core emits antineutrinos through the beta decay of fission products from the considered isotopes.

The FluxEstimation class models the antineutrino flux to be expected in the detector,
coming from the nuclear plants nearby:

- First we evaluate the characteristic flux for a given isotope:    S_iso(string isotope)
  and the contribution from inverse beta decay (approx. from [5]):    IBD_cross_section()
- Compute the ideal flux (without detector effects):                flux_tot()
- Then apply the detector effects:
    Energy resolution of each photomultiplier tube:                 energy_resolution(string PMT_key)
        Maps emitted energies to visible detector energies:            detector_response(string, PMT_key)
        Apply all effects to get a realistic antineutrino flux:        simulate_visible_spectrum(flux, string PMT_key)
    the variable PMT_key is a string to identify photomultiplier tubes

The result of this analysis is a visible antineutrino flux which gives about 50 antineutrinos detected per day.
The visible energy resolution is defined by a gaussian.
It takes 3 contributions [8]:
        - 1 driven by the Poisson statistic of the number of detected photoelectrons
        - 1 dominated by the detector's spatial non uniformity
        - 1 dominated by dark noise

This analysis returns three plots:
- the number of antineutrino per fission for each isotope
- the ideal flux with and without oscillations
- the observed flux for each type of photomultiplier tube and hierarchy, compared with the total ideal flux

and returns the number of antineutrinos
- ideal (withouth oscillations, for both hierarchies)
- observed (both hierarchies, both PMTs)



# 6 - Resolution Estimation #
Then, the goal of the project involves determining the relative precision obtained by JUNO for each
oscillation parameter of interest after some data-taking time.

The class ResolutionEstimation is used to compute such relative precisions.

We initialize a study with the constructor:
    ResolutionEstimation(float time_step_days, array/list S_expect, dictionnary data, dict Cores, dict Isotopes):
where time_step_days is the data-taking time in days at which we will evaluate resolutions.


There are two types of such studies:

1) The first method assumes uncorrelated parameters at a certain time self.time_step_days (float, initialized in the constructor):
    Res_estim(self, param, hiera, nb_points,
        coarse=False, ep_mm=None, ep_mp=None, ep_pm=None, ep_pp=None, needplot=True, GUI=False)

This method aims at finding the one sigma uncertainty on a given parameter using a Chi2 procedure
between an expected binned flux and a shifted parameter flux.
To do so it needs to explore different shifted parameter values.

It is called twice in (PerformanceStudy.full_resolution_computation):
        - a first time with coarse = True. It explores a small linspace to approximately locate the values of the
          shift epsilon (ep) such that DeltaChi2=1 crossing points
        - then, ep_mm, ep_mp, ep_pm and ep_pp are filled with the coarse estimates
          and we reapply the method with a double gaussian adaptative linspace centered on the previously
          estimated ep points.
This allows to use a small amount of points while assuring a good accuracy.

Returns:
the bounding and mean values of:
- the negative uncertainty
- the positive uncertainty
and two plots (if needplot=True):
- a plot to compare the flux obtained with PDG values with the flux obtained by shifting one parameter by a constant amount
- a plot showing the DeltaChi2 parabola, to compare it to the one obtained with PDG uncertainties and those from the article [2]


2) The second method allows to verify the correlation factors between parameters:
     Res_estim_multi(self, hiera, nb_points, param1_name, param2_name,
                         needplot=False, GUI=False):

    - we select a pair of parameters, the other two are fixed to their best-known values
    - we make the two parameters vary and draw a Chi2 map allowing to represent the different significance regions
        according to 2D Wilks' theorem
    - then we relate it to the correlation factor through the computation of the Hessian matrix

    Returns (if needplot=True) a two-dimensional oscillation parameter space and confidence levels contours assuming
    pairs of correlated parameters.

    Two oscillation parameters vary simultaneously, allowing to minimize the Chi2 function
    and project it onto two-dimensional spaces.

    From the map obtained, the computation of the Hessian matrix allows to go back to the correlation
    factor between pairs of parameters.
    This way, we can quantify how correlated those parameters are.

    Handleling correlated parameters requires modifying the Chi2 method.
    We need to introduce a covariance matrix, filled with both correlated and uncorrelated uncertainties (statistics and systematics).
    Those uncertainties refer to the flux, detection and background parts.
        - correlated flux uncertainties: 2%
            - uncorrelated flux:             0.8%    
            - detection:                     1%
            - background:                    1%
    To these systematics, we need to add statistical uncertainties corresponding to the standard deviation of count bin.



# 7 - Complete analysis #
The class PerformanceStudy automates different studies:

- full_resolution_computation(str param, bool hiera, float time, int nb_points=20,
                                 needplot=True, GUI=False):
    Automates the estimation of the relative precision for a given parameter and time using Res_estim.
    1) estimates expected flux with detector effects,
    2) realizes a coarse resolution estimation
    3) then a precise resolution estimation,
         
    Regarding nb_points, it's the number of points used in the precise estimation.
    The one used for the coarse estimation is deduced from it using f(x) = 20(1-exp(-nb_points*ln(2)/20)
    so that, if nb_points = 20: launches a coarse study of 10 points into a precise study of 20 points (empirical).

    Returns:
    - the estimates of the resolution returned by the precise Res_estim call
    - if needplot=True:
        - the flux comparison plot returned by the coarse Res_estim call
        - the plot with the parabolas returned by the precise Res_estim call


- loop_time(self, param, hiera, TIME, nb_points,
                  needplot=True, GUI=False):

    Runs full_resolution_computation for several data-taking times
    to get the time evolution of the relative precision of one parameter.

    Returns:
        - the list of resolutions, each returned by a call to full_resolution_computation for
          a different data-taking time
        - if needplot=True: the plot showing the evolution of the relative precision
          through time for the chosen parameter


- loop_params(self, Lparam, hiera, TIME, nb_points,
                    needplot=True, GUI=False):

    Runs loop_time for all parameters of the given hierarchy
    to get the time evolution of the relative precision for all parameters.
    
    if needplot=True: returns the plot of the evolution of the relative precisions
    through time for all parameters.


- multiparam_study(self, param1, param2, hiera, time,
                         nbpoints=10, needplot=True, GUI=False):
    
    Automates the call to Res_estim_multi
    by computing first the expected flux with detector effects.
    This is the only case where nbpoints=10 by default instead of 20 because it takes 4 times longer, but it can be modified and runs.
    
    Returns the two-dimensional oscillation parameter space and confidence level contours assuming
    pairs of correlated parameters, returned by Res_estim_multi.




- def analyze_hierarchy_distinction(self, time_min, time_max, nb_points,
                      series=5, GUI=False):

    The final goal of the project.
    Estimates the minimum data-taking time required for JUNO to be able to
    distinguish between both hierarchies with at least 3 sigmas.

    The method works as follows:
        - we assume NH to be True (although the same reasoning with IH yields the same results)
                and we compute all the associated binned spectra.
        - then we compute all IH spectra and add a gaussian noise in order to mimick statistical fluctuations
            - we compare NH and IH noisy through a Chi2 method
            - we get the p-value and associated significance thresholds using Wilks' theorem  
            - we loop over time to obtain the significance evolution
            - we repeat the process 'series' times to average over all simulations and get a more stable result        

    This method is the only using random functions in the project so we introduced a
    series argument that is the number of runs to average on.
    
    Returns:
    - the data-taking time it takes to detect a difference between both hierarchies at 3 and 4 sigma.
    - if needplot=True: the plot of the evolution of the significance through time.



# 8 - An example of the influence of hierarchies: the neutrinoless double beta decay #
The class DoubleBetaExample is used to provide an example of application of the hierarchy
parameter constraints and the difference that hierarchies make.

Studying neutrinoless double beta decay allows to probe the nature of neutrinos (Dirac or Majorana fermions).

The parameter of interest is the effective Majorana mass which is a combination of the neutrino
mass eigenstates and the neutrino mixing matrix terms.

The effective mass can be expressed as a function of the lightest neutrino mass [17] and relates to the inverse
of the processe's lifetime.

beta_graph(GUI=False):
returns the plot of the accessible parameter space depending on the hierarchy,
providing an example of consequence from Juno's hierarchy measurements and parameter constraints.


##############

## Utils ##
Utils.py contains utility functions that are called in other functions and methods in diverse locations
of the project.


## JUNO APP ##

To more easily check the results, we made a Graphic User Interface (GUI) which can be run from the JunoApp.py
file or using the JunoApp.exe. The latter option is useful because it doesn't require to install libraries
used for the GUI. Both should be supported on Windows and Linux. The assets folder contains the files
required for the GUI to work.

JunoApp is a Python application developed using the Customtkinter library.

When using the GUI, most plots will display almost instantly (except Flux Estimation and Double Beta Example).
Indeed we have pregenerated plots that would take too much time to generate so that the GUI would be a
swift way to quickly see the plots. If you want to check that the code works, prefer using main.py.
However, as there may be too many possible parameter combinations, we won't have generated every plot possible.
Moreover when using analyze_hierarchy_distinction, it is random so a dropdown offers to either load an
existing plot or, by choosing None, the App will generate a new plot which may take a long time.
Whenever a plot takes too long (outside of the two special cases mentioned above), check the terminal
(either you executed the program from a terminal, or if you launched the executable by double clicking
on it, a terminal should have appeared) if the program seems doing somme kind of loop and you can stop it
with CTRL+C in said terminal if you want. If you let it run, when it's over it will save the data and
if you try later the same combination of parameters, it should be displayed instantly because you ran it
previously.

We have pre-generated most of the graphs with default values,i.e. with 20 points (10 for multiparametric study). 
So if you want to generate a new graph, pick a different number of points.

/!\ When using the GUI, consequent warnings tend to appear in the terminal, they are caused by somme 
incompatibilities between graphical libraries and do not affect the program whatsoever.



# In case the executable is not supported by your device
Requirements (to run JunoApp.py)
Python 3.x

altgraph==0.17.4
chardet==5.2.0
contourpy==1.3.1
CTkMessagebox==2.7
CTkToolTip==0.8
customtkinter==5.2.2
cycler==0.12.1
darkdetect==0.8.0
dill==0.3.9
fonttools==4.56.0
kiwisolver==1.4.8
matplotlib==3.10.0
numpy==2.2.3
packaging==24.2
pandas==2.2.3
pefile==2023.2.7
pillow==11.1.0
pip==24.3.1
pyparsing==3.2.1
python-dateutil==2.9.0.post0
pytz==2025.1
pywin32-ctypes==0.2.3
reportlab==4.3.1
scipy==1.15.1
setuptools==75.8.0
six==1.17.0
tzdata==2025.1


Install the required libraries using pip:
pip install matplotlib numpy pandas pillow scipy customtkinter openpyxl CTkMessagebox CTkToolTip



# License
This software is provided as-is without any warranty. You are free to use, modify, and distribute it as needed.

## Authors ##
Simon Beaudoin, Nathan Lebas and Pierre Louis--Cistac
Université de Strasbourg, IPHC, 23 rue du Loess
67037 Strasbourg, France
