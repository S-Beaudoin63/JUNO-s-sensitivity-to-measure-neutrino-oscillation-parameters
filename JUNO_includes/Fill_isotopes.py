import numpy as np
import pandas as pd

class isotope:
    """
    Loads the isotopes data from the CSV
    """
    def __init__(self,di,index):
        self.name = di["name"][index]
        self.Ereleased = di["Ereleased"][index]
        self.FissRate = di["FissRate"][index]
        self.FuelFrac = di["FuelFrac"][index]
        self.FluxesPerFission = np.array([di["Nb_nu_params1"][index],di["Nb_nu_params2"][index],di["Nb_nu_params3"][index],
                                          di["Nb_nu_params4"][index],di["Nb_nu_params5"][index],di["Nb_nu_params6"][index]])

def Dic_of_isotopes():
    """
    Returns the dictionnary containing the isotopes object for all considered isotopes
    """
    Isotopes={} # dictionnary of core objects
    di = pd.read_csv("./data/isotopes.csv")
    for iso_index in range(len(di)):
        Isotopes[di["name"][iso_index]] = isotope(di,iso_index)
    return Isotopes