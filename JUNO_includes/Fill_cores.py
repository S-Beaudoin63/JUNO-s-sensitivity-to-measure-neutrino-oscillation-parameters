import numpy as np
import pandas as pd

class core:
    """
    Loads the cores data from the CSV
    """
    def __init__(self,dc,index):
        self.name = dc["name"][index]
        self.Power = dc["Power"][index]
        self.Baseline = dc["Baseline"][index]
        self.PowerMeVd = (dc["MeV/1e25"][index]*24)*1e25
    
def Dic_of_cores():
    """
    Returns the dictionnary containing the cores object for all considered cores.
    We exclude TS-C3, TS-C4, HZ because they won't be operationnal.
    """
    Cores={} # dictionnary of core objects
    dc = pd.read_csv("./data/cores.csv")
    for core_index in range(len(dc)):
        if dc["name"][core_index] not in ["TS-C3","TS-C4","HZ"]: # those cores won't be operationnal
            Cores[dc["name"][core_index]] = core(dc,core_index)
    return Cores

