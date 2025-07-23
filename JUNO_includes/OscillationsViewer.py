import matplotlib.pyplot as plt
import numpy as np

from JUNO_includes.Hierarchy import hierarchy

class neutrino_oscillations:
    """
    This additional code allows to visualize global properties of neutrino oscillations 
    as an introduction to the subject.
    It is better used from the GUI but it is still possible to use from the main with:
        osc_interface = neutrino_oscillations(NH,IH)
        osc_interface.plot_oscillations_interface()  
    """

    def __init__(self, NH, IH):
        self.NH = NH
        self.IH = IH
        self.flavours = ["e", "µ", "τ"]

    def oscillation_verif(self, hiera, ratio, a, b) :
        """
        Computes the oscillation probability from some initial neutrino flavor a to a final one b, 
        for a given energy/baseline ratio and hierarchy.
        This is a general method allowing to consider all flavor pairs. 

        This method is similar to Hierarchy.oscillation_anti but we only use it for oscillation visualization
        and not for other computations. 
        Thus, its structure is slightly simplified 
        """
        U = hiera.PMNS

        if hiera.inverted == 0 :   # NH
            mass = [[0, -hiera.m_21, -hiera.m_31_NH], [hiera.m_21, 0, -hiera.m_32_NH], [hiera.m_31_NH, hiera.m_32_NH, 0]]
        elif hiera.inverted == 1 : # IH
            mass = [[0, -hiera.m_21, hiera.m_13_IH], [hiera.m_21, 0, hiera.m_23_IH], [-hiera.m_13_IH, -hiera.m_23_IH, 0]]
        else :
            print("Please select a hierarchy")

        P_list = []     # List of oscillation probabilities
        den = (U @ np.conj(U.T))[a][a]*(U @ np.conj(U.T))[b][b]     # Denominator of the probability (normalisation factor)
        for r in ratio :
            P = 0
            for i in range(0, 3) :
                P += np.abs(U[a][i])**2*np.abs(U[b][i])**2
                for j in range(i+1, 3) :
                    P += 2*np.real(U[a][i]*U[b][j]*np.conj(U[a][j])*np.conj(U[b][i]))*np.cos(mass[i][j]*r*2.54) - 2*np.imag(U[a][i]*U[b][j]*np.conj(U[a][j])*np.conj(U[b][i]))*np.sin(mass[i][j]*r*2.54)
            P_list.append(P/den)
        return P_list


    def plot_oscillations_interface(self):
        """
        Launches a text interface in the terminal for the user to choose the oscillation pattern
        """
        color = ['green', 'blue', 'red']
        while True:
            print("Do you want to:")
            print("1: Explore each oscillation individually")
            print("2: Compare hierarchies")
            print("3: Quit")
            mode = int(input("Enter your choice (1, 2, or 3): "))

            if mode == 3:
                print("Exiting the program. Goodbye!")
                break

            print("Select initial flavour:")
            for i, flavour in enumerate(self.flavours) :
                print(f"{i}: {flavour}")
            initial = int(input("Enter the index for the initial flavour:"))

            if mode == 1:
                print("Select final flavour:")
                for i, flavour in enumerate(self.flavours):
                    print(f"{i}: {flavour}")
                final = int(input("Enter the index for the final flavour (or -1 to plot all):"))
                hierarchy = int(input("Enter hierarchy (0 for NH, 1 for IH):"))

                h = self.NH if hierarchy==0 else self.IH

                ratio = np.linspace(0, 36000, 4096)

                if final == -1:
                    # Plot all final flavours for the selected initial flavour
                    plt.figure()
                    for final in range(0, 3):
                        label = f"P({self.flavours[initial]} → {self.flavours[final]})"
                        plt.plot(ratio, self.oscillation_verif(h, ratio, initial, final), label=label, color = color[final], alpha = 0.6)
                    plt.title(f"Neutrino oscillation probabilities starting from {self.flavours[initial]}")
                    plt.xlabel("L/E [km/GeV]")
                    plt.ylabel("Probability")
                    plt.legend()
                    plt.show()

                else:
                    # Plot a single transition
                    label = f"P({self.flavours[initial]} → {self.flavours[final]})"
                    plt.plot(ratio, self.oscillation_verif(h, ratio, initial, final), label=label, color = color[final], alpha = 0.6)
                    plt.title(f"Neutrino oscillation probability: {self.flavours[initial]} → {self.flavours[final]}")
                    plt.xlabel("L/E [km/GeV]")
                    plt.ylabel("Probability")
                    plt.legend()
                    plt.show()

            elif mode == 2:
                print("Select final flavour:")
                for i, flavour in enumerate(self.flavours):
                    print(f"{i}: {flavour}")
                final = int(input("Enter the index for the final flavour: "))
                ratio = np.linspace(0, 36000, 4096)

                # Compare hierarchies for the selected initial and final flavours
                for hierarchy in [0, 1]:
                    h = {0:self.NH, 1:self.IH}[hierarchy]
                    color = 'blue' if hierarchy == 0 else 'red'
                    label = f"P({self.flavours[initial]} → {self.flavours[final]}) ({'NH' if hierarchy == 0 else 'IH'})"
                    plt.plot(ratio, self.oscillation_verif(h, ratio, initial, final), label=label, color = color, alpha = 0.6)

                plt.title(f"Neutrino oscillation probabilities: {self.flavours[initial]} → {self.flavours[final]} (NH vs IH)")
                plt.xlabel("L/E [km/GeV]")
                plt.ylabel("Probability")
                plt.legend()
                plt.show()
        else:
            print("Invalid choice. Please select 1, 2, or 3.")

    def plot_oscillations_GUI(self, mode, hiera, init_flavor, final_flavor):
        """
        Same function but addapted for the GUI 
        """
        color = ['green', 'blue', 'red']

        fig = plt.figure()
        ax = plt.gca()

        if mode == 1:
             
            h = {0:self.NH,1:self.IH}[hiera]
            ratio = np.linspace(0, 36000, 4096)

            if final_flavor == -1:
                # Plot all final flavours for the selected initial flavour
                for final in range(0, 3):
                    label = f"P({self.flavours[init_flavor]} → {self.flavours[final]})"
                    ax.plot(ratio, self.oscillation_verif(h, ratio, init_flavor, final), label=label, color = color[final], alpha = 0.6, linewidth=1.5)
                ax.set_title(f"Neutrino oscillation probabilities starting from {self.flavours[init_flavor]}", fontsize=14)

            else:
                # Plot a single transition
                label = f"P({self.flavours[init_flavor]} → {self.flavours[final_flavor]})"
                ax.plot(ratio, self.oscillation_verif(h, ratio, init_flavor, final_flavor), label=label, color = color[final_flavor], alpha = 0.6, linewidth=1.5)
                ax.set_title(f"Neutrino oscillation probability: {self.flavours[init_flavor]} → {self.flavours[final_flavor]}", fontsize=14)


        elif mode == 2:
            ratio = np.linspace(0, 36000, 4096)

            # Compare hierarchies for the selected initial and final flavours
            for final_hierarchy in [0, 1]:
                h = {0:self.NH, 1:self.IH}[final_hierarchy]
                color = {0:"blue",1:"red"}[final_hierarchy]
                label = f"P({self.flavours[init_flavor]} → {self.flavours[final_flavor]}) ({'NH' if final_hierarchy == 0 else 'IH'})"
                ax.plot(ratio, self.oscillation_verif(h, ratio, init_flavor, final_flavor), label=label, color = color, alpha = 0.6, linewidth=1.5)

            ax.set_title(f"Neutrino oscillation probabilities: {self.flavours[init_flavor]} → {self.flavours[final_flavor]} (NH vs IH)", fontsize=14)
 
        ax.tick_params(axis='both', labelsize=10)
        ax.set_xlabel("L/E [km/GeV]", fontsize=12)
        ax.set_ylabel("Probability", fontsize=12)
        ax.legend(fontsize=10)
        
        return fig