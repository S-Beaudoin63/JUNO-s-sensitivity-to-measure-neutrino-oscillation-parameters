from tkinter import PhotoImage
from tkinter import Canvas as tCanvas
import customtkinter as ctk
import matplotlib.pyplot as plt
import pandas as pd
from CTkToolTip import CTkToolTip
from CTkMessagebox import CTkMessagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
import os
import platform #check os
import turtle
import pickle
import matplotlib.figure
import matplotlib.gridspec
matplotlib.gridspec.SubplotParams = matplotlib.figure.SubplotParams
import sys

# import custom modules
if "./JUNO_includes" not in sys.path:
    sys.path.append("./JUNO_includes")
if "./JUNO_includes/GUI_includes" not in sys.path:
    sys.path.append("./JUNO_includes/GUI_includes")

from JUNO_includes.Fill_cores import Dic_of_cores
from JUNO_includes.Fill_isotopes import Dic_of_isotopes
from JUNO_includes.Hierarchy import hierarchy
import JUNO_includes.Utils as JI
from JUNO_includes.FluxEstimation import FluxEstimation
from JUNO_includes.PerformanceStudy import PerformanceStudy
from JUNO_includes.OscillationsViewer import neutrino_oscillations
from JUNO_includes.DoubleBetaExample import DoubleBetaExample

from assets.GUI_includes.LazyImageLoader import LazyImageLoader
from assets.GUI_includes.SlidePanel import SlidePanel

# import App's assets
if "./assets" not in sys.path:
    sys.path.append("/assets")
from assets.dico import dico_Textes

language = "English"


class JunoApp:
    """
    This class is used to launch a GUI to more easily view the graphs.
    As it is an additional part of the project and only deals with ploting data
    and saving it, it is not of major interest.
    """
    
    def __init__(self, root, animated_panel,user_os,appearance_mode = "dark"):
        self.user_os = user_os
        self.root = root
        self.animated_panel = animated_panel
        self.animated_panel.do_not_delete = True
        self.animated_panel.do_not_delete2 = True
        self.B_state = False
        self.appearance_mode = appearance_mode 
        self.base_fg_color = self.root.cget("fg_color")
        self.base_animPanel_fg_color = self.animated_panel.cget("fg_color")

        self.NH = hierarchy(0,["./data/hierarchies.csv"]) # Normal hierarchy
        self.IH = hierarchy(1,["./data/hierarchies.csv"])  # Inverted hierarchy
        self.year = 365.25

        # Information about the cores : Name, Power [GW], Baseline [km]
        self.Cores = Dic_of_cores()

        # There are 4 species of interest : U(235), U(238), Pu(239) and Pu(241), representing 99.7% of the thermal power
        # Information about the isotopes : Name, Energy released per fission (MeV, Kopeikin), fission rate (/s, from [358], [arXiv:1101.2266v2]), Fraction of each fuel component
        # The neutrino flux per energy per fission for an isotope can be approximated as the exponential of a quartic polynomial
        self.Isotopes = Dic_of_isotopes()

        self.initialize_ui()

        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("<Return>", self.enter_press)
        self.root.bind('<Control-l>', self.switch_lang_shortcut)
        self.root.bind('<Control-o>', self.open_option_press)
        self.root.bind('<Escape>', self.escape_option_press)
        self.root.bind('<Control-w>', self.on_close)
        self.root.bind('<Control-b>', self.release_B)

    def initialize_ui(self):
        self.language = language
        self.prec_lang = self.language
        self.root.title("TIPP Juno 2025")
        self.create_widgets()  
        self.root.after(100, self.maximize)  
        self.scale_ui(self.root)
    
    def maximize(self):
        self.root.update()
        window_width, window_height = self.root.winfo_reqwidth(), self.root.winfo_reqheight()
        self.screen_width = self.root.winfo_width()
        self.screen_height = self.root.winfo_height()

        self.root.geometry(f"{window_width}x{window_height}+{int((self.screen_width - window_width) / 2)}+{int((self.screen_height - window_height) / 2)}")
        
        if self.user_os == "Windows":
            self.root.state('zoomed')
        elif self.user_os == "Linux":
            self.root.attributes("-zoomed", True)

    def on_close(self,event=None):
    # Safe cleanup for canvas
        if hasattr(self, 'canvas'):
            try:
                self.canvas.get_tk_widget().destroy()
            except Exception as e:
                    print(f"Error destroying canvas: {e}")

    # Check and destroy the loading window if it exists
        if getattr(self, 'loading_window', None) is not None:
            try:
                self.loading_window.destroy()
            except Exception as e:
                print(f"Error destroying loading window: {e}")

    # Check and destroy the animated panel if it exists
        if getattr(self, 'animated_panel', None) is not None:
            try:
                self.animated_panel.destroy()
            except Exception as e:
                print(f"Error destroying animated panel: {e}")

    # Quit the main loop before destroying the root window
        try:
            self.root.quit()
            self.root.destroy()
        except Exception as e:
            print(f"Error shutting down application: {e}")

    def tooltip(self):
        CTkToolTip(self.animated_panel_button, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : CTRL+O or ECHAP")
        CTkToolTip(self.quit_button, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : CTRL+W")
        CTkToolTip(self.Bmode_button, bg_color = "yellow",border_color="lightgreen",border_width=2, message=dico_Textes[self.language][0]+" : CTRL+B")
        CTkToolTip(self.language_dropdown, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : CTRL+L")
        
    def init_colors(self):
        self.hover_color = "gray28" if self.appearance_mode=="Dark" else "#009440" if self.B_state else "white"
        self.text_color = "white" if self.appearance_mode=="Dark" else "#302681" if self.B_state else "black"
        self.label_text_color = "white" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "black"
        self.fg_color_trans="#ffcb00" if self.B_state else "transparent" 
        self.border_color="gold" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "gray"
        self.fg_color_dropdown = "gray28" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "#f0f0f0"
        self.button_color_dropdown = "goldenrod" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white"
        self.button_hover_color_dropdown = "dark goldenrod" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white"
        self.dropdown_hover_color = "gray25" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white"
        self.dropdown_fg_color = "gray28" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white"

    def create_widgets(self,reset=False):
        self.init_colors()
        self.create_images()
        self.create_labels()
        self.create_dropdowns()
        self.create_buttons(reset=reset)
        self.font_label()
        self.tooltip()

    def scale_ui(self, root):
        self.screen_width = root.winfo_screenwidth()
        self.screen_height = root.winfo_screenheight()

        self.full_screen_width = self.screen_width
        self.full_screen_height = self.screen_height

        if self.screen_width != 1920:
          self.scale_factor_width = self.screen_width / (1920 / 1.25)  
          self.scale_factor_height = self.screen_height / (1080 / 1.25)

        else :
            self.scale_factor_width = self.screen_width / (1920 * 1.25)  # Suppose that 1920x1080 à 125% is standard conception
            self.scale_factor_height = self.screen_height / (1080 * 1.25)

    def create_images(self):
        path_assets = ".\\assets\\logo_IPHC.png"
        path_juno = ".\\assets\\juno.png"
        path_Bmode = ".\\assets\\Bmode.png"
        #self.appearance_mode = "Light" #ctk.get_appearance_mode()

        path_opt_modified = ".\\assets\\opt-modified.png" if self.appearance_mode=="Dark" else ".\\assets\\opt-modified_black.png"
        
        if self.user_os == "Linux":
            path_assets = JI.convert_path_linux(path_assets)
            path_juno = JI.convert_path_linux(path_juno)
            path_opt_modified = JI.convert_path_linux(path_opt_modified)
            path_Bmode = JI.convert_path_linux(path_Bmode)

        self.img_iphc = LazyImageLoader(path_assets)
        self.img_juno = LazyImageLoader(path_juno)
        self.opt = LazyImageLoader(path_opt_modified)
        self.ImBmode = LazyImageLoader(path_Bmode)

    def create_labels(self):
        self.image_label_juno = ctk.CTkLabel(
            self.root, text="", image=self.img_juno.image)
        self.image_label_juno.place(x=10, y=10)
        self.image_label_juno.do_not_delete = True

        self.image_label_iphc = ctk.CTkLabel(
            self.root, text="", image=self.img_iphc.image)
        self.image_label_iphc.place(x=150, y=10)
        self.image_label_iphc.do_not_delete = True

        self.authors_label = self.most_labels("Simon BEAUDOIN - Pierre LOUIS--CISTAC - Nathan LEBAS",1100,765)
        self.authors_label.do_not_delete = True
        self.authors_label.do_not_delete2 = True

    def create_dropdowns(self): 
        # Creates and places dropdown menus on the UI for selecting material, data, and test options.
        self.language_dropdown = ctk.CTkOptionMenu(self.animated_panel,
                                                   values= ["English","Francais","Português"],#,"Deutsch","Italiano","Espanol","Chinese"], 
                                                   command = self.change_lang,
                                                   fg_color= "gray28" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "#f0f0f0",
                                                   button_color= "goldenrod" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white",
                                                   button_hover_color= "dark goldenrod" if self.appearance_mode=="Dark" else "#009440" if self.B_state else "white",
                                                   dropdown_hover_color= "gray25" if self.appearance_mode=="Dark" else "#009440" if self.B_state else "white",
                                                   dropdown_fg_color = self.dropdown_fg_color,
                                                   text_color="white" if self.appearance_mode=="Dark" else "#302681" if self.B_state else "black",
                                                   )
        self.language_dropdown.place(x=10,y=340)

    def most_buttons(self,cmd,txt,x,y):
        button = ctk.CTkButton(self.root, 
                               text=txt,
                               command= cmd,
                               fg_color=self.fg_color_trans,
                               hover_color= self.hover_color, 
                               text_color=self.text_color,
                               border_color=self.border_color,  
                               border_width=2,
                               )
        button.place(x=x, y=y)
        return button

    def most_dropdowns(self,vals,x,y):
        dropdown = ctk.CTkOptionMenu(self.root,
                                     values = vals, 
                                     fg_color = self.fg_color_dropdown,
                                     button_color = self.button_color_dropdown,
                                     button_hover_color = self.button_hover_color_dropdown,
                                     dropdown_hover_color = self.dropdown_hover_color,
                                     dropdown_fg_color = self.dropdown_fg_color,
                                     text_color = self.text_color,
                                     )
        dropdown.place(x=x, y=y)
        return dropdown

    def most_labels(self,txt,x,y):
        label = ctk.CTkLabel(self.root,
                             text=txt,
                             fg_color="transparent",
                             text_color = self.label_text_color,
                             )
        label.place(x=x, y=y)
        return label

    def labels_special_fonts(self,txt,x,y):
        label = ctk.CTkLabel(self.root,
                             text=txt,
                             fg_color="transparent",
                             text_color = self.label_text_color,
                             font=("Arial",16),
                             width=200,
                             anchor="w",
                             justify="left"
                             )
        label.place(x=x, y=y)
        return label

    def create_buttons(self, reset=False):
        # Creates and places buttons on the UI for various actions like opening files and analyzing data.
        # Quit button 
        if not hasattr(self,"quit_button") or reset:
            self.quit_button = self.most_buttons(self.on_close,
                                                     dico_Textes[self.language][11],
                                                     0,600)
            self.quit_button.do_not_delete = True
        
            # Button for displaying options animated panel 
            self.animated_panel_button = ctk.CTkButton(self.root,
                                                       image=self.opt.image,
                                                       text="Options",
                                                       command=self.animated_panel.animate,
                                                       fg_color=self.fg_color_trans,
                                                       hover_color= "black" if self.appearance_mode=="Dark" else "white",
                                                       text_color= self.text_color,
                                                       width=30)
            self.animated_panel_button.place(x=0, y=500)
            self.animated_panel_button.do_not_delete = True

            # Button to open the readme file
            self.open_readme_button = ctk.CTkButton(self.animated_panel,
                                                    text="Open README",
                                                    command=self.open_readme_window, 
                                                    fg_color=self.fg_color_trans,
                                                    hover_color= self.hover_color, 
                                                    text_color=self.text_color,
                                                    border_color=self.border_color,
                                                    border_width=2,)
            self.open_readme_button.pack()
            self.open_readme_button.do_not_delete = True
            self.open_readme_button.do_not_delete2 = True

            self.Bmode_button = ctk.CTkButton(self.animated_panel,
                                            image=self.ImBmode.image, text="B mode.exe",
                                            command=self.release_B,
                                            fg_color=self.fg_color_trans,
                                            hover_color= "yellow",
                                            text_color="darkgreen",
                                            border_color="darkgreen",
                                            width=30
                                            )
            self.Bmode_button.place(x=10, y=390)
            self.Bmode_button.do_not_delete = True
            self.Bmode_button.do_not_delete2 = True

        # Button to plot neutrino oscillations
        self.oscillations_button = self.most_buttons(self.init_oscillations,
                                                     dico_Textes[self.language][10],
                                                     130,100)


        # Button to study flux (nb nu/fission ; flux ; visible flux...)
        self.FluxStudy_button = self.most_buttons(self.launch_FluxStudy,
                                                     dico_Textes[self.language][9],
                                                     130,150)

        # Button to do a single epsilon study
        self.SingleEpsilon_button = self.most_buttons(self.init_SingleEpsilon,
                                                     dico_Textes[self.language][16],
                                                     130,200)

        # Button to study precisions over time for 1 parameter
        self.TimeStudy_button = self.most_buttons(self.init_TimeStudy,
                                                     dico_Textes[self.language][7],
                                                     130,250)

        # Button to study precisions for all parameters
        self.ParamStudy_button = self.most_buttons(self.init_ParamStudy,
                                                     dico_Textes[self.language][8],
                                                     130,300)
        
        # Button to study 2 parameters at once 
        self.MultiParamStudy_button = self.most_buttons(self.init_MultiParamStudy,
                                                     dico_Textes[self.language][31],
                                                     130,350)

        # Button to compare hierachies
        self.HieraDecision_button = self.most_buttons(self.init_HieraDecision,
                                                     dico_Textes[self.language][35],
                                                     130,400)
        
        # Button to show hierachies influence on beta beta invariant mass
        self.betabeta_button = self.most_buttons(self.init_betabeta,
                                                     dico_Textes[self.language][42],
                                                     130,450)

    def change_lang(self,Event = None):
        self.prec_lang = self.language
        self.language = self.language_dropdown.get()
        self.reset()

    def change_lang_custom(self,lang):
        self.prec_lang = self.language
        self.language = lang
        self.reset()

    def switch_lang_shortcut(self,event=None):
        if event.keysym == 'l' and event.state & 0x4 and not (event.state & 0x1) :
            current_value = self.language 
            current_index = ["Francais","English","Português"].index(current_value)#"Deutsch","Italiano","Espanol","Chinese"].index(current_value)
            next_index = (current_index + 1) % len(["Francais","English","Português"])#"Deutsch","Italiano","Espanol","Chinese"])
            self.language_dropdown.set(["Francais","English","Português",][next_index])#"Deutsch","Italiano","Espanol","Chinese"][next_index])
            self.change_lang()

    def clear_frame(self):
        for widget in self.root.winfo_children():
            if not hasattr(widget, 'do_not_delete'):
                widget.destroy()

    def back_to_title_screen(self):
        self.clear_frame()
        self.create_buttons()
        self.oscillations_analysis = None
        self.SingleEpsilon_analysis = None
        self.TimeStudy_analysis = None
        self.ParamStudy_analysis = None
        self.MultiParamStudy_analysis = None
        self.HieraDecision_analysis = None
        
        if self.B_state:
            self.call_flag()

    def back_to_title_screen_button(self):
        # Back to main screen button 
        self.back_main = self.HieraDecision_button = self.most_buttons(self.back_to_title_screen,
                                                                       dico_Textes[self.language][1],
                                                                       0,550)

    def hierarchy_checkbox(self):
        self.hiera_checkbox_widg = ctk.CTkCheckBox(self.root, text =dico_Textes[self.language][15],
                                                     fg_color="goldenrod" if self.appearance_mode=="Dark" else "#ffcb00" if self.B_state else "white",
                                                     hover_color = self.hover_color, 
                                                     text_color = self.label_text_color,
                                                     border_color="gold" if self.appearance_mode=="Dark" else "gray", 
                                                     border_width=2,)
        self.hiera_checkbox_widg.place(x=129,y=200)

    def param_dropdown(self):
        # Param dropdown
        self.param_dropdown_widg = self.most_dropdowns(["s12","s13","m21","m31_NH","m13_IH"],130,150)
        self.param_label = self.most_labels(dico_Textes[self.language][3],50,150)

    def time_dropdown(self):
        # Time dropdown
        self.time_dropdown_widg = self.most_dropdowns(["6 years","100 days","20 years"],130,250)
        self.time_label = self.most_labels(dico_Textes[self.language][5],50,250)

    def nbPoints_entry(self,default_val=20):
        # Nb points entry (fine study)
        self.nbPoints_entry_widg = ctk.CTkEntry(
           self.root,
           fg_color=self.fg_color_trans,
           text_color = self.text_color,
           border_color = self.border_color,
           border_width=2,
           )
        self.nbPoints_entry_widg.place(x=129, y=300)
        # NbPoints label
        self.nbPoints_label = self.most_labels(dico_Textes[self.language][4]+f" ({default_val})",10,300)
        
    def plot_fig(self,fig,x,y):
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=self.root)
        widget = canvas.get_tk_widget()
        widget.place(x=x, y=y)
        canvas.draw()

    def SingleEpsilon_widgets(self):
        self.param_dropdown()
        self.hierarchy_checkbox()
        self.time_dropdown()
        self.nbPoints_entry()

        # Start SingleEpsilon Analysis Button
        self.SingleEpsilon_analysis = self.most_buttons(self.launch_SingleEpsilon,
                                                     dico_Textes[self.language][2],
                                                     150,550)
        CTkToolTip(self.SingleEpsilon_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def TimeStudy_widgets(self):
        self.param_dropdown()
        self.hierarchy_checkbox()
        self.nbPoints_entry()

        # Start SingleEpsilon Analysis Button
        self.TimeStudy_analysis = self.most_buttons(self.launch_TimeStudy,
                                                     dico_Textes[self.language][6],
                                                     150,550)
        CTkToolTip(self.TimeStudy_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def ParamStudy_widgets(self):
        self.hierarchy_checkbox()
        self.nbPoints_entry()

        # Start All params Analysis Button
        self.ParamStudy_analysis = self.most_buttons(self.launch_ParamStudy,
                                                     dico_Textes[self.language][13],
                                                     150,550)
        CTkToolTip(self.ParamStudy_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def MultiParamStudy_widgets(self):
        self.param_dropdown()
        self.hierarchy_checkbox()
        self.time_dropdown()
        self.nbPoints_entry(default_val=10)

        # Second Param dropdown
        self.param_dropdown2_widg = self.most_dropdowns(["s12","s13","m21","m31_NH","m13_IH"],380,150)
        self.param2_label = self.most_labels(dico_Textes[self.language][3]+" 2",300,150)

        # Start multiparametric Analysis Button
        self.MultiParamStudy_analysis = self.most_buttons(self.launch_MultiParamStudy,
                                                     dico_Textes[self.language][32],
                                                     150,550)
        CTkToolTip(self.MultiParamStudy_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def HieraDecision_widgets(self):
        # nbseries dropdown
        self.nbSeries_dropdown = self.most_dropdowns(["5","10"],170,150)
        self.nbSeries_label = self.most_labels(dico_Textes[self.language][33],20,150)
        
        # Display existing graphs
        if self.user_os == "Windows":
            self.existingGraph_path = ".\\assets\\saves\\HieraDecision\\"
        elif self.user_os == "Linux":
            self.existingGraph_path = "./assets/saves/HieraDecision/"
        self.existing_df = pd.read_csv(self.existingGraph_path+"HieraDecision.csv")
        if len(self.existing_df)>0:
            self.existingGraphs_dropdown = self.most_dropdowns(["None"]+[str(k) for k in range(len(self.existing_df))],200,200)
            self.existingGraphs_label = self.most_labels(dico_Textes[self.language][36],50,200)

        # Start evolution of chi2
        self.HieraDecision_analysis = self.most_buttons(self.launch_HieraDecision,
                                                     dico_Textes[self.language][34],
                                                     150,550)
        CTkToolTip(self.HieraDecision_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def init_oscillations(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.hierarchy_checkbox()

        # mode dropdown
        self.osc_mode_dropdown_widg = self.most_dropdowns([dico_Textes[self.language][21],dico_Textes[self.language][22]],130,150)
        self.osc_mode_label = self.most_labels(dico_Textes[self.language][20],50,150)

        # init flavor dropdown
        self.init_flavor_dropdown_widg = self.most_dropdowns(["e", "µ", "τ"],130,250)
        self.init_flavor_label = self.most_labels(dico_Textes[self.language][23],50,250)

        # final flavor dropdown
        self.final_flavor_dropdown_widg = self.most_dropdowns(["e", "µ", "τ", "all (individual mode only)"],130,300)
        self.final_flavor_label = self.most_labels(dico_Textes[self.language][24],50,300)

        # Button to start ploting oscillations
        self.oscillations_analysis = self.most_buttons(self.launch_oscillations,
                                                     dico_Textes[self.language][19],
                                                     150,550)
        CTkToolTip(self.oscillations_analysis, bg_color = "#f0f0f0",border_color="gold",border_width=2, message=dico_Textes[self.language][0]+" : ENTER")

    def launch_oscillations(self):
        hiera = self.hiera_checkbox_widg.get()
        mode = self.osc_mode_dropdown_widg.get()
        mode = {dico_Textes[self.language][21]:1,dico_Textes[self.language][22]:2}[mode]

        init_flavor = self.init_flavor_dropdown_widg.get()
        init_flavor = {"e":0, "µ":1, "τ":2}[init_flavor]

        final_flavor = self.final_flavor_dropdown_widg.get()
        final_flavor = {"e":0, "µ":1, "τ":2, "all (individual mode only)":-1}[final_flavor]
        
        print("\n\n\n",mode,hiera,init_flavor,final_flavor)
        
        if mode==2 and final_flavor==-1:
            CTkMessagebox(title=dico_Textes[self.language][25], icon="warning",
                          message=dico_Textes[self.language][26])
            return 

        if not hasattr(self, "osc_interface"):
            self.osc_interface = neutrino_oscillations(self.NH,self.IH)

        fig_oscillations = self.osc_interface.plot_oscillations_GUI(mode,hiera,init_flavor,final_flavor)
        self.plot_fig(fig_oscillations,500,115)

    def launch_FluxStudy(self):
        self.clear_frame()
        self.back_to_title_screen_button()

        if not hasattr(self, "fig_nbnufission"):
            FluxStudy = FluxEstimation(self.NH,self.IH,self.Cores,self.Isotopes,6*self.year)
            self.fig_nbnufission=FluxStudy.plot_Nb_nu_pfission(GUI=True)
            self.fig_RawFlux, self.integral_no_osc,self.integral_NH,self.integral_IH = FluxStudy.plot_flux(needplot=True,GUI=True)
            # Detector response
            self.fig_VisibleFlux, self.i_NH_L, self.i_IH_L, self.i_NH_S, self.i_IH_S = FluxStudy.plot_visible_spectrum(needplot=True,GUI=True)
       
        self.plot_fig(self.fig_nbnufission,2,115)
        size_inches_nb_nufission = self.fig_nbnufission.get_size_inches()
        dpi_nbnufission = self.fig_nbnufission.dpi

        self.plot_fig(self.fig_RawFlux,4+dpi_nbnufission*size_inches_nb_nufission[0],115)
        size_inches_RawFlux = self.fig_nbnufission.get_size_inches()
        dpi_RawFlux = self.fig_nbnufission.dpi

        self.plot_fig(self.fig_VisibleFlux,6+ dpi_nbnufission*size_inches_nb_nufission[0] + dpi_RawFlux*size_inches_RawFlux[0],115)
        self.fig_VisibleFlux.tight_layout()
        
        self.res_flux_label = self.labels_special_fonts(f"Theoretical Flux: \nIntegral without oscillations: {self.integral_no_osc:.4}/day  \nIntegral NH: {self.integral_NH:.3}/day      \nIntegral IH: {self.integral_IH:.3}/day        \n\nMeasured Flux: \nIntegral NH Large PM:  {self.i_NH_L:.3}/day    Integral NH Small PM: {self.i_NH_S:.3}/day       \nIntegral IH Large PM:   {self.i_IH_L:.3}/day    Integral IH Small PM:  {self.i_IH_S:.3}/day       ",
                                                        300,500)

    def launch_betabeta(self):
        if not hasattr(self, "fig_betabeta"):
            betaStudy = DoubleBetaExample()
            self.fig_betabeta = betaStudy.DoubleBetaExample_graph(GUI=True)
        
        self.plot_fig(self.fig_betabeta,500,115)
        self.beta_label = self.labels_special_fonts(f"Example of one measurement whose outcome is directly related to the mass hierarchy.\nDetermining it with JUNO would help constraining the lifetime of neutrinoless double beta decay.",
                                                               500,600)

    def incompability_params(self,param,hiera):
        if ((param=="m31_NH" and hiera==1) or (param=="m13_IH" and hiera==0)):
            msgbox = CTkMessagebox(title=dico_Textes[self.language][37],
                                   option_1=dico_Textes[self.language][40],
                                   icon="warning",
                                   message=dico_Textes[self.language][38]+param+dico_Textes[self.language][39]+str(hiera),
                                   )
            return True

    def insufficient_nbpoints(self,nb_points):
        if nb_points < 5:
            msgbox = CTkMessagebox(title=dico_Textes[self.language][43],
                                   option_1=dico_Textes[self.language][40],
                                   icon="warning",
                                   message= dico_Textes[self.language][44] + str(nb_points) + dico_Textes[self.language][45],
                                   )
            return True

    def check_save_graph(self,study,params):  
        # loads saved data 
        if self.user_os == "Windows":
            path = ".\\assets\\saves\\"+study+"\\"
        elif self.user_os == "Linux":
            path = "./assets/saves/"+study+"/"
        dr = pd.read_csv(path + study+".csv")
        
        dt = dr.copy()
        test = True
        for k,param in enumerate(params):
            col = list(dt.keys())[k]
            if param not in list(dt[col]):
                test = False
            else:
                dt = dt[dt[col]==params[k]].reset_index(drop=True)
        if not test:
            return(test)
        else:
            keys = dr.keys()
            for k in range(len(params)):
                dr = dr[dr[keys[k]]==params[k]].reset_index(drop=True)
            res = dr["res"][0]

            if study in ["TimeStudy","ParamStudy","MultiParam"]:
                res = JI.str2list(res)
                res=res[0][1:]
                res = res[0:len(res)-1]

                with open(path+res+".pkl", "rb") as f:
                    fig = pickle.load(f)
                if study != "MultiParam":
                    return(fig)
                else:
                    return fig, dr["rho"][0] 
            
            elif study in ["SingleEpsilon"]:
                res = JI.str2list(res)
                res1 = res[0][1:]
                res1 = res1[0:len(res1)-1]
                res2 = res[1][1:]
                res2 = res2[1:]
                res2 = res2[0:len(res2)-1]

                with open(path+res1+".pkl", "rb") as f1:
                    fig1 = pickle.load(f1)

                with open(path+res2+".pkl", "rb") as f2:
                    fig2 = pickle.load(f2)

                dep = pd.read_csv(path+"EpsilonVals.csv")
                dep = dep[dep["single_ep_param"]==params[0]].reset_index(drop=True)
                dep = dep[dep["single_ep_hiera"]==params[1]].reset_index(drop=True)
                dep = dep[dep["single_ep_time"]==params[2]].reset_index(drop=True)
                dep = dep[dep["single_ep_nbPoints"]==params[3]].reset_index(drop=True)   
          
                return(dep["ep_mm"][0], dep["ep_mp"][0], dep["ep_m_mean"][0], dep["ep_pm"][0], dep["ep_pp"][0], dep["ep_p_mean"][0],fig1, fig2)

    def save_data(self,study,params,res,rho=None):
        if self.user_os == "Windows":
            path = ".\\assets\\saves\\"+study+"\\"
        elif self.user_os == "Linux":
            path = "./assets/saves/"+study+"/"
        df = pd.read_csv(path + study+".csv")
        
        dict={}
        for k in range(len(params)):
            key = df.keys()[k]
            dict[key] =[params[k] ]

        L = []
        if study in ["TimeStudy","ParamStudy","MultiParam"]:
            x="_".join([str(k) for k in params])
            with open(path+x+".pkl", "wb") as f:
                pickle.dump(res, f)
            L.append(x)
            res = tuple(L)
            dict["res"]=[res]

            if study=="MultiParam":
                dict["rho"] = [rho]
                

        elif study == "SingleEpsilon":
            x1 = "_".join([str(k) for k in params])
            x2 = x1
            x1 += "f1"
            x2 += "f2"
            with open(path+x1+".pkl","wb") as f:
                pickle.dump(res[-2],f)
            with open(path+x2+".pkl","wb") as f:
                pickle.dump(res[-1],f)
            L.append(x1)
            L.append(x2)
            L.append('&')

            dep = pd.read_csv(path+"EpsilonVals.csv")
            new_rowep = pd.DataFrame({"single_ep_param":params[0], "single_ep_hiera":params[1], "single_ep_time":params[2], "single_ep_nbPoints":params[3], "ep_mm":[res[0]], "ep_mp":[res[1]], "ep_m_mean":[res[2]], "ep_pm":[res[3]], "ep_pp":[res[4]], "ep_p_mean":[res[5]]})
            dep = pd.concat([dep, new_rowep], ignore_index=True)
            dep.to_csv(path+"EpsilonVals.csv",index=False)
            res = tuple(L)
            dict["res"]=[res]
        
        elif study == "HieraDecision":
            x1 = "_".join([str(k) for k in params])
            c=0
            while x1+"f1.pkl" in os.listdir(path):
                x1+=str(c)
                c+=1
            x2 = x1
            x1 += "f1"
            x2 += "f2"
            with open(path+x1+".pkl","wb") as f:
                pickle.dump(res[0],f)
            with open(path+x2+".pkl","wb") as f:
                pickle.dump(res[1],f)
            dict["fig1"] = [x1]
            dict["fig2"] = [x2]

        new_row = pd.DataFrame(dict)

        df = pd.concat([df, new_row], ignore_index=True)
        df.to_csv(path+ study+".csv",index=False)

    def init_SingleEpsilon(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.SingleEpsilon_widgets()

    def init_TimeStudy(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.TimeStudy_widgets()

    def init_ParamStudy(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.ParamStudy_widgets()

    def init_MultiParamStudy(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.MultiParamStudy_widgets()

    def init_HieraDecision(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.HieraDecision_widgets()
    
    def init_betabeta(self):
        self.clear_frame()
        self.back_to_title_screen_button()
        self.launch_betabeta()

    def launch_SingleEpsilon(self):
        single_ep_param = self.param_dropdown_widg.get()
        single_ep_hiera = self.hiera_checkbox_widg.get()
        single_ep_time = {"100 days":100,"6 years":6*self.year,"20 years":20*self.year}[self.time_dropdown_widg.get()]
        single_ep_nbPoints = JI.str2nbPoints(self.nbPoints_entry_widg.get())
        if self.incompability_params(single_ep_param,single_ep_hiera):
            return
        if self.insufficient_nbpoints(single_ep_nbPoints):
            return
        print(single_ep_param, single_ep_hiera, single_ep_time, single_ep_nbPoints)

        check_save = self.check_save_graph("SingleEpsilon",[single_ep_param,single_ep_hiera,single_ep_time,single_ep_nbPoints])
        
        if not hasattr(self, 'EpsilonStudies') and (type(check_save)==bool):
            self.EpsilonStudies = PerformanceStudy(self.NH,self.IH,self.Cores,self.Isotopes)
    
        if  type(check_save)==bool:
            self.SingleEpsilon_res = self.EpsilonStudies.full_resolution_computation(single_ep_param,
                                                                                  single_ep_hiera,
                                                                                  single_ep_time,
                                                                                  single_ep_nbPoints,
                                                                                  needplot=False,GUI=True)
            self.save_data("SingleEpsilon",[single_ep_param,single_ep_hiera,single_ep_time,single_ep_nbPoints],self.SingleEpsilon_res)
            self.fig_FluxComparison, self.fig_Chi2Parabolas = self.SingleEpsilon_res[-2],self.SingleEpsilon_res[-1]
        
        elif type(check_save)!=bool:  # if a save exists for the chosen parameters, we don't compute the graph
            self.SingleEpsilon_res = check_save 

        ep_mm, ep_mp, ep_m_mean, ep_pm, ep_pp, ep_p_mean, fig_FluxComparison, fig_Chi2Parabolas = self.SingleEpsilon_res

        self.plot_fig(fig_FluxComparison,500,115)
        size_inches_FluxComparison = fig_FluxComparison.get_size_inches()
        dpi_FluxComparison = fig_FluxComparison.dpi

        self.plot_fig(fig_Chi2Parabolas,500+20+dpi_FluxComparison*size_inches_FluxComparison[0],115)

        self.res_epsilon_label = self.res_flux_label = self.labels_special_fonts(f"Parameter shifts at 1 sigma: lower bound <= estimated uncertainty <= upper bound \n\nepsilon - : {ep_mm:.3e}  <=  {ep_m_mean:.3e}  <=  {ep_mp:.3e}       \nepsilon + :  {ep_pm:.3e}  <=   {ep_p_mean:.3e}  <=   {ep_pp:.3e}        ",
                                                                                 500,600)

    def launch_TimeStudy(self):
        TimeStud_param = self.param_dropdown_widg.get()
        TimeStud_hiera = self.hiera_checkbox_widg.get()
        TimeStud_nbPoints = JI.str2nbPoints(self.nbPoints_entry_widg.get())

        if self.incompability_params(TimeStud_param,TimeStud_hiera):
            return
        if self.insufficient_nbpoints(TimeStud_nbPoints):
            return
        print(TimeStud_param,TimeStud_hiera,TimeStud_nbPoints)
        
        check_save = self.check_save_graph("TimeStudy",[TimeStud_param,TimeStud_hiera,TimeStud_nbPoints])

        if not hasattr(self, 'EpsilonStudies') and (type(check_save)==bool):
            self.EpsilonStudies = PerformanceStudy(self.NH,self.IH,self.Cores,self.Isotopes)

        if type(check_save)==bool:
            self.fig_TimeStudy = self.EpsilonStudies.loop_time(TimeStud_param,
                                                               TimeStud_hiera,
                                                               [100, 6*365.25, 20*365.25],
                                                               TimeStud_nbPoints,
                                                               needplot=False,GUI=True)
            self.save_data("TimeStudy",[TimeStud_param,TimeStud_hiera,TimeStud_nbPoints],self.fig_TimeStudy)

        elif type(check_save)!=bool:
            self.fig_TimeStudy = check_save
        self.plot_fig(self.fig_TimeStudy,500,115)

    def launch_ParamStudy(self):
        ParamStud_hiera = self.hiera_checkbox_widg.get()
        nbstr = self.nbPoints_entry_widg.get()
        ParamStud_nbPoints = JI.str2nbPoints(self.nbPoints_entry_widg.get())
        if self.insufficient_nbpoints(ParamStud_nbPoints):
            return
        print(ParamStud_hiera,ParamStud_nbPoints)

        check_save = self.check_save_graph("ParamStudy",[ParamStud_hiera,ParamStud_nbPoints])

        if not hasattr(self, 'EpsilonStudies') and (type(check_save)==bool):
            self.EpsilonStudies = PerformanceStudy(self.NH,self.IH,self.Cores,self.Isotopes)

        if type(check_save)==bool:
            self.fig_ParamStudy = self.EpsilonStudies.loop_params({0:["s12" , "s13" , "m21", "m31_NH"] , 1:["s12" , "s13" , "m21", "m13_IH"]},
                                                                  ParamStud_hiera,
                                                                  [100, 6*365.25, 20*365.25],
                                                                  ParamStud_nbPoints,
                                                                  needplot=False,GUI=True)
            self.save_data("ParamStudy",[ParamStud_hiera,ParamStud_nbPoints],self.fig_ParamStudy)

        elif type(check_save)!=bool:
            self.fig_ParamStudy = check_save
        self.plot_fig(self.fig_ParamStudy,500,115)


    def launch_MultiParamStudy(self):
        Multi_param1 = self.param_dropdown_widg.get()
        Multi_param2 = self.param_dropdown2_widg.get()
        Multi_hiera = self.hiera_checkbox_widg.get()
        Multi_time = {"100 days":100,"6 years":6*self.year,"20 years":20*self.year}[self.time_dropdown_widg.get()]        
        nbstr = self.nbPoints_entry_widg.get()
        Multi_nbPoints = JI.str2nbPoints(self.nbPoints_entry_widg.get(), default=10)

        if self.incompability_params(Multi_param1,Multi_hiera) or self.incompability_params(Multi_param2,Multi_hiera):
            return
        if self.insufficient_nbpoints(Multi_nbPoints):
            return
        print(Multi_param1, Multi_param2, Multi_hiera, Multi_time, Multi_nbPoints)

        check_save = self.check_save_graph("MultiParam",[Multi_param1,Multi_param2,Multi_hiera,Multi_time,Multi_nbPoints])
        if not hasattr(self, 'EpsilonStudies') and (type(check_save)==bool):
            self.EpsilonStudies = PerformanceStudy(self.NH,self.IH,self.Cores,self.Isotopes)
    
        if  type(check_save)==bool:
            self.fig_MultiParam, self.rho = self.EpsilonStudies.multiparam_study(Multi_param1,
                                                                       Multi_param2,
                                                                       Multi_hiera,
                                                                       Multi_time,
                                                                       nbpoints = Multi_nbPoints,
                                                                       needplot=True,
                                                                       GUI=True)
            self.save_data("MultiParam", [Multi_param1,Multi_param2,Multi_hiera,Multi_time,Multi_nbPoints], self.fig_MultiParam, self.rho)

        elif type(check_save)!=bool:
            self.fig_MultiParam, self.rho = check_save
        self.plot_fig(self.fig_MultiParam,600,115)

        self.rho_label = self.labels_special_fonts(f"Estimated correlation factor: {self.rho:.4}                   ", 500, 600)


    def launch_HieraDecision(self):
        nbSeries = int(self.nbSeries_dropdown.get())
        print(nbSeries)

        if hasattr(self, "existingGraphs_dropdown"):
            graph_to_plot = self.existingGraphs_dropdown.get()
            if graph_to_plot != "None":
                graph_to_plot = int(graph_to_plot)
                x1 = self.existing_df["fig1"][graph_to_plot]
                x2 = self.existing_df["fig2"][graph_to_plot]

                with open(self.existingGraph_path+x1+".pkl","rb") as f1:
                    fig1 = pickle.load(f1)
                with open(self.existingGraph_path+x2+".pkl","rb") as f2:
                    fig2 = pickle.load(f2)
                graph_to_plot = [self.existing_df["nbSeries"][graph_to_plot],self.existing_df["t_3sigma"][graph_to_plot],self.existing_df["t_4sigma"][graph_to_plot],fig1,fig2]
            else:
                graph_to_plot=None
        else:
            graph_to_plot = None

        if not hasattr(self, 'EpsilonStudies') and (type(graph_to_plot)!=bool):
            self.EpsilonStudies = PerformanceStudy(self.NH,self.IH,self.Cores,self.Isotopes)

        if  type(graph_to_plot)!=list:
            self.HieraDecision_res = self.EpsilonStudies.analyze_hierarchy_distinction(100,
                                                                                       10*self.year,
                                                                                       20,
                                                                                       series=nbSeries,
                                                                                       GUI=True)
            
            self.save_data("HieraDecision", [nbSeries, self.HieraDecision_res[0], self.HieraDecision_res[1]], [self.HieraDecision_res[-2],self.HieraDecision_res[-1]])
            t_3sigma, t_4sigma, fig_flux_noise, fig_significance = self.HieraDecision_res

        else:
            self.HieraDecision_res = graph_to_plot
            nbSeries, t_3sigma, t_4sigma, fig_flux_noise, fig_significance = self.HieraDecision_res
        

        self.plot_fig(fig_flux_noise,500,115)
        size_inches_flux_noise = fig_flux_noise.get_size_inches()
        dpi_flux_noise = fig_flux_noise.dpi
        self.plot_fig(fig_significance,500+20+dpi_flux_noise*size_inches_flux_noise[0],115)

        self.res_time_significance = self.labels_special_fonts(f"3-sigma significance expected after at least:  {t_3sigma/self.year:.3} years      \n4-sigma significance expected after at least:  {t_4sigma/self.year:.3} years      ",
                                                                                 500,600)


    def call_flag(self):
        width = 600
        height=600
        self.frame_Bmode = ctk.CTkFrame(self.root, width=width, height=height)
        self.frame_Bmode.place(x=900, y=100)

        self.canvas_Bmode = tCanvas(self.frame_Bmode, width=width, height=height, bg="white")
        self.canvas_Bmode.pack()

        self.turtle_screen = turtle.TurtleScreen(self.canvas_Bmode)
        self.turtle_screen.bgcolor("#ffcb00")
        self.turtle = turtle.RawTurtle(self.turtle_screen)
        JI.B_mode(GUIturtle=self.turtle)

    def release_B(self,event=None):
        if not self.B_state:
            print("B intensifies")
            msgbox = CTkMessagebox(title=dico_Textes[self.language][27],
                        icon="warning",
                        message=dico_Textes[self.language][28],
                        option_1=dico_Textes[self.language][29],
                        option_2=dico_Textes[self.language][30]
                        )
            response = ""
            while response=="":
                response = msgbox.get()
            if response == dico_Textes[self.language][30]:
                return
            self.B_state = True
            self.root.configure(fg_color="#009440")
            self.animated_panel.set_color("#009440")
            self.change_lang_custom("Português")
            self.call_flag()
        
        # get out of B_mode 
        else:
            self.B_state = False
            self.change_lang_custom(self.prec_lang)
            self.root.configure(fg_color = self.base_fg_color)
            self.animated_panel.set_color(self.base_animPanel_fg_color)

    def reset(self):
        for widget in self.root.winfo_children():
            if not hasattr(widget, 'do_not_delete2'):
                widget.destroy()
        for widget in self.animated_panel.winfo_children() :
            widget.destroy()
        self.create_widgets(reset=True)
   
   
    def font_label(self):
        for widget in root.winfo_children():
            if isinstance(widget, ctk.CTkLabel) or isinstance(widget, ctk.CTkButton) or isinstance(widget, ctk.CTkOptionMenu) or  isinstance(widget, ctk.CTkEntry):
                widget.configure(font=("Roboto",16))


    def open_readme_window(self):
        # Function to create a new window and display the contents of readme.md.
        readme_window = ctk.CTkToplevel(self.root)
        readme_window.title("README")
        readme_window.geometry("600x700")  # Adjust size as needed

        # Make the new window modal
        readme_window.grab_set()

        # Read the content of readme.md
        if self.user_os == "Windows":
            with open(".\\readme.md", "r") as file:
                readme_content = file.read()
        elif self.user_os == "Linux":
            with open("./readme.md", "r") as file:
                readme_content = file.read()

        # Display the content in a label or text box
        readme_label = ctk.CTkLabel(readme_window, text=readme_content, wraplength=500)
        readme_label.pack(pady=10, padx=10)


    def enter_press(self,event=None):
        # to launch studies with the enter key
        if hasattr(self, "oscillations_analysis"):
            if self.oscillations_analysis != None:
                self.launch_oscillations()
        if hasattr(self, "SingleEpsilon_analysis"):
            if self.SingleEpsilon_analysis != None:
                self.launch_SingleEpsilon()
        if hasattr(self, "TimeStudy_analysis"):
            if self.TimeStudy_analysis != None:
                self.launch_TimeStudy()
        if hasattr(self, "ParamStudy_analysis"):
            if self.ParamStudy_analysis != None:
                self.launch_ParamStudy()
        if hasattr(self, "MultiParamStudy_analysis"):
            if self.MultiParamStudy_analysis != None:
                self.launch_MultiParamStudy()
        if hasattr(self, "HieraDecision_analysis"):
            if self.HieraDecision_analysis != None:
                self.launch_HieraDecision()

    def open_option_press(self,event=None):
            self.animated_panel.animate()

    def escape_option_press(self,event=None):
        self.animated_panel.animate_backwards()
        
if __name__ == "__main__":
    plt.style.use("fivethirtyeight") #'dark_background' , 'default'
    print(os.getcwd())

    root = ctk.CTk()
    animated_panel = SlidePanel(root, 1.0, 0.7)
    user_os = platform.system()
    appearance_mode = "Light" #Dark , system
    app = JunoApp(root, animated_panel,user_os,appearance_mode)
    ctk.set_appearance_mode(appearance_mode) # dark,system

    if user_os == "Windows":
        root.iconbitmap(r'.\assets\iconJUNO.ico')
    elif user_os == "Linux":
        img = Image.open('./assets/iconJUNO.png').convert("RGBA") 
        icon = ImageTk.PhotoImage(img)  
        root.iconphoto(True, icon)    

    root.mainloop()