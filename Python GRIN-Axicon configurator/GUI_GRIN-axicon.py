import tkinter as tk
from tkinter import ttk
import tkinter.font as font
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib as mpl

# ===========================================================
#                   GRIN-Axicon functions
# ===========================================================

import numpy as np
from numpy import arange
from numpy import meshgrid

from scipy.integrate import odeint
from matplotlib import pyplot as plt

from typing_extensions import dataclass_transform

import csv
import shutil
import os

import pandas as pd

from scipy.integrate import odeint
from matplotlib import pyplot as plt

from typing_extensions import dataclass_transform

from operator import pos
from scipy.optimize import curve_fit

def T(alpha_):
    # GRIN length
    
    return np.pi/(alpha_*2)

def alpha(n1_, n0_, a0_):
    # GRIN convergence
    
    return np.arccosh(n1_/n0_)/a0_

def r1(r0_, n1_, n0_, a0_, z_):
    # radial position of a ray at the end of the GRIN 
    
    alpha_ = alpha(n1_, n0_, a0_)
    if type(z_) == np.array:
        return np.reshape([(np.arcsinh(np.sinh(alpha_*(r0_-a0_))*np.cos(alpha_*z_)))/alpha_ + a0_], [np.size(z_)]) 
    
    else:
        return (np.arcsinh(np.sinh(alpha_*(r0_-a0_))*np.cos(alpha_*z_)))/alpha_ + a0_

def traceDeRayons(r0_, n1_, n0_, a0_, z_, nb_points):
    alpha_ = alpha(n1_, n0_, a0_)
    return 1/alpha_*np.arcsinh(np.sinh(alpha_*(r0_ - a0_))*np.cos(alpha_*z_))+a0_

def nr(r0_, n1_, n0_, a0_, z_):
    # refractive index with respect to the radial position (index profile)
    
    r_ = r1(r0_, n1_, n0_, a0_, z_)
    alpha_ = alpha(n1_, n0_, a0_)
    return n1_/np.cosh(alpha_*(r_ - a0_))

def dr_dz(r0_, n1_, n0_, a0_, z_):
    # slope a ray at the end of the GRIN
    
    alpha_ = alpha(n1_, n0_, a0_)
    return -(np.sin(alpha_*(r0_-a0_))*np.sin(alpha_*z_))/np.cosh(alpha_*(r1(r0_, n1_, n0_, a0_, z_) - a0_))

def snell(n_i, theta_i, n_j = 1):
    # calculate the refracted angle after the GRIN-air interface 
    
    return np.arcsin((n_i/n_j)*np.sin(theta_i))

def theta_int(r0_, n1_, n0_, a0_, z_):
    # angle of a ray before the GRIN-air interface
    
    return np.arctan(dr_dz(r0_, n1_, n0_, a0_, z_))

def theta1(r0_, n1_, n0_, a0_, z_):
    # angle of a ray after the GRIN-air interface
    
    alpha_ = alpha(n1_, n0_, a0_)
    return snell(nr(r0_, n1_, n0_, a0_, z_), theta_int(r0_, n1_, n0_, a0_, z_))

def GRIN_focRing(r0_, n1_, n0_, a0_, z_):
    # GRIN to focal ring distance
    
    return (a0_ - r1(r0_, n1_, n0_, a0_, z_))/theta1(r0_, n1_, n0_, a0_, z_)

def ExitAngle(a0_, f_):
    return -(180/np.pi)*(a0_/f_)

def f_ga(theta2_, a0_):                
    # Focal distance needed for theta2_ exit angle and 4a0 entrance pupil
    
    return -(180/np.pi)*(a0_/theta2_)

def r2(r0_, n1_, n0_, a0_, z_, f_):
    return a0_+ f_*theta1(r0_, n1_, n0_, a0_, z_)

def D2fit(r0_, n1_, Dn_, a0_, z_, f_): 
    # Exit diameter function to fit
    
    n0_ = n1_ - Dn_
    return 2*a0_+ 2*f_*theta1(r0_, n1_, n0_, a0_, z_)

def ax_ExitAngle(alpha, n):
    return -(np.arcsin(n*np.sin(alpha*(np.pi/180))) - (np.pi/180)*alpha)*(180/np.pi)

def ax_k_constant(alpha_param):
    return -(1+(np.tan(np.pi/2 - alpha_param*(np.pi/180)))**2)

def r_GRIN_lens(r1_, theta1_, z1_, z_array):
    return np.reshape([np.tan(theta1_)*(z_array - z1_)+r1_], [np.size(z_array)])

def r_lens_air(r2_, theta2_, z1_, z_array):
    return np.reshape([np.tan(theta2_)*(z_array - z1_)+r2_], [np.size(z_array)])
   
# ===========================================================
#                       GUI setup
# ===========================================================

class DynamicPlotApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GRIN-Axicon configurator")

        self.font = font.Font(family='Helvetica', size=14)
        self.subtitlefont = font.Font(family='Helvetica', size=18, weight='normal', underline=True)
        self.titlefont = font.Font(family='Helvetica', size=20, weight='bold')
        
        self.alpha_ax_val = 2.5
        self.n_ax_val = 1.5168
        self.beam_diam_val = 5.854128
        self.exitangle_ax_val = ax_ExitAngle(self.alpha_ax_val, self.n_ax_val)
        self.ax_k_constant_val = ax_k_constant(self.alpha_ax_val)      
        
        self.min_n0 = 1
        self.max_n0 = 2
        self.init_cursor_n0 = (self.max_n0+self.min_n0)/2
        
        self.min_Dn = 0.001
        self.max_Dn = 0.05
        self.init_cursor_Dn = (self.max_Dn+self.min_Dn)/2
        
        self.init_val_n1 = self.init_cursor_n0 + self.init_cursor_Dn
        
        self.min_z = 0
        self.max_z = T(alpha(self.init_val_n1, self.init_cursor_n0, self.beam_diam_val/4))
        self.init_cursor_z = (self.max_z + self.min_z)/2 
        
        self.f_gax_val = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
        self.grin_lens_distance_val = GRIN_focRing(1E-07, self.init_val_n1, self.init_cursor_n0, self.beam_diam_val/4, self.init_cursor_z) + self.f_gax_val
        
        self.grin_ax_exit_pupil_val = D2fit(1E-07, self.init_val_n1, self.init_cursor_Dn, self.beam_diam_val/4, self.init_cursor_z, self.f_gax_val)
        
        self.create_widgets()

    def create_widgets(self):
        dpi = root.winfo_fpixels('1i') # Found here: https://stackoverflow.com/questions/42961810/detect-dpi-scaling-factor-in-python-tkinter-application
        screen_width_in_inches = root.winfo_screenwidth() / dpi
        screen_height_in_inches = root.winfo_screenheight() / dpi

        mpl.rcParams['axes.titlesize'] = 'x-small'
        mpl.rcParams['axes.labelsize'] = 'x-small'
        mpl.rcParams['xtick.labelsize'] = 'x-small'
        mpl.rcParams['ytick.labelsize'] = 'x-small'

        self.fig, self.ax = plt.subplots(figsize=(screen_width_in_inches * 0.1, screen_height_in_inches * 0.15) )
        self.fig1, self.ax1 = plt.subplots(figsize=(screen_width_in_inches * 0.1, screen_height_in_inches * 0.15))
        self.fig2, self.ax2 = plt.subplots(figsize=(screen_width_in_inches * 0.2, 0.5))


        # Create a frame for the GRIN-axicon section
        self.grin_axicon_frame = tk.LabelFrame(self.root, text="GRIN-Axicon", font=self.titlefont)
        self.grin_axicon_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create athe GRIN-axicon configuration panel
        self.grin_ax_config_panel = tk.LabelFrame(self.grin_axicon_frame, text='Configuration panel', font=self.subtitlefont)
        self.grin_ax_config_panel.pack(side=tk.TOP, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create a frame inside the GRIN-axicon section for the GRIN-axicon properties
        self.grin_ax_param_frame = tk.LabelFrame(self.grin_axicon_frame, text='GRIN-Axicon setup', font=self.subtitlefont)
        self.grin_ax_param_frame.pack(side=tk.BOTTOM, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create a frame for the Axicon section
        self.axicon_frame = tk.LabelFrame(self.root, text="Axicon", font=self.titlefont)
        self.axicon_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create axicon configuration panel
        self.axicon_config_panel = tk.LabelFrame(self.axicon_frame, text="Configuration panel", font=self.subtitlefont)
        self.axicon_config_panel.pack(side=tk.TOP, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create a frame inside the Axicon section for the axicon properties
        self.axicon_param_frame = tk.LabelFrame(self.axicon_frame, text="Axicon properties", font=self.subtitlefont)
        self.axicon_param_frame.pack(side=tk.BOTTOM, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create a matplotlib canvas to display the plot inside the GRIN-axicon section
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.grin_ax_config_panel)
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=3, padx=10, pady=20)
        
        # Create a matplotlib canvas to display the plot inside the GRIN-axicon section
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.grin_ax_config_panel)
        self.canvas1.get_tk_widget().grid(row=0, column=3, columnspan=3, padx=10, pady=20)
        
        # Create a matplotlib canvas to display the plot inside the GRIN-axicon section
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.grin_ax_param_frame)
        self.canvas2.get_tk_widget().grid(row=2, column=0, columnspan=6, padx=10, pady=20)

        # Create tkinter widgets for sliders in GRIN-axicon section
        self.n0_label = ttk.Label(self.grin_ax_config_panel, text="n0:", font=self.font)
        self.n0_label.grid(row=1, column=0, padx=5, pady=5)
        self.n0_slider = ttk.Scale(self.grin_ax_config_panel, from_= self.min_n0, to = self.max_n0, length=300, orient=tk.HORIZONTAL, value=self.init_cursor_n0)
        self.n0_slider.grid(row=1, column=1, padx=5, pady=5)

        self.n0_min_input_label = ttk.Label(self.grin_ax_config_panel, text="Min:", font=self.font)
        self.n0_min_input_label.grid(row=1, column=2, padx=5, pady=5)
        self.n0_min_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
        self.n0_min_input.insert(0, str(self.min_n0))
        self.n0_min_input.grid(row=1, column=3, padx=5, pady=5)

        self.n0_max_input_label = ttk.Label(self.grin_ax_config_panel, text="Max:", font=self.font)
        self.n0_max_input_label.grid(row=1, column=4, padx=5, pady=5)
        self.n0_max_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
        self.n0_max_input.insert(0, str(self.max_n0))
        self.n0_max_input.grid(row=1, column=5, padx=5, pady=5)

        self.Dn_label = ttk.Label(self.grin_ax_config_panel, text="Dn:", font=self.font)
        self.Dn_label.grid(row=2, column=0, padx=5, pady=5)
        self.Dn_slider = ttk.Scale(self.grin_ax_config_panel, from_ = self.min_Dn, to = self.max_Dn, length=300, orient=tk.HORIZONTAL, value=self.init_cursor_Dn)
        self.Dn_slider.grid(row=2, column=1, padx=5, pady=5)

        self.Dn_min_input_label = ttk.Label(self.grin_ax_config_panel, text="Min:", font=self.font)
        self.Dn_min_input_label.grid(row=2, column=2, padx=5, pady=5)
        self.Dn_min_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
        self.Dn_min_input.insert(0, str(self.min_Dn))
        self.Dn_min_input.grid(row=2, column=3, padx=5, pady=5)

        self.Dn_max_input_label = ttk.Label(self.grin_ax_config_panel, text="Max:", font=self.font)
        self.Dn_max_input_label.grid(row=2, column=4, padx=5, pady=5)
        self.Dn_max_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
        self.Dn_max_input.insert(0, str(self.max_Dn))
        self.Dn_max_input.grid(row=2, column=5, padx=5, pady=5)

        self.z_label = ttk.Label(self.grin_ax_config_panel, text="z:", font=self.font)
        self.z_label.grid(row=3, column=0, padx=5, pady=5)
        self.z_slider = ttk.Scale(self.grin_ax_config_panel, from_ = self.min_z, to = self.max_z, length=300, orient=tk.HORIZONTAL, value=self.init_cursor_z)
        self.z_slider.grid(row=3, column=1, padx=5, pady=5)

        self.z_min_input_label = ttk.Label(self.grin_ax_config_panel, text="Min:", font=self.font)
        self.z_min_input_label.grid(row=3, column=2, padx=5, pady=5)
        self.z_min_value_label = ttk.Label(self.grin_ax_config_panel, text=str(self.min_z), font=self.font)
        self.z_min_value_label.grid(row=3, column=3, padx=5, pady=5)

        self.z_max_input_label = ttk.Label(self.grin_ax_config_panel, text="Max:", font=self.font)
        self.z_max_input_label.grid(row=3, column=4, padx=5, pady=5)
        self.z_max_value_label = ttk.Label(self.grin_ax_config_panel, text=str(self.max_z), font=self.font)
        self.z_max_value_label.grid(row=3, column=5, padx=5, pady=5)

        # Create input boxes in Axicon section
        self.alpha_ax = ttk.Label(self.axicon_config_panel, text=f"Alpha parameter [{chr(176)}]:", font=self.font)
        self.alpha_ax.grid(row=0, column=0, padx=5, pady=5)
        self.alpha_ax = tk.Entry(self.axicon_config_panel, font=self.font)
        self.alpha_ax.insert(0, str(self.alpha_ax_val))
        self.alpha_ax.grid(row=0, column=1, padx=5, pady=5)

        self.n_ax = ttk.Label(self.axicon_config_panel, text="Refractive index [-]:", font=self.font)
        self.n_ax.grid(row=1, column=0, padx=5, pady=5)
        self.n_ax = tk.Entry(self.axicon_config_panel, font=self.font)
        self.n_ax.insert(0, str(self.n_ax_val))
        self.n_ax.grid(row=1, column=1, padx=5, pady=5)

        self.beam_diam = ttk.Label(self.axicon_config_panel, text="Beam diameter [mm]:", font=self.font)
        self.beam_diam.grid(row=2, column=0, padx=5, pady=5)
        self.beam_diam = tk.Entry(self.axicon_config_panel, font=self.font)
        self.beam_diam.insert(0, str(self.beam_diam_val))
        self.beam_diam.grid(row=2, column=1, padx=5, pady=5)

        # Create the axicon properties labels
        self.exitangle_ax = ttk.Label(self.axicon_param_frame, text=f"Exit angle: {self.exitangle_ax_val:.5f}", font=self.font)
        self.exitangle_ax.grid(row=0, column=0, padx=5, pady=5)

        self.beam_diam_ax = ttk.Label(self.axicon_param_frame, text=f"Beam diameter: {self.beam_diam_val:.5f}", font=self.font)
        self.beam_diam_ax.grid(row=1, column=0, padx=5, pady=5)
        
        self.ax_k_cnt = ttk.Label(self.axicon_param_frame, text=f"K constant: {self.ax_k_constant_val:.5f}", font=self.font)
        self.ax_k_cnt.grid(row=2, column=0, padx=5, pady=5)
        
        # Create the grin axicon properties labels
        self.exitangle_grin_ax = ttk.Label(self.grin_ax_param_frame, text=f"n1: {self.n0_slider.get()+self.Dn_slider.get():.5f}", font=self.font)
        self.exitangle_grin_ax.grid(row=0, column=0, padx=5, pady=5)
        
        self.n0_grin_ax = ttk.Label(self.grin_ax_param_frame, text=f"n0: {self.n0_slider.get():.5f}", font=self.font)
        self.n0_grin_ax.grid(row=0, column=1, padx=5, pady=5)
        
        self.grin_length = ttk.Label(self.grin_ax_param_frame, text=f"GRIN length [mm]: {self.z_slider.get():.5f}", font=self.font)
        self.grin_length.grid(row=0, column=2, padx=5, pady=5)
        
        self.lens_focal = ttk.Label(self.grin_ax_param_frame, text=f"Lens' focal distance [mm]: {self.f_gax_val:.5f}", font=self.font)
        self.lens_focal.grid(row=0, column=3, padx=5, pady=5)
        
        self.grin_lens_distance = ttk.Label(self.grin_ax_param_frame, text=f"GRIN-lens distance [mm]: {self.grin_lens_distance_val:.5f}", font=self.font)
        self.grin_lens_distance.grid(row=1, column=2, padx=5, pady=5)
        
        self.grin_ax_entrance_pupil = ttk.Label(self.grin_ax_param_frame, text=f"Input beam diameter [mm]: {self.beam_diam_val:.5f}", font=self.font)
        self.grin_ax_entrance_pupil.grid(row=1, column=0, padx=5, pady=5)
        
        self.grin_ax_exit_pupil = ttk.Label(self.grin_ax_param_frame, text=f"Output beam diameter [mm]: {self.grin_ax_exit_pupil_val:.5f}", font=self.font)
        self.grin_ax_exit_pupil.grid(row=1, column=1, padx=5, pady=5)
        
        
        
        


        # Bind the update_plot function to the "<Motion>" event of the sliders in Axicon section
        self.n0_slider.bind("<Motion>", self.update_plot)
        self.Dn_slider.bind("<Motion>", self.update_plot)
        self.z_slider.bind("<Motion>", self.update_plot)
        
        self.alpha_ax.bind('<Return>', self.update_ax_param)
        self.n_ax.bind('<Return>', self.update_ax_param)
        self.beam_diam.bind('<Return>', self.update_ax_param)
        
        self.n0_min_input.bind("<Return>", self.update_slider_range)
        self.n0_max_input.bind("<Return>", self.update_slider_range)
        
        self.Dn_min_input.bind("<Return>", self.update_slider_range)
        self.Dn_max_input.bind("<Return>", self.update_slider_range)
        
        # Draw initial plot
        self.update_plot(None)

    def update_ax_param(self, event):
        self.alpha_ax_val = float(self.alpha_ax.get())
        self.n_ax_val = float(self.n_ax.get())
        self.beam_diam_val = float(self.beam_diam.get())
        self.exitangle_ax_val = ax_ExitAngle(self.alpha_ax_val, self.n_ax_val)
        self.f_gax_val = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
        
        self.exitangle_ax['text'] = f"Exit angle [{chr(176)}]: {self.exitangle_ax_val:.5f}"
        self.beam_diam_ax['text'] = f"Beam diameter [mm]: {self.beam_diam_val:.5f}"
    
        self.ax_k_cnt['text'] = f"K constant: {ax_k_constant(self.alpha_ax_val):.5f}"
        
        self.update_plot(event)
        
    def update_grin_ax_param(self, event, n1, n0, Dn, D, z_cursor):
        self.exitangle_grin_ax['text'] = f"n1: {self.n0_slider.get()+self.Dn_slider.get():.5f}"
        self.n0_grin_ax['text'] = f"n0: {self.n0_slider.get():.5f}"
        self.grin_length['text'] = f"GRIN length: {self.z_slider.get():.5f}"
        self.lens_focal['text'] = f"Lens' focal distance [mm]: {self.f_gax_val:.5f}"
        
        self.grin_lens_distance_val = GRIN_focRing(1E-07, n1, n0, D/4, z_cursor) + self.f_gax_val
        self.grin_lens_distance['text'] = f"GRIN-lens distance [mm]: {self.grin_lens_distance_val:.5f}"
        
        self.grin_ax_entrance_pupil['text'] = f"Input beam diameter [mm]: {self.beam_diam_val:.5f}"
        self.grin_ax_exit_pupil_val = D2fit(1E-07, n1, Dn, D/4, z_cursor, self.f_gax_val)
        self.grin_ax_exit_pupil['text'] = f"Output beam diameter [mm]: {self.grin_ax_exit_pupil_val:.5f}"
        
    def index_profile_viewer(self, event, n1_, n0_, a0_):
        
        self.ax1.clear()
        ray_color = (0.18, 0.45, 0.71)
        alpha_ = alpha(n1_, n0_, a0_)
        
        r = np.linspace(-2*a0_, 2*a0_, 100)
        
        nr_ = []
        for i in r:
            if i >= 0:
                nr_.append(n1_/np.cosh(alpha_*(i - a0_)))
            else:
                nr_.append(n1_/np.cosh(alpha_*(i + a0_)))

        self.ax1.plot(r, nr_, color=ray_color)
        self.ax1.set_ylabel('Refractive index $n$ [-]')
        self.ax1.set_xlabel('Radial position $r$ [mm]')
        self.ax1.set_title('GRIN index profile')
        plt.tight_layout()
        self.canvas1.draw()
        
    
    def grin_ax_viewer(self, event, n1_, n0_, a0_, z_, f_, nb_r0, nb_points):
    
        T_grin = T(alpha(n1_, n0_, a0_))
        L_grin_lens = GRIN_focRing(1E-07, n1_, n0_, a0_, z_)
        exit_diam = self.grin_ax_exit_pupil_val
        
        r0_ = np.linspace(2*a0_, -2*a0_, nb_r0)
        l1 = z_
        l2 = L_grin_lens + f_ + l1
        l3 = exit_diam/(2*np.tan(-(np.pi/180)*self.exitangle_ax_val)) + l2
        z1 = np.linspace(0, l1, int(nb_points/2))
        z2 = np.linspace(l1, l2, int(nb_points/4))
        z3 = np.linspace(l2, l3, int(nb_points/4))
        full_z = []
        full_z.extend(z1)
        full_z.extend(z2)
        full_z.extend(z3)
        
        
        self.ax2.clear()
        ray_color = (0.18, 0.45, 0.71)
        
        for i in r0_:
            if i > 0:
                beam_list = []
                
                r_z1 = r1(i, n1_, n0_, a0_, z1)
                beam_list.extend(r_z1)
                
                r1_ = r1(i, n1_, n0_, a0_, z_)
                theta1_ = theta1(i, n1_, n0_, a0_, z_)
                
                r_z2 = r_GRIN_lens(r1_, theta1_, z_, z2)
                beam_list.extend(r_z2)
                
                r2_ = r2(i, n1_, n0_, a0_, z_, f_)
                theta2_ = (np.pi/180)*ExitAngle(a0_, f_)
                
                r_z3 = r_lens_air(r2_, theta2_, (z_ + L_grin_lens + f_), z3)
                beam_list.extend(r_z3)
            
                self.ax2.plot(full_z, beam_list, color=ray_color)
                self.ax2.plot(full_z, [ -x for x in beam_list], color=ray_color)
            elif i == 0:
                beam_list = []
                
                r_z1 = r1(i, n1_, n0_, a0_, z1)
                beam_list.extend(r_z1)
                
                r1_ = r1(i, n1_, n0_, a0_, z_)
                theta1_ = theta1(i, n1_, n0_, a0_, z_)
                
                r_z2 = r_GRIN_lens(r1_, theta1_, z_, z2)
                beam_list.extend(r_z2)
                
                r2_ = r2(i, n1_, n0_, a0_, z_, f_)
                theta2_ = (np.pi/180)*ExitAngle(a0_, f_)
                
                r_z3 = r_lens_air(r2_, theta2_, (z_ + L_grin_lens + f_), z3)
                beam_list.extend(r_z3)
                
                self.ax2.plot(full_z, beam_list, color=ray_color)
                self.ax2.plot(full_z, [ -x for x in beam_list], color=ray_color)
            else:
                break
        self.ax2.set_xlabel('Axial position [mm]')
        self.ax2.set_ylabel('Radial position [mm]')
        self.ax2.set_title('GRIN-axicon setup scheme')
        self.ax2.vlines(x = [0, z_, l2], ymin = [-2*a0_, -2*a0_, -exit_diam/2], ymax = [2*a0_, 2*a0_, exit_diam/2],
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='solid')
        self.ax2.hlines([-2*a0_, 2*a0_], 0, z_,
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='solid')
        self.canvas2.draw()
    
    # Function to update the plot based on the slider values
    def update_plot(self, event):        
        n0 = self.n0_slider.get()
        Dn = self.Dn_slider.get()
        n1 = n0 + Dn
        D = self.beam_diam_val
        a0 = D/4
        f = self.f_gax_val
                
        # Update the z slider range    
        self.update_z_slider_range(event)    
        
        # Get cursor position
        z_cursor = self.z_slider.get() 
        y_cursor = D2fit(1E-07, n1, Dn, a0, z_cursor, self.f_gax_val)
        
        # Update the grin axicon setup values
        self.update_grin_ax_param(event, n1, n0, Dn, D, z_cursor)   
        
        z = np.linspace(0, self.max_z, 50)
        y = D2fit(1E-07, n1, Dn, a0, z, self.f_gax_val)
        
        self.ax.clear()
        self.ax.plot(z, y)
        
        self.ax.vlines(x = z_cursor, ymin = 0, ymax = y_cursor,
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='dashed')
        
        self.ax.hlines(y_cursor, 0, z_cursor,
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='dashed')
        
        self.ax.set_xlabel('GRIN length [mm]', fontsize="x-small")
        self.ax.set_ylabel('Exit beam diameter [mm]', fontsize="x-small")
        self.ax.set_title('GRIN exit beam diameter adjusting curve', fontsize="x-small")
        self.canvas.draw()

        # Update slider value labels
        self.n0_label.config(text=f"n0: {n0:.5f}")
        self.Dn_label.config(text=f"Dn: {Dn:.5f}")
        self.z_label.config(text=f"z: {z_cursor:.5f}")
        
        self.grin_ax_viewer(event, n1, n0, a0, z_cursor, f, 5, 50)
        self.index_profile_viewer(event, n1, n0, a0)

    # Function to update the slider range based on input box values
    def update_slider_range(self, event):           
                
        self.min_n0 = float(self.n0_min_input.get())
        self.max_n0 = float(self.n0_max_input.get())
        
        self.min_Dn = float(self.Dn_min_input.get())
        self.max_Dn = float(self.Dn_max_input.get())
        
        self.update_slider(self.n0_slider, self.n0_min_input, self.n0_max_input)
        self.update_slider(self.Dn_slider, self.Dn_min_input, self.Dn_max_input)
        
        self.update_plot(event)

    # Function to update the z slider's range based on the n0, Dn and Beam Diameter
    def update_z_slider_range(self, event):
        n0 = self.n0_slider.get()
        Dn = self.Dn_slider.get()
        n1 = n0 + Dn
        D = self.beam_diam_val
        
        self.max_z = T(alpha(n1, n0, D/4))
        self.z_max_value_label['text'] = str(self.max_z)
        
        try:
            if self.min_z >= self.max_z:
                raise ValueError
            self.z_slider.config(from_=self.min_z, to=self.max_z)
            
            if self.z_slider['value'] > self.max_z:
                self.z_slider['value'] = self.max_z
        except ValueError:
            pass




    # Function to update a single slider's range based on input box values
    def update_slider(self, slider, min_input, max_input):
        try:
            min_value = float(min_input.get())
            max_value = float(max_input.get())

            if min_value >= max_value:
                raise ValueError
            slider.config(from_=min_value, to=max_value)
        except ValueError:
            pass


if __name__ == "__main__":
    root = tk.Tk()
    app = DynamicPlotApp(root)
    root.mainloop()


