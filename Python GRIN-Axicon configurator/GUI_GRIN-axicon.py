import tkinter as tk
from tkinter import ttk
import tkinter.font as font
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib as mpl
from PIL import ImageTk, Image


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
    # Focal distance needed for theta2_ Exit angle (Beta) and 4a0 entrance pupil
    
    return -(180/np.pi)*(a0_/theta2_)

def r2(r0_, n1_, n0_, a0_, z_, f_):
    # radial position of a ray after the lens
    
    return a0_+ f_*theta1(r0_, n1_, n0_, a0_, z_)

def Exit_D(r0_, n1_, Dn_, a0_, z_, f_): 
    # diameter of the exit beam
    
    n0_ = n1_ - Dn_
    return 2*a0_+ 2*f_*theta1(r0_, n1_, n0_, a0_, z_)

def ax_ExitAngle(alpha, n):
    # axicon Exit angle (Beta)
    
    return -(np.arcsin(n*np.sin(alpha*(np.pi/180))) - (np.pi/180)*alpha)*(180/np.pi)

def ax_k_constant(alpha_param):
    # axicon k constant 
    
    return -(1+(np.tan(np.pi/2 - alpha_param*(np.pi/180)))**2)

def r_GRIN_lens(r1_, theta1_, z1_, z_array):
    # create an array with the radial information of a ray through the GRIN
    
    return np.reshape([np.tan(theta1_)*(z_array - z1_)+r1_], [np.size(z_array)])

def r_lens_air(r2_, theta2_, z1_, z_array):
     # create an array with the radial information of a ray after the lens
    
    return np.reshape([np.tan(theta2_)*(z_array - z1_)+r2_], [np.size(z_array)])
   
# ===========================================================
#                       GUI setup
# ===========================================================

class DynamicPlotApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GRIN-Axicon configurator")

        self.fontname = 'Helvetica'

        # focal length init check box val
        self.n0_CheckVar = tk.IntVar(value=1)
        self.Dn_CheckVar = tk.IntVar(value=0)
        self.z_CheckVar = tk.IntVar(value=1)
        self.f_CheckVar = tk.IntVar(value=1)

        self.font = font.Font(family=self.fontname, size=14)
        self.subtitlefont = font.Font(family=self.fontname, size=18, weight='normal', underline=True)
        self.titlefont = font.Font(family=self.fontname, size=20, weight='bold')
        
        self.alpha_ax_val = 2.5
        self.n_ax_val = 1.5168
        self.beam_diam_val = 5.854128
        self.exitangle_ax_val = ax_ExitAngle(self.alpha_ax_val, self.n_ax_val)
        self.ax_k_constant_val = ax_k_constant(self.alpha_ax_val)      
        
        # Default values
        self.default_n0 = 1.45535
        self.default_Dn = 0.01
        self.default_n1 = self.default_n0 + self.default_Dn
        self.default_Z = T(alpha(self.default_n1, self.default_n0, self.beam_diam_val/4))
        
        self.n0 = self.default_n0
        self.Dn = self.default_Dn
        self.z = self.default_Z
        self.f = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
        
        self.grin_lens_distance_val = GRIN_focRing(1E-07, self.default_n1, self.default_n0, self.beam_diam_val/4, self.default_Z) + self.f
        self.grin_ax_exit_pupil_val = Exit_D(1E-07, self.default_n1, self.default_Dn, self.beam_diam_val/4, self.default_Z, self.f)
        self.depth_of_focus_val = self.grin_ax_exit_pupil_val/(2*np.tan(-(np.pi/180)*self.exitangle_ax_val))
        
        self.create_widgets()

    def create_widgets(self):
        # dpi = root.winfo_fpixels('1i') # Found here: https://stackoverflow.com/questions/42961810/detect-dpi-scaling-factor-in-python-tkinter-application
        # screen_width_in_inches = root.winfo_screenwidth() / dpi
        # screen_height_in_inches = root.winfo_screenheight() / dpi

        # mpl.rcParams['axes.titlesize'] = 'x-small'
        # mpl.rcParams['axes.labelsize'] = 'x-small'
        # mpl.rcParams['xtick.labelsize'] = 'x-small'
        # mpl.rcParams['ytick.labelsize'] = 'x-small'

        # self.fig, self.ax = plt.subplots(figsize=(screen_width_in_inches * 0.1, screen_height_in_inches * 0.15) )
        # self.fig1, self.ax1 = plt.subplots(figsize=(screen_width_in_inches * 0.1, screen_height_in_inches * 0.15))
        # self.fig2, self.ax2 = plt.subplots(figsize=(screen_width_in_inches * 0.2, 0.5))

        self.fig, self.ax = plt.subplots(figsize=(10, 4))

        # Create a frame for the GRIN-axicon section
        self.grin_axicon_frame = tk.LabelFrame(self.root, text="GRIN-Axicon", font=self.titlefont)
        self.grin_axicon_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.BOTH, expand=True)
        
        # Create a GRIN-axicon configuration panel
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
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=6, padx=10, pady=20)
        
        self.n0_label = ttk.Label(self.grin_ax_config_panel, text="n0:", font=self.font)
        self.n0_label.grid(row=1, column=0, padx=5, pady=5)
        
        self.n0_label2 = ttk.Label(self.grin_ax_config_panel, text=f"{self.default_n0:.4f}", font=self.font)
        self.n0_label2.grid(row=1, column=1, padx=5, pady=5)
        
        self.n0_check_box = ttk.Checkbutton(self.grin_ax_config_panel, text='default', variable=self.n0_CheckVar, command=self.replace_n0)
        self.n0_check_box.grid(row=1, column=2, padx=5, pady=5, sticky='W')

        self.Dn_label = ttk.Label(self.grin_ax_config_panel, text="Dn:", font=self.font)
        self.Dn_label.grid(row=2, column=0, padx=5, pady=5)
        
        self.Dn_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
        self.Dn_input.insert(0, f"{self.default_Dn:.4f}")
        self.Dn_input.grid(row=2, column=1, padx=5, pady=5)
        self.Dn_input.bind('<Return>', self.update_values)
        
        self.Dn_check_box = ttk.Checkbutton(self.grin_ax_config_panel, text='auto solve', variable=self.Dn_CheckVar)
        self.Dn_check_box.grid(row=2, column=2, padx=5, pady=5, sticky='W')

        self.z_label = ttk.Label(self.grin_ax_config_panel, text="L:", font=self.font)
        self.z_label.grid(row=3, column=0, padx=5, pady=5)
        
        self.z_check_box = ttk.Checkbutton(self.grin_ax_config_panel, text='auto pitch', variable=self.z_CheckVar,  command=self.replace_z)
        self.z_check_box.grid(row=3, column=2, padx=5, pady=5, sticky='W')
        
        self.z_label2 = ttk.Label(self.grin_ax_config_panel, text=f"{self.default_Z:.4f}", font=self.font)
        self.z_label2.grid(row=3, column=1, padx=5, pady=5)
        
        
        
        # Focal length input label, input box and check box
        self.f_input_label = ttk.Label(self.grin_ax_config_panel, text='f:', font=self.font)
        self.f_input_label.grid(row=4, column=0, padx=5, pady=5)
        self.f_input_label2 = ttk.Label(self.grin_ax_config_panel, text=f'{self.f:.4f}', font=self.font)
        self.f_input_label2.grid(row=4, column=1, padx=5, pady=5)

        self.f_check_box = ttk.Checkbutton(self.grin_ax_config_panel, text='auto', variable=self.f_CheckVar, command=self.replace_f)
        self.f_check_box.grid(row=4, column=2, padx=5, pady=5, sticky='W')

        
        

        # Create input boxes in Axicon section
        self.alpha_ax = ttk.Label(self.axicon_config_panel, text=f"Alpha parameter [{chr(176)}]:", font=self.font)
        self.alpha_ax.grid(row=1, column=0, padx=5, pady=5)
        self.alpha_ax = tk.Entry(self.axicon_config_panel, font=self.font)
        self.alpha_ax.insert(0, str(self.alpha_ax_val))
        self.alpha_ax.grid(row=1, column=1, padx=5, pady=5)
        

        self.n_ax = ttk.Label(self.axicon_config_panel, text="Refractive index [-]:", font=self.font)
        self.n_ax.grid(row=2, column=0, padx=5, pady=5)
        self.n_ax = tk.Entry(self.axicon_config_panel, font=self.font)
        self.n_ax.insert(0, str(self.n_ax_val))
        self.n_ax.grid(row=2, column=1, padx=5, pady=5)

        self.beam_diam = ttk.Label(self.axicon_config_panel, text="Beam diameter [mm]:", font=self.font)
        self.beam_diam.grid(row=3, column=0, padx=5, pady=5)
        self.beam_diam = tk.Entry(self.axicon_config_panel, font=self.font)
        self.beam_diam.insert(0, str(self.beam_diam_val))
        self.beam_diam.grid(row=3, column=1, padx=5, pady=5)

        # Create axicon image in configuration panel
        scale_factor = 0.5
        original_width = 893
        original_height = 515
        
        width = round(original_width * scale_factor)
        height = round(original_height * scale_factor)
        
        # Load the image
        img = Image.open("axicon_image.png")

        # Resize the image if needed
        img = img.resize((width, height), Image.LANCZOS)

        # Create a Tkinter-compatible image
        img_tk = ImageTk.PhotoImage(img)
        # Create a label to display the image
        self.image_label = tk.Label(self.axicon_config_panel, image=img_tk)
        self.image_label.grid(row=0, column=0, padx=10, pady=10, columnspan=2)
        self.image_label.image = img_tk  # Keep a reference to prevent garbage collection
                
        

        # Create the axicon properties labels
        self.exitangle_ax = ttk.Label(self.axicon_param_frame, text=f"Exit angle (Beta): {self.exitangle_ax_val:.4f}", font=self.font)
        self.exitangle_ax.grid(row=0, column=0, padx=5, pady=5)

        self.beam_diam_ax = ttk.Label(self.axicon_param_frame, text=f"Beam diameter: {self.beam_diam_val:.4f}", font=self.font)
        self.beam_diam_ax.grid(row=1, column=0, padx=5, pady=5)
        
        self.ax_k_cnt = ttk.Label(self.axicon_param_frame, text=f"K constant: {self.ax_k_constant_val:.4f}", font=self.font)
        self.ax_k_cnt.grid(row=2, column=0, padx=5, pady=5)
        
        # Create the grin axicon properties labels
        self.exitangle_grin_ax = ttk.Label(self.grin_ax_param_frame, text=f"n1: {self.default_n1:.4f}", font=self.font)
        self.exitangle_grin_ax.grid(row=0, column=0, padx=5, pady=5)
        
        self.n0_grin_ax = ttk.Label(self.grin_ax_param_frame, text=f"n0: {self.default_n0:.4f}", font=self.font)
        self.n0_grin_ax.grid(row=0, column=1, padx=5, pady=5)
        
        self.grin_length = ttk.Label(self.grin_ax_param_frame, text=f"GRIN length [mm]: {self.default_Z:.4f}", font=self.font)
        self.grin_length.grid(row=0, column=2, padx=5, pady=5)
        
        self.lens_focal = ttk.Label(self.grin_ax_param_frame, text=f"Lens' focal distance [mm]: {self.f:.4f}", font=self.font)
        self.lens_focal.grid(row=0, column=3, padx=5, pady=5)
        
        self.grin_ax_entrance_pupil = ttk.Label(self.grin_ax_param_frame, text=f"Input beam diameter [mm]: {self.beam_diam_val:.4f}", font=self.font)
        self.grin_ax_entrance_pupil.grid(row=1, column=0, padx=5, pady=5)
        
        self.grin_ax_exit_pupil = ttk.Label(self.grin_ax_param_frame, text=f"Output beam diameter [mm]: {self.grin_ax_exit_pupil_val:.4f}", font=self.font)
        self.grin_ax_exit_pupil.grid(row=1, column=1, padx=5, pady=5)
        
        self.grin_lens_distance = ttk.Label(self.grin_ax_param_frame, text=f"GRIN-lens distance [mm]: {self.grin_lens_distance_val:.4f}", font=self.font)
        self.grin_lens_distance.grid(row=1, column=2, padx=5, pady=5)
        
        self.depth_of_focus = ttk.Label(self.grin_ax_param_frame, text=f"Depth of focus (DOF) [mm]: {self.grin_ax_exit_pupil_val:.4f}", font=self.font)
        self.depth_of_focus.grid(row=1, column=3, padx=5, pady=5)

        self.alpha_ax.bind('<Return>', self.update_ax_param)
        self.n_ax.bind('<Return>', self.update_ax_param)
        self.beam_diam.bind('<Return>', self.update_ax_param)

        # Draw initial plot
        self.update_all(None)






    def replace_n0(self):
        status = self.n0_CheckVar.get()
        
        if status == 0:
            self.n0_label2.destroy()
            
            self.n0_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
            self.n0_input.insert(0, f"{self.default_n0:.4f}")
            self.n0_input.grid(row=1, column=1, padx=5, pady=5)                
            self.n0_input.bind('<Return>', self.update_values)
            
        else:
            self.n0_input.destroy()
            
            self.n0_label2 = ttk.Label(self.grin_ax_config_panel, text=f"{self.default_n0:.4f}", font=self.font)
            self.n0_label2.grid(row=1, column=1, padx=5, pady=5)
            
        self.update_all(None)
        
    def replace_z(self):
        status = self.z_CheckVar.get()
        
        if status == 0:
            self.z_label2.destroy()
            
            self.z_input = tk.Entry(self.grin_ax_config_panel, font=self.font)
            self.z_input.insert(0, f"{self.z:.4f}")
            self.z_input.grid(row=3, column=1, padx=5, pady=5)              
            self.z_input.bind('<Return>', self.update_values)
            
        else:
            self.z = T(alpha(self.n0 + self.Dn, self.n0, self.beam_diam_val/4))
            
            self.z_input.destroy()
            
            self.z_label2 = ttk.Label(self.grin_ax_config_panel, text=f'{self.z:.4f}', font=self.font)
            self.z_label2.grid(row=3, column=1, padx=5, pady=5)
            
        self.update_all(None)
        
    def replace_f(self):
        status = self.f_CheckVar.get()
        
        if status == 0:
            self.f_input_label2.destroy()
            
            self.f_input = ttk.Entry(self.grin_ax_config_panel, font=self.font)
            self.f_input.insert(0, f'{self.f:.4f}')
            self.f_input.grid(row=4, column=1, padx=5, pady=5)              
            self.f_input.bind('<Return>', self.update_values)
            
        else:
            self.f = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
            
            self.f_input.destroy()
            
            self.f_input_label2 = ttk.Label(self.grin_ax_config_panel, text=f'{self.f:.4f}', font=self.font)
            self.f_input_label2.grid(row=4, column=1, padx=5, pady=5)
            
        self.update_all(None)
        
        
    def update_values(self, event):      
        
        if self.n0_CheckVar.get() == 0:
            self.n0 = float(self.n0_input.get())
        
        if self.Dn_CheckVar.get() == 0:
            self.Dn = float(self.Dn_input.get())         
        
        
        
        if self.f_CheckVar.get() == 0:
            self.f = float(self.f_input.get())
        
        else:
            self.f = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
            self.f_input_label2['text'] = f'{self.f:.4f}'
            
        
        
        if self.z_CheckVar.get() == 0:
            self.z = float(self.z_input.get())
        
        else:
            self.z = T(alpha(self.n0 + self.Dn, self.n0, self.beam_diam_val/4))
            self.z_label2['text'] = f'{self.z:.4f}'
        
        
        self.update_all(None)
       
        
        
        
        
    def update_ax_param(self, event):
        self.alpha_ax_val = float(self.alpha_ax.get())
        self.n_ax_val = float(self.n_ax.get())
        self.beam_diam_val = float(self.beam_diam.get())
        self.exitangle_ax_val = ax_ExitAngle(self.alpha_ax_val, self.n_ax_val)
        self.f = f_ga(self.exitangle_ax_val, self.beam_diam_val/4)
        
        self.exitangle_ax['text'] = f"Exit angle (Beta) [{chr(176)}]: {self.exitangle_ax_val:.4f}"
        self.beam_diam_ax['text'] = f"Beam diameter [mm]: {self.beam_diam_val:.4f}"
    
        self.ax_k_cnt['text'] = f"K constant: {ax_k_constant(self.alpha_ax_val):.4f}"
        
        self.update_values(None)

        
        
        
        
        
    def update_grin_ax_param(self, event, n1, n0, Dn, D, z_cursor):

        self.exitangle_grin_ax['text'] = f"n1: {self.default_n0+self.default_Dn:.4f}"
        self.n0_grin_ax['text'] = f"n0: {self.default_n0:.4f}"
        self.grin_length['text'] = f"GRIN length: {self.default_Z:.4f}"
        self.lens_focal['text'] = f"Lens' focal distance [mm]: {self.f:.4f}"
        
        self.grin_lens_distance_val = GRIN_focRing(1E-07, n1, n0, D/4, z_cursor) + self.f
        self.grin_lens_distance['text'] = f"GRIN-lens distance [mm]: {self.grin_lens_distance_val:.4f}"
        
        self.grin_ax_entrance_pupil['text'] = f"Input beam diameter [mm]: {self.beam_diam_val:.4f}"
        self.grin_ax_exit_pupil_val = Exit_D(1E-07, n1, Dn, D/4, z_cursor, self.f)
        self.grin_ax_exit_pupil['text'] = f"Output beam diameter [mm]: {self.grin_ax_exit_pupil_val:.4f}"
        
        self.depth_of_focus_val = self.grin_ax_exit_pupil_val/(2*np.tan(-(np.pi/180)*self.exitangle_ax_val))
        self.depth_of_focus['text'] = f"Depth of focus (DOF) [mm]: {self.depth_of_focus_val:.4f}"
        
    # def index_profile_viewer(self, event, n1_, n0_, a0_):
        
    #     self.ax1.clear()
    #     ray_color = (0.18, 0.45, 0.71)
    #     alpha_ = alpha(n1_, n0_, a0_)
        
    #     r = np.linspace(-2*a0_, 2*a0_, 50)
        
    #     nr_ = []
    #     for i in r:
    #         if i >= 0:
    #             nr_.append(n1_/np.cosh(alpha_*(i - a0_)))
    #         else:
    #             nr_.append(n1_/np.cosh(alpha_*(i + a0_)))

    #     self.ax1.plot(r, nr_, color=ray_color)
    #     self.ax1.set_ylabel('Refractive index $n$ [-]')
    #     self.ax1.set_xlabel('Radial position $r$ [mm]')
    #     self.ax1.set_title('GRIN index profile')
    #     plt.tight_layout()
    #     self.canvas1.draw()
    
    def grin_ax_viewer(self, event, n1_, n0_, a0_, z_, f_, nb_r0, nb_points):
    
        L_grin_lens = GRIN_focRing(1E-07, n1_, n0_, a0_, z_)
        exit_diam = self.grin_ax_exit_pupil_val
        
        r0_ = np.linspace(2*a0_, -2*a0_, nb_r0)
        l1 = z_
        l2 = L_grin_lens + f_ + l1
        l3 = self.depth_of_focus_val + l2
        z1 = np.linspace(0, l1, int(nb_points/2))
        z2 = np.linspace(l1, l2, int(nb_points/4))
        z3 = np.linspace(l2, l3, int(nb_points/4))
        full_z = []
        full_z.extend(z1)
        full_z.extend(z2)
        full_z.extend(z3)
        
        
        self.ax.clear()
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
            
                self.ax.plot(full_z, beam_list, color=ray_color)
                self.ax.plot(full_z, [ -x for x in beam_list], color=ray_color)
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
                
                self.ax.plot(full_z, beam_list, color=ray_color)
                self.ax.plot(full_z, [ -x for x in beam_list], color=ray_color)
            else:
                break
        self.ax.set_xlabel('Axial position [mm]')
        self.ax.set_ylabel('Radial position [mm]')
        self.ax.set_title('GRIN-axicon setup scheme')
        self.ax.vlines(x = [0, z_, l2], ymin = [-2*a0_, -2*a0_, -exit_diam/2], ymax = [2*a0_, 2*a0_, exit_diam/2],
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='solid')
        self.ax.hlines([-2*a0_, 2*a0_], 0, z_,
           colors = 'black',
           label = 'vline_multiple - full height',
           linestyles='solid')
        self.canvas.draw()
        
    def update_all(self, event):        
        n0 = self.n0
        Dn = self.Dn
        n1 = n0 + Dn
        z = self.z
        
        D = self.beam_diam_val
        a0 = D/4
        f = self.f
        
                
        # Update the grin axicon setup values
        self.update_grin_ax_param(event, n1, n0, Dn, D, self.z)   
        
        # Update plot
        self.grin_ax_viewer(self, n1, n0, a0, z, f, 7, 200)
        


if __name__ == "__main__":
    root = tk.Tk()
    app = DynamicPlotApp(root)
    #Get the current screen width and height
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    #Print the screen size
    print("Screen width:", screen_width)
    print("Screen height:", screen_height)
    
    root.geometry(f"{screen_width}x{screen_height}")
    
    root.mainloop()


# add focal rinf distance from GRIN
