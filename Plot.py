import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math
import matplotlib.ticker as plticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.patches import Polygon
import matplotlib.ticker as plticker
from matplotlib import rc, rcParams
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rcParams.update({'font.size': 18})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

Input = np.genfromtxt("Input/ip.dat")

if(Input[22]==0):
    label = "E*F(E)"
elif (np.abs(Input[23]) == 1):
    label = "Real Comp"
elif (np.abs(Input[23]) == 2):
    label = "Imaginary Comp"
elif (np.abs(Input[23]) == 3 or np.abs(Input[23]) == 5):
    label = "Modulus"   
elif (np.abs(Input[23]) == 4 or np.abs(Input[23]) == 6 or Input[23] == 7):
    label = "Lag(s)"    
else:
    label = "Unclear"   
  
Total = np.genfromtxt("Output/Total.dat") 
PivotingPL = np.genfromtxt("Output/PivPL.dat")
LightTravelTime = np.genfromtxt("Output/Reverb.dat") 
PivotingReflection = np.genfromtxt("Output/PivRef.dat")
IonVariations = np.genfromtxt("Output/IonVar.dat")
Continuum = np.genfromtxt("Output/Continuum_spec.dat")

if (Input[23] != 7):
    min_x_timing = 0.3
    max_x_timing = 10.

    min_x_spectrum = 0.1
    max_x_spectrum = 100.

    max_y_timing = -1e20
    min_y_timing= 1e20

    max_y_spectrum = -1e20
    min_y_spectrum = 1e20

    #this looks for the maximum/minimum values of the model output of the code in the x-axis interval selected, for the timing modes
    for i in range(len(Total.T[0])):
        if (Total.T[1][i] < min_y_timing and Total.T[0][i] > min_x_timing and Total.T[0][i] < max_x_timing):
            min_y_timing = Total.T[1][i]
        if (Total.T[1][i] > max_y_timing and Total.T[0][i] > min_x_timing and Total.T[0][i] < max_x_timing):
            max_y_timing = Total.T[1][i]    
    ''''
    for i in range(len(LightTravelTime.T[0])):
        if (LightTravelTime.T[1][i] < min_y_timing and LightTravelTime.T[0][i] > min_x_timing and LightTravelTime.T[0][i] < max_x_timing):
            min_y_timing = LightTravelTime .T[1][i]
        if (LightTravelTime.T[1][i] > max_y_timing and LightTravelTime.T[0][i] > min_x_timing and LightTravelTime.T[0][i] < max_x_timing):
            max_y_timing = LightTravelTime.T[1][i]
    for i in range(len(PivotingPL.T[0])):
        if (PivotingPL.T[1][i] < min_y_timing and PivotingPL.T[0][i] > min_x_timing and PivotingPL.T[0][i] < max_x_timing):
            min_y_timing = PivotingPL .T[1][i]
        if (PivotingPL.T[1][i] > max_y_timing and PivotingPL.T[0][i] > min_x_timing and PivotingPL.T[0][i] < max_x_timing):
            max_y_timing = PivotingPL.T[1][i]
    for i in range(len(PivotingReflection.T[0])):
        if (PivotingReflection.T[1][i] < min_y_timing and PivotingReflection.T[0][i] > min_x_timing and PivotingReflection.T[0][i] < max_x_timing):
            min_y_timing = PivotingReflection .T[1][i]
        if (PivotingReflection.T[1][i] > max_y_timing and PivotingReflection.T[0][i] > min_x_timing and PivotingReflection.T[0][i] < max_x_timing):
            max_y_timing = PivotingReflection.T[1][i]
    for i in range(len(IonVariations.T[0])):
        if (IonVariations.T[1][i] < min_y_timing and IonVariations.T[0][i] > min_x_timing and IonVariations.T[0][i] < max_x_timing):
            min_y_timing = IonVariations .T[1][i]
        if (IonVariations.T[1][i] > max_y_timing and IonVariations.T[0][i] > min_x_timing and IonVariations.T[0][i] < max_x_timing):
            max_y_timing = IonVariations.T[1][i]
    '''
    if (min_y_timing < 0):
        min_y_timing = 1.5*min_y_timing
    else:
        min_y_timing = 0.5*min_y_timing
    max_y_timing = 1.2*max_y_timing 
            
    #same but for the spectra:
    for i in range(len(Total.T[0])):
        if (Total.T[1][i]*(Total.T[0][i])**2 < min_y_spectrum and Total.T[0][i] > min_x_spectrum and Total.T[0][i] < max_x_spectrum):
            min_y_spectrum = Total.T[1][i]*(Total.T[0][i])**2
        if (Total.T[1][i]*(Total.T[0][i])**2 > max_y_spectrum and Total.T[0][i] > min_x_spectrum and Total.T[0][i] < max_x_spectrum):
            max_y_spectrum = Total.T[1][i]*(Total.T[0][i])**2      

    min_y_spectrum = min_y_spectrum/1.5
    max_y_spectrum = 1.5*max_y_spectrum 
else:
    if (Input[17] > 1.e4):
        min_x_timing = 0.7e-5
        max_x_timing = 3.e-2
    elif (Input[23] == 7):
        min_x_timing = 0.1
        max_x_timing = 450.        
    else:
        min_x_timing = 0.1
        max_x_timing = 300.

    max_y_timing = -1e20
    min_y_timing= 1e20

    #this looks for the maximum/minimum values of the model output of the code in the x-axis interval selected, for the timing modes
    for i in range(len(Total.T[0])):
        if (Total.T[1][i] < min_y_timing and Total.T[0][i] > min_x_timing and Total.T[0][i] < max_x_timing):
            min_y_timing = Total.T[1][i]
        if (Total.T[1][i] > max_y_timing and Total.T[0][i] > min_x_timing and Total.T[0][i] < max_x_timing):
            max_y_timing = Total.T[1][i]   
    max_y_timing = 1.5*max_y_timing
    if (min_y_timing < 0):
        min_y_timing = 1.5*min_y_timing
    else:
        min_y_timing = 0.9*min_y_timing
    #this is just to make the plot a bit prettier - at low frequencies it's super messy because the continuum lags become huge
    #min_x_timing = 0.9*min_x_timing

fig, (ax1) = plt.subplots(1,1,figsize=(7.5,4.5))

dashed_line = np.zeros(50)
line_array = np.logspace(np.log10(min_x_timing),np.log10(max_x_timing),50)

#timing plot
if(Input[22] != 0 and Input[23] != 7):
    ax1.plot(Total.T[0],Total.T[1],linewidth=3.0,label='Total',color=colors[0],zorder=2)
    ax1.plot(LightTravelTime.T[0],LightTravelTime.T[1],linewidth=2.0,label='Reverberation',color=colors[1],zorder=1)
    ax1.plot(PivotingPL.T[0],PivotingPL.T[1],linewidth=2.0,label='Continuum',color=colors[2],zorder=1)
    ax1.plot(PivotingReflection.T[0],PivotingReflection.T[1],linewidth=2.0,label='Piv. Reflection',color=colors[3],zorder=1)
    ax1.plot(IonVariations.T[0],IonVariations.T[1],linewidth=2.0,label='Ionisation Variations',color=colors[4],zorder=1)
    ax1.plot(line_array,dashed_line,linestyle='dashed',linewidth=1.0,color='black')
    ax1.legend(loc='best',fontsize=14)
    # ax1.set_xlim([min_x_timing,max_x_timing])
    # ax1.set_ylim([min_y_timing,max_y_timing])
    ax1.set_xscale('log', base=10)
    if (Input[23] == 3 or Input[23] == 5):
        ax1.set_yscale('log', base=10)    
    ax1.set_xlabel('Energy (kev)',fontsize=22)
    ax1.set_ylabel(label,fontsize=22)
#lag-frequency plot
elif(Input[23] == 7):
    ax1.plot(Total.T[0],Total.T[1],linewidth=2.5,label='Total',color=colors[1],zorder=2)    
    ax1.plot(line_array,dashed_line,linestyle='dashed',linewidth=1.0,color='black')
    #ax1.legend(loc='best',fontsize=14)
    #ax1.set_xlim([min_x_timing,max_x_timing])
    #ax1.set_ylim([min_y_timing,max_y_timing])
    ax1.set_xscale('log', base=10)
    ax1.set_xlabel('Frequency (Hz)',fontsize=22)
    ax1.set_ylabel(label,fontsize=22)
#spectrum plot
else:
    ax1.plot(Total.T[0],Total.T[1]*(Total.T[0])**2,linewidth=2.5,label='Total')
    ax1.plot(Continuum.T[0],Continuum.T[1]*(Continuum.T[0])**2,linewidth=2.5,label='Continuum')
    ax1.set_yscale('log', base=10)
    ax1.set_xlim([min_x_spectrum,max_x_spectrum])
    ax1.set_ylim([min_y_spectrum,max_y_spectrum]) 
    ax1.legend(loc='best',fontsize=14)
    ax1.set_xscale('log', base=10)
    ax1.set_xlabel('Energy (kev)',fontsize=22)
    ax1.set_ylabel(label,fontsize=22)

plt.tight_layout()
plt.show()
