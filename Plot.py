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

if(Input[21]==0):
    label = "E*F(E)"
elif (np.abs(Input[22]) == 1):
    label = "Real Comp"
elif (np.abs(Input[22]) == 2):
    label = "Imaginary Comp"
elif (np.abs(Input[22]) == 3 or np.abs(Input[22]) == 5):
    label = "Modulus"   
elif (np.abs(Input[22]) == 4 or np.abs(Input[22]) == 6 or Input[22] == 7):
    label = "Lag(s)"    
else:
    label = "Unclear"   
    
#check if second LP is worth including: if h1 is different from h2 AND if eta and eta_0 are not 0, set nlp to two, otherwise to 1

nlp = 1
if ( (Input[0] != Input [1]) and (Input[12] != 0) and (Input[13] !=0)):
    nlp = 2
   
Total = np.genfromtxt("Output/Total.dat") 
LightTravelTime = np.genfromtxt("Output/TravelTime.dat") 
TotalReflection = np.genfromtxt("Output/TotReflect.dat") 
Continuum = np.genfromtxt("Output/Continuum_spec.dat")
PivotingPL = np.zeros((nlp,len(Total.T[0]),2))
PivotingReflection =  np.zeros((nlp,len(Total.T[0]),2))

#print(PivotingPL)

if (Input[22] != 7):
    #loop over LPs: 
    for i in range(nlp):
        PivotingString = "Output/PivotingPL_" + str(i+1) + ".dat"
        ReflectionString = "Output/PivotingRE_" + str(i+1) + ".dat"
        PivotingPL_file = np.genfromtxt(PivotingString)
        PivotingReflection_file = np.genfromtxt(ReflectionString)
        for j in range(len(Total.T[0])):
            for m in range(2):
                PivotingPL[i][j][m] = PivotingPL_file[j][m]
                PivotingReflection[i][j][m] = PivotingReflection_file[j][m]
     
    min_x_timing = 0.5
    max_x_timing = 10.

    min_x_spectrum = 0.5
    max_x_spectrum = 80.

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
    for i in range(len(LightTravelTime.T[0])):
        if (LightTravelTime.T[1][i] < min_y_timing and LightTravelTime.T[0][i] > min_x_timing and LightTravelTime.T[0][i] < max_x_timing):
            min_y_timing = LightTravelTime .T[1][i]
        if (LightTravelTime.T[1][i] > max_y_timing and LightTravelTime.T[0][i] > min_x_timing and LightTravelTime.T[0][i] < max_x_timing):
            max_y_timing = LightTravelTime.T[1][i]
    for i in range(len(TotalReflection.T[0])):
        if (TotalReflection.T[1][i] < min_y_timing and TotalReflection.T[0][i] > min_x_timing and TotalReflection.T[0][i] < max_x_timing):
            min_y_timing = TotalReflection .T[1][i]
        if (TotalReflection.T[1][i] > max_y_timing and TotalReflection.T[0][i] > min_x_timing and TotalReflection.T[0][i] < max_x_timing):
            max_y_timing = TotalReflection.T[1][i]  
    
    for j in range(nlp):                
        for i in range(len(PivotingPL[j].T[0])):
            if (PivotingPL[j].T[1][i] < min_y_timing and PivotingPL[j].T[0][i] > min_x_timing and PivotingPL[j].T[0][i] < max_x_timing):
                min_y_timing = PivotingPL[j].T[1][i]
            if (PivotingPL[j].T[1][i] > max_y_timing and PivotingPL[j].T[0][i] > min_x_timing and PivotingPL[j].T[0][i] < max_x_timing):
                max_y_timing = PivotingPL[j].T[1][i]             
        for i in range(len(PivotingReflection[j].T[0])):
                if (PivotingReflection[j].T[1][i] < min_y_timing and PivotingReflection[j].T[0][i] > min_x_timing and PivotingReflection[j].T[0][i] < max_x_timing):
                    min_y_timing = PivotingReflection[j].T[1][i]
                if (PivotingReflection[j].T[1][i] > max_y_timing and PivotingReflection[j].T[0][i] > min_x_timing and PivotingReflection[j].T[0][i] < max_x_timing):
                    max_y_timing = PivotingReflection[j].T[1][i] 
     
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
    if (Input[16] > 1.e4):
        min_x_timing = 0.7e-5
        max_x_timing = 3.e-2
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
if(Input[21] != 0 and Input[22] != 7):
    ax1.plot(Total.T[0],Total.T[1],linewidth=3.0,label='Total',color=colors[0],zorder=2)
    ax1.plot(LightTravelTime.T[0],LightTravelTime.T[1],linewidth=2.0,label='Reverberation',color=colors[1],zorder=1)
    ax1.plot(TotalReflection.T[0],TotalReflection.T[1],linewidth=2.0,label='Total reflection',color=colors[2],zorder=1)
    for i in range(nlp):          
        ax1.plot(PivotingPL[i].T[0],PivotingPL[i].T[1],linewidth=2.0,label='Pivoting PL '+str(i+1),color=colors[3+i],zorder=1)
        #this is skipped by default to avoid the plot being a horrible mess 
        ax1.plot(PivotingReflection[i].T[0],PivotingReflection[i].T[1],linewidth=2.0,label='Pivoting Ref '+str(i+1),color=colors[5+i],zorder=1)
    ax1.plot(line_array,dashed_line,linestyle='dashed',linewidth=1.0,color='black')
    ax1.legend(loc='best',fontsize=14)
    ax1.set_xlim([min_x_timing,max_x_timing])
    ax1.set_ylim([min_y_timing,max_y_timing])
    ax1.set_xscale('log', base=10)
    if (Input[22] == 3 or Input[22] == 5):
        ax1.set_yscale('log', base=10)    
    ax1.set_xlabel('Energy (kev)',fontsize=22)
    ax1.set_ylabel(label,fontsize=22)
#lag-frequency plot
elif(Input[22] == 7):
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
