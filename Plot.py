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

if (np.abs(Input[19]) == 1):
    label = "Real Comp"
elif (np.abs(Input[19]) == 2):
    label = "Imaginary Comp"
elif (np.abs(Input[19]) == 3 or np.abs(Input[19]) == 5):
    label = "Modulus"   
elif (np.abs(Input[19]) == 4 or np.abs(Input[19]) == 6):
    label = "Lag(s)"    
else:
    label = "Unclear"   
   
Total = np.genfromtxt("Output/Total.dat")  
PivotingPL = np.genfromtxt("Output/PivotingPL.dat")
LightTravelTime = np.genfromtxt("Output/LightTravelTime.dat")
PivotingReflection = np.genfromtxt("Output/PivotingReflection.dat")

min_x = 0.5
max_x = 10.

max_y = -1e20
min_y= 1e20

#this looks for the maximum/minimum values of the total output of the code in the x-axis interval selected
for i in range(len(Total.T[0])):
    if (Total.T[1][i] < min_y and Total.T[0][i] > min_x and Total.T[0][i] < max_x):
        min_y = Total.T[1][i]
    if (Total.T[1][i] > max_y and Total.T[0][i] > min_x and Total.T[0][i] < max_x):
        max_y = Total.T[1][i]    

if (min_y < 0):
    min_y = 1.5*min_y
else:
    min_y = 0.5*min_y
max_y = 1.5*max_y

fig, (ax1) = plt.subplots(1,1,figsize=(7.5,4.5))

ax1.plot(PivotingPL.T[0],PivotingPL.T[1],linewidth=2.5,label='Pivoting PL')
ax1.plot(Total.T[0],Total.T[1],linewidth=2.5,label='Total')
ax1.plot(PivotingReflection.T[0],PivotingReflection.T[1],linewidth=2.5,label='Pivoting Ref')
ax1.plot(LightTravelTime.T[0],LightTravelTime.T[1],linewidth=2.5,label='Reverberation')

ax1.set_xscale('log', base=10)
ax1.set_xlim([min_x,max_x])
ax1.set_ylim([min_y,max_y])

#ax1.set_yscale('log', base=10)
ax1.set_xlabel('Energy (kev)',fontsize=22)
ax1.set_ylabel(label,fontsize=22)
ax1.legend(loc='best',fontsize=14)

plt.tight_layout()
plt.show()
