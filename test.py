import numpy as np
from  matplotlib import pyplot as plt
import f2py_interface as ib

Emin = 0.1
Emax = 500.0
ear = np.logspace(np.log10(Emin), np.log10(Emax), 1000, dtype = np.float32)

param = np.zeros(21, dtype = np.float32)

param[0]  = 6.0     #h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
param[1]  = 0.998   #a     !BH spin
param[2]  = 30.0    #inc   !Inclination angle in degrees
param[3]  = -1.0    #rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
param[4]  = 1e3     #rout  !Disk outer radius in Rg - will probably hardwire this
param[5]  = 0.0     #zcos  !Cosmological redshift
param[6]  = 2.0     #Gamma !Photon index
param[7]  = 3.0     #logxi !log10xi - ionisation parameter
param[8]  = 1.0     #Afe   !Iron abundance      
param[9]  = 15      #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[10] = 60.0    #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[11] = 0.0     #Nh
param[12] = 1.0     #1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
param[13] = 10.0    #M     !BH mass in solar masses
param[14] = 0.0     #flo   !Lowest frequency in band (Hz)
param[15] = 0.0     #fhi   !Highest frequency in band (Hz)
param[16] = 1.0     #ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
param[17] = 0.0     #DelA
param[18] = 0.0     #DelAB
param[19] = 0.0     #gamma
param[20] = 1       #telescope response

photar = ib.reltransDCp(ear, param)

print()
# for phot in photar:
#     print(phot)

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
font = 20 

E = (ear[1:] + ear[:-1]) * 0.5 

ax.plot(E, photar, lw = 5, color = 'k', label = 'reltrans spectrum')

font = 25
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'keV * photons / $cm^2$ / s / keV ', fontsize = font )
ax.set_xlabel(r'Energy [keV]', fontsize = font )
ax.tick_params(which='major', width=2, length=8, labelsize = font, pad=10)
ax.tick_params(which='minor', width=2, length=5, labelsize = font, pad=10)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
ax.yaxis.set_ticks_position('both')
# ax.set_xlim(0.1, 1e3)
# ax.set_ylim(0.99, 1.01)

plt.show()
