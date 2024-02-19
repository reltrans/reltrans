import numpy as np
from  matplotlib import pyplot as plt
import f2py_interface as ib
import os

#-----------------------------#
#set env variables for tests
os.environ["REV_VERB" ] = "2"
os.environ["TEST_RUN" ] = "1"
os.environ["MU_ZONES" ] = "1"
os.environ["ION_ZONES"] = "10"
os.environ["A_DENSITY"] = "0"
os.environ["EMIN_REF" ] = "0.3"
os.environ["EMAX_REF" ] = "10.0"
os.environ["EMIN_REF2"] = "0.3"
os.environ["EMAX_REF2"] = "10.0"
os.environ["SEED_SIM" ] = "-2851043"
os.environ["BACKSCL"  ] = "1.0"
os.environ["RMF_SET"  ] = "./Benchmarks/resp_matrix/nicer-rmf6s-teamonly-array50.rmf"
os.environ["ARF_SET"  ] = "./Benchmarks/resp_matrix/nicer-consim135p-teamonly-array50.arf"
#-----------------------------#

Emin = 0.1
Emax = 500.0
ne = 5001

ear = np.logspace(np.log10(Emin), np.log10(Emax), ne, dtype = np.float32)

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

# parameters = np.zeros(21, dtype = np.float32)

source = ['xrb']
frange = ['0,12_0,25', '0,31_0,73', '0,80_2,10', '2,10_5,80', '5,80_16,0']
mode   = ['Lags', 'Mods', 'Imag', 'Real', 'Spec'] 
mode_par   = [4, 3, 2, 1, 0]
model_type = ['Total']
for _source in source:
    for _frange in frange:
        name = './Benchmarks/'+str(_source)+'/ip_'+str(_frange)+'.dat'
        print (f'reading input parameters in {name} file ')
        parameters = np.genfromtxt(name, dtype = np.float32)

        for i, _mode in enumerate(mode):
            print(f'running model for {_source} in {_mode} mode')
            if (mode_par[i] == 0):
                parameters[14] = 0.0     #flo   !Lowest frequency in band (Hz)
                parameters[15] = 0.0     #fhi   !Highest frequency in band (Hz)
                parameters[16] = 1.0     #ReIm  !1=Re, 2=Im, 3=modulus, 4=time
            else:
                parameters[16] = mode_par[i]
            print(f' for {_mode} mode the freqneucy range is [\
            {parameters[14]}-{parameters[15]}] and ReIm is {parameters[16]}')
            photar = ib.reltransDCp(ear, parameters)

            if (mode_par[i] == 0):
                name_archive = './Benchmarks/'+str(_source)+'/'+str(_mode)+'/'\
                    +str(_model_type)+'_'+str(_frange)+'.dat'
                print(f'Then comparing with the stored data: {name_archive}')

                data = np.genfromtxt(name_archive)

                # Print the two models 
                fig, ax = plt.subplots(1, 1, figsize=(10, 8))
                font = 20

                E = (ear[1:] + ear[:-1]) * 0.5 
                dE = (ear[1:] - ear[:-1]) 

                ax.plot(E, photar/dE, lw = 5, color = 'red', label = 'reltrans spectrum')
                ax.plot(data[:,0], data[:,1], linewidth = 3, linestyle = '--', color = 'grey'   , label = 'phot (old)')

                # ax.set_xscale('log')
                # ax.set_yscale('log')
                ax.set_ylabel(r'keV * photons / $cm^2$ / s / keV ', fontsize = font )
                ax.set_xlabel(r'Energy [keV]', fontsize = font )
                ax.tick_params(which='major', width=2, length=8, labelsize = font, pad=10)
                ax.tick_params(which='minor', width=2, length=5, labelsize = font, pad=10)
                for axis in ['top','bottom','left','right']:
                    ax.spines[axis].set_linewidth(3)
                    ax.yaxis.set_ticks_position('both')
                    # ax.set_xlim(0.1, 1e3)
                    # ax.set_ylim(0.99, 1.01)
                    # ax.set_ylim(1e-2, 1e3)

# plt.show()
