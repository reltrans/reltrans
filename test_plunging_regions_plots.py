import numpy as np
from  matplotlib import pyplot as plt
import f2py_interface as ib
import os



#-----------------------------#
#set env variables for tests
os.environ["REV_VERB" ] = "1"
os.environ["TEST_RUN" ] = "0"
os.environ["MU_ZONES" ] = "1"
os.environ["ION_ZONES"] = "1"
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
Emax = 1000.0
ne = 1000
# ear = np.logspace(np.log10(Emin), np.log10(Emax), ne, dtype = np.float32
ear = np.zeros(ne, dtype = np.float32)
for i in range(ne):
    ear[i] = Emin * (Emax/Emin)**(i/ne)

# ear = np.genfromtxt('Benchmarks/xrb/Spec/Total.dat').T[0]
    
param = np.zeros(21, dtype = np.float32)


param[0]  = 3.0001    #h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
param[1]  = 0.0     #a     !BH spin
param[2]  = 80.0    #inc   !Inclination angle in degrees
param[3]  = 1.0    #rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
param[4]  = 2e4     #rout  !Disk outer radius in Rg - will probably hardwire this
param[5]  = 0.0     #zcos  !Cosmological redshift
param[6]  = 2.0     #Gamma !Photon index
param[7]  = 3.0     #logxi !log10xi - ionisation parameter
param[8]  = 1.0     #Afe   !Iron abundance     
param[9]  = 15      #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[10] = 60.0    #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[11] = 0.0     #Nh
param[12] = 1.0     #1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
param[13] = 4.6e7    #M     !BH mass in solar masses
param[14] = 0.0     #flo   !Lowest frequency in band (Hz)
param[15] = 0.0     #fhi   !Highest frequency in band (Hz)
param[16] = 1.0     #ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
param[17] = 0.0     #DelA
param[18] = 0.0     #DelAB
param[19] = 0.0     #gamma
param[20] = 1       #telescope response  
              

print()
print()
print()

parameters = param
print('*********************************************************')
print('Start runnign reltrans')            

h_list = [3.0001, 5.0, 10.0, 25.0, 100.]
spin_list = [0.0, 0.5, 0.99]

# for a in spin_list:
#     for h in h_list:
        # param[0] = h
        # param[1] = a
        # photar_test = ib.reltransDCp(ear, parameters)
        # namefile_GR = f'GR_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rinRh.dat'
        # namefile_NW = f'NW_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rinRh.dat'
        # os.system(f'mv GR_loop.dat {namefile_GR}')
        # os.system(f'mv NW_loop.dat {namefile_NW}')
        # print()
        # print(f'Parameters that change: spin = {param[1]} and h = {param[0]}')

plt.ion()
# #Print the two models 
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
font = 20

param[0] = 5.0 #h
param[1] = 0.0 #a
param[2]  = 80.0    #inc  
param[3]  = 1.0    #rin   
param[4]  = 6.0 #Rout
photar_test = ib.reltransDCp(ear, parameters)
namefile_GR = f'data_test_GR_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rin{int(param[3])}_rout{int(param[4])}.dat'
namefile_NW = f'data_test_NW_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rin{int(param[3])}_rout{int(param[4])}.dat'
os.system(f'mv data_test_GR_loop.dat {namefile_GR}')
os.system(f'mv data_test_NW_loop.dat {namefile_NW}')
print()
print(f'Parameters that change: spin = {param[1]} and h = {param[0]}')

ax.plot(ear[:-1], photar_test* ear[:-1]**2 , lw = 3, ls = '-' , label = f'h{int(param[0])}_a{param[1]:.1f} incl{int(param[2])} rin{int(param[3])} rout{int(param[4])}')

param[0] = 5.0 #h
param[1] = 0.0 #a
param[2]  = 80.0    #inc  
param[3]  = 6.0    #rin   
param[4]  = 2000 #Rout
photar_test2 = ib.reltransDCp(ear, parameters)
namefile_GR = f'data_test_GR_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rin{int(param[3])}_rout{int(param[4])}.dat'
namefile_NW = f'data_test_NW_quantities_h{int(param[0])}_a{param[1]:.1f}_gamma{param[6]:.1f}_incl{int(param[2])}_rin{int(param[3])}_rout{int(param[4])}.dat'
os.system(f'mv data_test_GR_loop.dat {namefile_GR}')
os.system(f'mv data_test_NW_loop.dat {namefile_NW}')
print()
print(f'Parameters that change: spin = {param[1]} and h = {param[0]}')


print()
print('RELTRANS finished!')
print('*********************************************************')
print()

# with open ('test_final_model.txt', 'w') as out:
#     E = (ear[1:] + ear[:-1]) * 0.5 
#     for i, phot in enumerate(photar_test):
#         out.writelines(str(E[i]) + ' ' + str(phot) + '\n')

ax.plot(ear[:-1], photar_test2* ear[:-1]**2, lw = 3, ls = '-' , label = f'h{int(param[0])}_a{param[1]:.1f} incl{int(param[2])} rin{int(param[3])} rout{int(param[4])}')

ax.set_ylim(0.05,20) 
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Energy [keV]', fontsize = font )
ax.tick_params(which='major', width=2, length=8, labelsize = font, pad=10)
ax.tick_params(which='minor', width=2, length=5, labelsize = font, pad=10)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
ax.yaxis.set_ticks_position('both')
ax.legend(fontsize = 10)


input('Press Enter')

