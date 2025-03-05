import numpy as np
from  matplotlib import pyplot as plt
import f2py_interface as ib
import os



#-----------------------------#
#set env variables for tests
os.environ["REV_VERB" ] = "2"
os.environ["TEST_RUN" ] = "0"
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
Emax = 1000.0
ne = 1000
# ear = np.logspace(np.log10(Emin), np.log10(Emax), ne, dtype = np.float32
ear = np.zeros(ne, dtype = np.float32)
for i in range(ne):
    ear[i] = Emin * (Emax/Emin)**(i/ne)

# ear = np.genfromtxt('Benchmarks/xrb/Spec/Total.dat').T[0]
    
param = np.zeros(21, dtype = np.float32)

param[0]  = 6.0     #h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
param[1]  = 0.998   #a     !BH spin
param[2]  = 30.0    #inc   !Inclination angle in degrees
param[3]  = -1.0    #rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
param[4]  = 1e3     #rout  !Disk outer radius in Rg - will probably hardwire this
param[5]  = 0.0     #zcos  !Cosmological redshift
param[6]  = 3.0     #Gamma !Photon index
param[7]  = 3.0     #logxi !log10xi - ionisation parameter
param[8]  = 1.0     #Afe   !Iron abundance      
param[9]  = 15      #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[10] = 60.0    #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
param[11] = 0.0     #Nh
param[12] = -1.0     #1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
param[13] = 4.6e7    #M     !BH mass in solar masses
param[14] = 0.0     #flo   !Lowest frequency in band (Hz)
param[15] = 0.0     #fhi   !Highest frequency in band (Hz)
param[16] = -1.0     #ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
param[17] = 0.0     #DelA
param[18] = 0.0     #DelAB
param[19] = 0.0     #gamma
param[20] = 1       #telescope response         
          
# param[0]  = 29.7014      #h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
# param[1]  = 0.5          #a     !BH spin
# param[2]  = 44.2792      #inc   !Inclination angle in degrees
# param[3]  = -1           #rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
# param[4]  = 1000         #rout  !Disk outer radius in Rg - will probably hardwire this
# param[5]  = 0            #zcos  !Cosmological redshift
# param[6]  = 1.79279      #Gamma !Photon index
# param[7]  = 1.44735      #logxi !log10xi - ionisation parameter
# param[8]  = 1            #Afe   !Iron abundance      
# param[9]  = 19.9874      #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
# param[10] = 60           #kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
# param[11] = 0.0          #Nh
# param[12] = 1            #1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
# param[13] = 10           #M     !BH mass in solar masses
# param[14] = 0.12207      #flo   !Lowest frequency in band (Hz)
# param[15] = 0.24414      #fhi   !Highest frequency in band (Hz)
# param[16] = 4            #ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
# param[17] = -0.0166535   #DelA
# param[18] = -0.686551    #DelAB
# param[19] = 0.15228      #gamma
# param[20] = 1            #telescope response
              

plt.ion()
print('')
print('')
print('')
print('---------------------------------------------------------')
# name_input = './Benchmarks/xrb/ip_0,12_0,25.dat'
name_input = './Benchmarks/test_parametes.dat'
model_type = 'xrb'
# name_input = './Benchmarks/test_par_rtdist.dat'
# model_type = 'rtdist'


print (f'reading input parameters in {name_input} file ')
#parameters = np.genfromtxt(name_input, dtype = np.float32)
parameters = param
print('')            
print('*********************************************************')
print(f'running model for {model_type}  mode')

# parameters = np.array([6,0.9,57,-2,2e4,0.024917,2.5,1e5,1,17,50.,5e-2,1,3e6,0.02,0,0,0,0,0,0,-0.8,0.3,2.2e-4,1,1.])


match model_type:
    case 'xrb':
        photar_test = ib.reltransDCp(ear, parameters)
    case 'dbl':
        photar_test = ib.reltransDbl(ear, parameters)
    case 'agn':
        photar_test = ib.reltransDCp(ear, parameters)
    case 'rtdist':
        photar_test = ib.rtdist(ear, parameters)

print('')
print('')
print('')
print('RELTRANS finished!')
print('*********************************************************')
print('')
# #Print the two models 
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
font = 20

ax.plot(ear[:-1], photar_test, lw = 3, ls = '-' , label = 'testing reltrans')
print(photar_test)
    
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Energy [keV]', fontsize = font )
ax.set_ylabel(r'Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$', fontsize = font )
ax.tick_params(which='major', width=2, length=8, labelsize = font, pad=10)
ax.tick_params(which='minor', width=2, length=5, labelsize = font, pad=10)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
ax.yaxis.set_ticks_position('both')
ax.legend(fontsize = 10)
plt.savefig("Output/spec.png")
plt.close()
# print(photar_test)

            # E = (ear[1:] + ear[:-1]) * 0.5 
            # dE = (ear[1:] - ear[:-1])

            # if (_mode != 'Spec'):        
            #     ax.set_title('Plot source: '+ str(_source)+', mode: '+ str(_mode) \
            #                  + ', freq: ' + str(_frange), fontsize = font)
            #     for _model_type in model_type:
            #         name_archive = './Benchmarks/' + str(_source) + '/' \
            #             + str(_mode) + '/' + str(_model_type) + '_' + \
            #             str(_frange) + '.dat'
            #         print(f'Then comparing with the stored data: {name_archive}')
            #         data_archive = np.genfromtxt(name_archive).T

            #         name_model = './Output/' + str(_model_type) + '.dat'
            #         data_model = np.genfromtxt(name_model).T
            #         plot(ax, data_archive, data_model, sub_type = _model_type, identifier=str(_mode)+'_'+str(_frange))
            #         # plot(ax, data_archive, photar_test)
            # else:
            #     name_archive = './Benchmarks/'+str(_source)+'/'+str(_mode)+'/Total.dat'
            #     print(f'Then comparing with the stored data: {name_archive}')
            #     data_archive = np.genfromtxt(name_archive).T       

            #     name_model = './Output/Total.dat'
            #     data_model = np.genfromtxt(name_model).T
            #     ax.set_title('Plot source: '+ str(_source)+', mode: '+ str(_mode) \
            #                  + ', freq: ' + str(_frange), fontsize = font)
            #     plot(ax, data_archive, data_model, identifier=str(_mode)+'_'+str(_frange))

                
            # ax.set_xlim(0.1, 1e3)
            # ax.set_ylim(0.99, 1.01)
            # ax.set_ylim(1e-2, 1e3)

depad = np.loadtxt("Output/depad_re.dat")
plt.plot(depad[:,1])
plt.xlabel("nex")
plt.ylabel("photons/cm^2/s/keV?")
plt.yscale("log")
plt.savefig("Output/depad.png")
plt.close()

conv_re = np.loadtxt("Output/convolution_re.dat")
conv_im = np.loadtxt("Output/convolution_im.dat")
plt.plot(conv_re[:,1],label="Real")
plt.plot(conv_im[:,1],label="Imag",ls="--")
plt.xlabel("nec")
plt.ylabel("units")
plt.legend()
plt.savefig("Output/convolution.png")
plt.close()

conv = conv_re[:,1]+conv_im[:,1]*1j
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
plt.plot(np.angle(conv),label="Argument",lw=0.1)
plt.xlabel("nec")
plt.ylabel("Argument")
plt.legend()
plt.savefig("Output/convolution_arg.png")
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
plt.plot(np.abs(conv),label="Modulus")
plt.xlabel("nec")
plt.ylabel("Modulus")
plt.legend()
plt.savefig("Output/convolution_mod.png")
plt.close()