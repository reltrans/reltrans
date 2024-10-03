import numpy as np
from  matplotlib import pyplot as plt
import f2py_interface as ib
import os


def plot(ax, data_archive, test_model, sub_type = 'Spectrum', identifier = ''):
    # Emin = 0.1
    # Emax = 500.0
    # ne = 5001
    # ear = np.logspace(np.log10(Emin), np.log10(Emax), ne, dtype = np.float32)
    # data = np.genfromtxt(name_archive)            
    ax.plot(test_model[0], test_model[1], lw = 3, ls = '-' , label = 'testing reltrans (type: ' + str(sub_type) + ')')
    
    ax.plot(data_archive[0], data_archive[1], lw = 3, linestyle = '--', label = 'archive (old) reltrans (type: ' + str(sub_type) + ')')
    
    ax.set_xscale('log')
    ax.set_xlabel(r'Energy [keV]', fontsize = font )
    ax.tick_params(which='major', width=2, length=8, labelsize = font, pad=10)
    ax.tick_params(which='minor', width=2, length=5, labelsize = font, pad=10)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3)
    ax.yaxis.set_ticks_position('both')
    ax.legend(fontsize = 10)

    with open('Output/compare_' + str(identifier) + '_' + str(sub_type) + '.dat', 'w') as file:
        for i, _ in enumerate(test_model[0]):
            line = f'{data_archive[0][i]} {data_archive[1][i]} {test_model[0][i]} {test_model[1][i]}'
            file.write(line)
            file.write('\n')

#-----------------------------#
#set env variables for tests
os.environ["REV_VERB" ] = "2"
os.environ["TEST_RUN" ] = "1"
os.environ["MU_ZONES" ] = "1"
os.environ["ION_ZONES"] = "20"
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
Emax = 200.0
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


source     = ['xrb','dbl','agn']
frange     = [['0,12_0,25', '0,31_0,73', '0,80_2,10', '2,10_5,80', '5,80_16,0'], \
              ['0,10_0,40', '0,50_0,60', '1,10_1,40', '3,00_4,20', '3,00_4,20'], \
              ['xmm_01_06', 'nus_01_06', 'xmm_06_20', 'nus_06_20', 'xmm_20_90', 'nus_20_90']]
mode       = ['Lags', 'Mods', 'Imag', 'Real', 'Spec'] 
model_type = ['Total', 'IonVar', 'PivPL', 'PivRef', 'Reverb']

# source = ['agn']
# frange = [['xmm_01_06', 'nus_01_06', 'xmm_06_20', 'nus_06_20', 'xmm_20_90', 'nus_20_90']]
# mode       = ['Lags', 'Mods', 'Imag', 'Real', 'Spec'] 
# model_type = ['Total']

# source = ['dbl']
# frange = [['0,10_0,40', '0,50_0,60', '1,10_1,40', '3,00_4,20', '3,00_4,20']]
# mode       = ['Lags', 'Mods', 'Imag', 'Real', 'Spec'] 
# model_type = ['Total', 'IonVar', 'PivPL', 'PivRef', 'Reverb']

# source = ['xrb']
# frange = [['0,12_0,25', '0,31_0,73', '0,80_2,10', '2,10_5,80', '5,80_16,0']]
# mode       = ['Lags', 'Mods', 'Imag', 'Real', 'Spec'] 
# model_type = ['Total', 'IonVar', 'PivPL', 'PivRef', 'Reverb']


plt.ion()
for j, _source in enumerate(source):
    for _frange in frange[j]:
        print('')
        print('')
        print('')
        print('---------------------------------------------------------')
        name_input = './Benchmarks/'+str(_source)+'/ip_'+str(_frange)+'.dat'
        print (f'reading input parameters in {name_input} file ')
        parameters = np.genfromtxt(name_input, dtype = np.float32)

        for i, _mode in enumerate(mode):

            # #Print the two models 
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            font = 20

            print('')            
            print('*********************************************************')
            print(f'running model for {_source} in {_mode} mode')
            match _mode:
                case 'Spec':
                    if (_source == 'dbl'):
                        parameters[18] = 0.0 #flo !Lowest  freq in band (Hz)
                        parameters[19] = 0.0 #fhi !Highest freq in band (Hz)
                        parameters[20] = 1.0 #ReIm !1=Re,2=Im,3=modulus, 4=time
                    else:
                        parameters[14] = 0.0 #flo !Lowest  freq in band (Hz)
                        parameters[15] = 0.0 #fhi !Highest freq in band (Hz)
                        parameters[16] = 1.0 #ReIm !1=Re,2=Im,3=modulus, 4=time
                    ax.set_yscale('log')
                    ax.set_ylabel(r'keV * photons / $cm^2$ / s / keV ', fontsize = font )
                    
                case 'Real':
                    if (_source == 'dbl'):
                        parameters[20] = 1.0
                    else:
                        parameters[16] = 1.0
                    # ax.set_yscale('log')
                    ax.set_ylabel(r'Real Part ', fontsize = font )
                case 'Imag':
                    if (_source == 'dbl'):
                        parameters[20] = 2.0
                    else:
                        parameters[16] = 2.0
                    ax.set_ylabel(r'Imaginary Part ', fontsize = font )
                case 'Mods':
                    if (_source == 'dbl'):
                        parameters[20] = 3.0
                    else:
                        parameters[16] = 3.0
                    # ax.set_yscale('log')
                    ax.set_ylabel(r'Modulus ', fontsize = font )
                case 'Lags':
                    if (_source == 'dbl'):
                        parameters[20] = 4.0
                    else:
                        parameters[16] = 4.0
                    ax.set_ylabel(r'Time Lag ', fontsize = font )
                    
            print('')            
            if (_source == 'dbl'):
                print(f' for {_mode} mode the freqneucy range is [\
                {parameters[18]} - {parameters[19]}] and ReIm is {parameters[20]}')
            else:
                print(f' for {_mode} mode the freqneucy range is [\
                {parameters[14]} - {parameters[15]}] and ReIm is {parameters[16]}')
            # os.environ["TEST_RUN" ] = "0"
            # photar = ib.reltransDCp(ear, parameters)
            print('Running RELTRANS with the above configuration...')
            print('')
            match _source:
                case 'xrb':
                    photar_test = ib.reltransDCp(ear, parameters)
                case 'dbl':
                    photar_test = ib.reltransDbl(ear, parameters)
                case 'agn':
                    photar_test = ib.reltransDCp(ear, parameters)

                    #the model has written the outputs in the folder Output
                    #since we are using REV_VERB = 2 
            print('')
            print('')
            print('')
            print('RELTRANS finished!')
            print('*********************************************************')
            print('')

            E = (ear[1:] + ear[:-1]) * 0.5 
            dE = (ear[1:] - ear[:-1])

            if (_mode != 'Spec'):        
                ax.set_title('Plot source: '+ str(_source)+', mode: '+ str(_mode) \
                             + ', freq: ' + str(_frange), fontsize = font)
                for _model_type in model_type:
                    name_archive = './Benchmarks/' + str(_source) + '/' \
                        + str(_mode) + '/' + str(_model_type) + '_' + \
                        str(_frange) + '.dat'
                    print(f'Then comparing with the stored data: {name_archive}')
                    data_archive = np.genfromtxt(name_archive).T

                    name_model = './Output/' + str(_model_type) + '.dat' #this file has been written during the run of the model
                    data_model = np.genfromtxt(name_model).T
                    
                    plot(ax, data_archive, data_model, sub_type = _model_type, identifier=str(_source) +'_'+ str(_mode)+'_'+str(_frange))
                    # plot(ax, data_archive, photar_test)
            else:
                name_archive = './Benchmarks/'+str(_source)+'/'+str(_mode)+'/Total.dat'
                print(f'Then comparing with the stored data: {name_archive}')
                data_archive = np.genfromtxt(name_archive).T       

                name_model = './Output/Total.dat'
                data_model = np.genfromtxt(name_model).T
                ax.set_title('Plot source: '+ str(_source)+', mode: '+ str(_mode) \
                             + ', freq: ' + str(_frange), fontsize = font)
                plot(ax, data_archive, data_model, identifier=str(_mode)+'_'+str(_frange))

                
            # ax.set_xlim(0.1, 1e3)
            # ax.set_ylim(0.99, 1.01)
            # ax.set_ylim(1e-2, 1e3)

input('Press Enter')

