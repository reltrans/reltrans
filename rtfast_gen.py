from joblib import dump, load
import f2py_interface as ib
import numpy as np
import os
import scipy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed

def lhc_filter(lhc,disk=False,spin=False,inner_r=False):
    """
    Removes unphysical parameter sets from the Latin Hypercube. This prevents
    overly bright sources from being generated as well as reducing time spent
    on generating model data for objects that we won't see

    Parameters
    ----------
    lhc : np.ndarray
        latin hypercube containing parameter sets for rtdist.

    Returns
    -------
    new_lhc : np.ndarray
        latin hypercube with unphysical parameter sets removed..

    """
    #bad sets indexes all parameter sets that don't fit the filter criteria
    #removes parameter sets that have a photon index higher than 3 AND a iron
    #solar abundance above 6 AND a electron density in the disk of higher than
    #10^19.
    if disk == True:
        bad_disks = np.nonzero((lhc[:,6]>2.75)&(10**lhc[:,7]>4)&(lhc[:,9]>17))
    #coronas within the blackhole horizon
    if spin == True:
        heights = 1+ np.sqrt(1-lhc[:,1]**2)
        bad_heights = np.nonzero(10**lhc[:,0]<1.5*heights)
    if inner_r == True:
        a = lhc[:,1]
        z1 = ( 1.0 - a**2.0 )**(1.0/3.0)
        z1 = z1 * ( (1.0+a)**(1.0/3.0)+(1.0-a)**(1.0/3.0))+1.0
        z2 = np.sqrt( 3.0 * a**2.0 + z1**2.0 )
        dISCO = 3.0 + z2 - np.sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
        bad_inner_r = np.nonzero(10**lhc[:,3]*dISCO>300)
    #collates all bad sets together
    if (spin==True)&(disk==True)&(inner_r==True):
        bad_sets = np.unique(np.concatenate((bad_disks,bad_heights,bad_inner_r),axis=None))
        new_lhc = np.delete(lhc,bad_sets,0)
    elif (spin==True)&(disk==True):
        bad_sets = np.unique(np.concatenate((bad_disks,bad_heights),axis=None))
        #removes all unphysical sets from the parameter sets
        new_lhc = np.delete(lhc,bad_sets,0)
    elif (spin==True)&(inner_r==True):
        bad_sets = np.unique(np.concatenate((bad_heights,bad_inner_r),axis=None))
        new_lhc = np.delete(lhc,bad_sets,0)
    elif (disk==True)&(inner_r==True):
        bad_sets = np.unique(np.concatenate((bad_disks,bad_inner_r),axis=None))
        new_lhc = np.delete(lhc,bad_sets,0)
    elif spin==True:
        new_lhc = np.delete(lhc,bad_heights,0)
    elif disk == True:
        new_lhc = np.delete(lhc,bad_disks,0)
    elif inner_r == True:
        new_lhc = np.delete(lhc,bad_inner_r,0)
    else:
        new_lhc = lhc
    return new_lhc

def lhc_generation(size,range_all):
    if size < 1e5:
        gen_size = int(1e5)
    else:
        gen_size = size
    lhc = lhc_cycle(gen_size, range_all)
    percent = 0
    while lhc.shape[0] < size:
        if percent <= ((lhc.shape[0]/size)*100 - 10):
            print(f"Currently {lhc.shape[0]}/{size}.")
            percent = np.round((lhc.shape[0]/size)*100,0)/10
            percent = np.floor(percent)*10
        lhc_temp = lhc_cycle(gen_size, range_all)
        lhc = np.concatenate((lhc,lhc_temp),axis=0)
    np.random.shuffle(lhc)
    lhc = lhc[:size]
    return lhc

def lhc_cycle(size,range_all):
    sampler = scipy.stats.qmc.LatinHypercube(d=len(range_all))
    sample = sampler.random(n=size)
    theta_lhc = scipy.stats.qmc.scale(sample, range_all[:,0], range_all[:,1])
    final_lhc = lhc_filter(theta_lhc,spin=True,disk=True,inner_r=True)
    return final_lhc

def nn_pars_to_rtdist(nn_pars, pars_list, negatives, logged):
    """
    Converts sampled paramters for neural network training into correct format
    for use in generating data and adds non-sampled parameters needed by the
    model. This is the generic form in which the list of sampled parameters 
    and manipulations are passed to the function

    Parameters
    ----------
    nn_pars : np.ndarray
        array of sampled parameters.
    ReIm : int
        value of ReIm from rtdist to determine the type of output rtdist 
        produces.
    pars_list : list
        list of indexes of parameters.
    negatives : list
        list of indexes of parameters that need to be turned positive.
    logged : list
        list of indexes of parameters that need to be transformed to power of
	10.

    Returns
    -------
    converted_pars : np.ndarray
        array of parameters to be used for parameter generation.

    """
    #set up base parameters which can be used to fix parameter sets to
    #reasonable values
    pars_base = [6,0.9,57,-1,2e4,0,2.45,3,1,17,50.,0,1,3e6,0,0,0,
                 0,0,0,1,1.]
    #set up base parameters to transform according to sampled parameters
    converted_pars = []
    for i in range(nn_pars.shape[0]):
        converted_pars.append(pars_base)
    converted_pars = np.asarray(converted_pars)
    #iterate through 
    for i, parameter in enumerate(pars_list):
        if parameter in logged:
            converted_pars[:,parameter] = 10**nn_pars[:,i]
        else:
            converted_pars[:,parameter] = nn_pars[:,i]
        if parameter in negatives:
            converted_pars[:,parameter] = -converted_pars[:,parameter]
    return converted_pars

def main():
    pars_list = [0,1,2,3,4,6,7,8,9,10]
    negatives = [3]
    logged = [0,2,3,4,10]
    
    height_range = [np.log10(3),np.log10(7e2)]
    spin_range = [0,0.998]
    inclination_range = [np.log10(5),np.log10(85)]
    r_inner_range = [np.log10(1),np.log10(200)]
    r_outer_range = [np.log10(400),np.log10(1e5)]
    Gamma_range = [1.4,2.7]
    logxi_range = [0,4.7]
    Afe_range = [0.5,10]
    logNe_range = [15,20]
    kte_range = [np.log10(30),np.log10(500)]
    
    
    range_all = [height_range,spin_range,inclination_range,r_inner_range,
             r_outer_range,Gamma_range,logxi_range,Afe_range,
             logNe_range,kte_range]
    range_all = np.asarray(range_all)
    
    n=int(5e6)
    lhc = lhc_generation(n,range_all)
    lhc = nn_pars_to_rtdist(lhc,pars_list,negatives,logged)
    lhc = lhc.astype(np.float32)
    
    #set env variables for tests
    os.environ["REV_VERB" ] = "0"
    os.environ["MU_ZONES" ] = "1"
    os.environ["ION_ZONES"] = "1"
    os.environ["A_DENSITY"] = "0"
    os.environ["BACKSCL"  ] = "1.0"
    os.environ["TEST_RUN" ] = "1"
    nn_dir = "/data/benr/rtfast-training/"
    
    Emin = 0.1
    Emax = 100.0
    ne = 1000
    egrid = np.zeros(ne, dtype = np.float32)
    for i in range(ne):
        egrid[i] = Emin * (Emax/Emin)**(i/ne)
    
    cpu_num = os.cpu_count()
    lhc_split = np.split(lhc,int(n/20000))
    start = 38
    if start != 0:
        spec_scaler = load(nn_dir+"scalers/scaler.bin")
        pca = load(nn_dir+"scalers/pca.bin")
    for i,lhc in enumerate(lhc_split[start:],start=start):
        with Parallel(n_jobs=cpu_num-2,verbose=1,backend="multiprocessing") as parallel:
            spectra =  parallel(delayed(ib.reltransDCp)(egrid,pars.astype(np.float32))
                                            for pars in lhc)
            spectra = np.asarray(spectra)
            
            mask = np.any(spectra <= 0,axis=1)
            spectra[spectra<=0] = 1e-11
            spectra = spectra[~mask]
            mask = np.any(np.isnan(spectra),axis=1)
            spectra = spectra[~mask]
            
            if i == 0:
                spec_scaler = StandardScaler()
                data = spec_scaler.fit_transform(np.log10(spectra))
                comps = 200
                pca = PCA(n_components=comps)
                pca_comps = pca.fit_transform(data)
                dump(spec_scaler,nn_dir+"scalers/scaler.bin")
                dump(pca,nn_dir+"scalers/pca.bin")
            else:
                data = spec_scaler.transform(np.log10(spectra))
                pca_comps = pca.transform(data)
            np.savetxt(nn_dir+f"data/pca_comps_{i}.txt",pca_comps)
            np.savetxt(nn_dir+f"data/pars_{i}.txt",lhc)


if __name__ == "__main__":
    main()
