#!/usr/bin/env python3
"""Convert reltransDCp MCMC chains to include distance parameter."""

import numpy as np
from astropy.io import fits
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'test'))
from wrapper import Reltrans


def calculate_distance(params, reltrans, xe=20, verbose=0):
    # TODO: Implement distance calculation
    # Requires: genreltrans(dset=0) -> logxir, genreltrans(dset=1) -> logxieff
    # Distance = 10^(0.5*(logxir-logxieff)) / boost
    raise NotImplementedError(
        "Distance calculation requires internal Fortran routines. "
        "Use Fortran post_process executable instead."
    )


def convert_chain(input_chain, output_chain, param_columns, fixed_params):
    with fits.open(input_chain) as hdul:
        chain_hdu = hdul['CHAIN']
        chain_data = chain_hdu.data
        n_steps = len(chain_data)
        
        reltrans = Reltrans()
        distances = np.zeros(n_steps)
        
        for i in range(n_steps):
            params = {}
            for pname, col_idx in param_columns.items():
                if col_idx > 0:
                    params[pname] = chain_data[i][col_idx - 1]
                else:
                    params[pname] = fixed_params[pname]
            
            try:
                distances[i] = calculate_distance(params, reltrans)
            except NotImplementedError:
                print("ERROR: Use Fortran post_process executable")
                sys.exit(1)
        col_dkpc = fits.Column(name='Dkpc', format='D', unit='kpc', array=distances)
        new_columns = chain_hdu.columns + fits.ColDefs([col_dkpc])
        new_hdu = fits.BinTableHDU.from_columns(new_columns, name='CHAIN')
        primary_hdu = hdul[0].copy()
        new_hdul = fits.HDUList([primary_hdu, new_hdu])
        new_hdul.writeto(output_chain, overwrite=True)


def main():
    input_chain = "input_chain.fits"
    output_chain = "output_chain_with_distance.fits"
    
    param_columns = {
        'h': 7, 'a': 0, 'inc': 8, 'rin': 9, 'rout': 0, 'zcos': 0,
        'Gamma': 10, 'logxi': 11, 'Afe': 12, 'lognep': 13, 'kTe': 14,
        'Nh': 0, 'boost': 15, 'Mass': 16, 'Anorm': 17
    }
    
    fixed_params = {'a': 0.998, 'rout': 2e4, 'zcos': 0.0, 'Nh': 0.0}
    
    if not os.path.exists(input_chain):
        print(f"ERROR: Input chain file not found: {input_chain}")
        sys.exit(1)
    
    convert_chain(input_chain, output_chain, param_columns, fixed_params)


if __name__ == "__main__":
    main()