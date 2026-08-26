# Reltrans

[![Test](https://github.com/reltrans/reltrans/actions/workflows/test.yml/badge.svg)](https://github.com/reltrans/reltrans/actions/workflows/test.yml)
[![Docs](https://img.shields.io/badge/docs-readthedocs-blue.svg?style=flat)](https://reltransdocs.readthedocs.io/en/)

Reltrans is an X-ray reverberation mapping model that accounts for the effects
of spacetime curvature near the black hole. It is designed to be applicable to
both active galactic nuclei (AGN) and black hole X-ray binaries (BHXRB). It can
be used to model both time averaged spectra, and energy-dependent cross
spectra. Depending on the model flavour, the model can be used to characterize
accretion flows and to estimate black hole mass, spin, and source distance.

Reltrans supports Linux and MacOS, and is open source under the MIT license.

## Quickstart

```bash
git clone https://github.com/reltrans/reltrans
cd reltrans

# Install dependencies
pip install xspectrampoline kerrz-lib

# And compile:
make
```

> [!IMPORTANT]
> To run the model you require re-normalised [xillver
> tables](https://sites.srl.caltech.edu/~javier/xillver/). Our testing tables
> can be fetched using:
> ```
> make fetch-tables
> ```
> Consult the [documentation](https://reltransdocs.readthedocs.io/en/) for more information.
>
> The `RELTRANS_TABLES` environment variable must then be configured to point
> to the directory containing the tables. For example:
> ```
> export RELTRANS_TABLES="$(pwd)/cache/tables"
> ```

### XSPEC

To compile for use with XSPEC, in addition to the above also run:
```bash
make xspec
```

The instructions for loading the model are then printed to the terminal.

### Python

To use the Python bindings, use
```bash
make python
```
and install into your Python environment with `pip install .` Example usage:

```python
import pyreltrans
import numpy as np

# Load the library
rt = pyreltrans.Reltrans()

params = pyreltrans.DCP_Parameters()
energy = np.geomspace(0.2, 10.0, 200)

# Call the model
flux = rt.dcp(energy, params)
```

## Documentation

You can find the documentation at [this
link](https://reltransdocs.readthedocs.io/en/). The documentation provides
detailed installation and setup in instructions, and describes the various
flavours of the model that are available to the users.

## Bugs and feature requests

We welcome comments and contributions from the community. The best way to get
in touch is via the [issues](https://github.com/reltrans/reltrans/issues) page,
which you can use both to report bugs or odd model behaviour, and to request or
discuss features.

## Citing Reltrans

If you use Reltrans in your research, we ask that you please cite the following papers:
- [Ingram et al., 2019 MNRAS, 488, 324](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488..324I), ([DOI](https://doi.org/10.1093/mnras/stz1720))
- [Mastroserio et al., 2021 MNRAS, 507, 55M](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507...55M/) ([DOI](https://doi.org/10.1093/mnras/stab2056))

If you use the double lamppost model flavour, please cite:

- [Lucchini et al., 2023 ApJ, 951, 19](https://ui.adsabs.harvard.edu/abs/2023ApJ...951...19L) ([DOI](https://doi.org/10.3847/1538-4357/acd24f))

Additionally, if you use the RTDist flavour we ask that you also cite:

- [Ingram et al., 2022 MNRAS, 509, 619](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..619I) ([DOI](https://doi.org/10.1093/mnras/stab2950))

All of the above papers are included in the [CITATION.bib](./CITATION.bib) file.
