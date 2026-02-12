|Docs| |License| |GitHub release|

Reltrans
========

Reltrans is an X-ray reverberation mapping model designed to be applicable to 
both AGN and black hole X-ray binaries. It can be used to model both time 
averaged spectra, and energy-dependent cross spectra. Depending on the model 
version, the model can be used to characterize accretion flows and to estimate 
black hole mass, spin, and source distance. 

Documentation
=============

You can find the documentation at this `link <https://reltransdocs.readthedocs.io/en/>`_. 
The documentation provides detailed installation and setup in instructions, and 
describes the various flavours of the model that are available to the users. 

Development Setup
=================

Reltrans requires HEASOFT/XSPEC libraries for native compilation.

If the shared library is not built, native tests will be skipped automatically.

To enable native tests:

1. Install and source HEASOFT.
2. Run ``make`` in the reltrans root directory.
3. Run ``pytest``.

If the shared library is missing, the test suite will provide
diagnostic information including the expected library path
and the detected ``HEADAS`` environment variable.

Reporting errors and requesting features
========================================

We welcome comments and contributions from the community. The best way to get in 
touch is via the `issues_` page, which you can use both to report bugs or odd
model behaviour, and to request or discuss features. 

Citing Reltrans
===============

If you use Reltrans in your research, we ask that you please cite the following papers:

.. raw:: html

    <script type="text/javascript">
        function copyIngramBib() {
            var bibtex = `@ARTICLE{2019MNRAS.488..324I,
                author = {{Ingram}, Adam and {Mastroserio}, Guglielmo and {Dauser}, Thomas and {Hovenkamp}, Pieter and {van der Klis}, Michiel and {Garc{\'\i}a}, Javier A.},
                 title = "{A public relativistic transfer function model for X-ray reverberation mapping of accreting black holes}",
               journal = {\mnras},
              keywords = {black hole physics, methods: data analysis, galaxies: active, X-rays: binaries, Astrophysics - High Energy Astrophysical Phenomena},
                  year = 2019,
                 month = sep,
                volume = {488},
                number = {1},
                 pages = {324-347},
                   doi = {10.1093/mnras/stz1720},
         archivePrefix = {arXiv},
                eprint = {1906.08310},
          primaryClass = {astro-ph.HE},
                adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.488..324I},
               adsnote = {Provided by the SAO/NASA Astrophysics Data System}
               }`;
            const el = document.createElement('textarea');
            el.value = bibtex;
            document.body.appendChild(el);
            el.select();
            document.execCommand('copy');
            document.body.removeChild(el);
       
        }
        
        function copyGulloBib() {
            var bibtex = `@ARTICLE{2021MNRAS.507...55M,
                author = {{Mastroserio}, Guglielmo and {Ingram}, Adam and {Wang}, Jingyi and {Garc{\'\i}a}, Javier A. and {van der Klis}, Michiel and {Cavecchi}, Yuri and {Connors}, Riley and {Dauser}, Thomas and {Harrison}, Fiona and {Kara}, Erin and {K{\"o}nig}, Ole and {Lucchini}, Matteo},
                 title = "{Modelling correlated variability in accreting black holes: the effect of high density and variable ionization on reverberation lags}",
               journal = {\mnras},
              keywords = {accretion, accretion discs, black hole physics, relativistic processes, X-rays: binaries, X-rays: galaxies, Astrophysics - High Energy Astrophysical Phenomena},
                  year = 2021,
                 month = oct,
                volume = {507},
                number = {1},
                 pages = {55-73},
                   doi = {10.1093/mnras/stab2056},
         archivePrefix = {arXiv},
                eprint = {2107.06893},
          primaryClass = {astro-ph.HE},
                adsurl = {https://ui.adsabs.harvard.edu/abs/2021MNRAS.507...55M},
               adsnote = {Provided by the SAO/NASA Astrophysics Data System}
               }`;
            const el = document.createElement('textarea');
            el.value = bibtex;
            document.body.appendChild(el);
            el.select();
            document.execCommand('copy');
            document.body.removeChild(el);       
        }        
    </script>

    <ul>
        <li>Ingram et al., 2019 MNRAS, 488, 324
            [<a href="https://doi.org/10.1093/mnras/stz1720">DOI</a>]
            [<a href="https://ui.adsabs.harvard.edu/abs/2019MNRAS.488..324I">ADS</a>]
            [<a onclick="copyIngramBib()">Copy BibTeX to clipboard</a>]</li>
            
        <li>Mastroserio et al., 2021 MNRAS, 507, 55M
            [<a href="https://doi.org/10.1093/mnras/stab2056">DOI</a>]
            [<a href="https://ui.adsabs.harvard.edu/abs/2021MNRAS.507...55M/">ADS</a>]
            [<a onclick="copyGulloBib()">Copy BibTeX to clipboard</a>]</li>
    </ul>
   
Additionally, if you use the RTDist flavour we ask that you also cite: 

.. raw:: html

    <script type="text/javascript">
        function copyRTdistBib() {
            var bibtex = `@ARTICLE{2022MNRAS.509..619I,
                author = {{Ingram}, Adam and {Mastroserio}, Guglielmo and {van der Klis}, Michiel and {Nathan}, Edward and {Connors}, Riley and {Dauser}, Thomas and {Garc{\'\i}a}, Javier A. and {Kara}, Erin and {K{\"o}nig}, Ole and {Lucchini}, Matteo and {Wang}, Jingyi},
                 title = "{On measuring the Hubble constant with X-ray reverberation mapping of active galactic nuclei}",
               journal = {\mnras},
              keywords = {black hole physics, methods: data analysis, galaxies: active, cosmological parameters, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Cosmology and Nongalactic Astrophysics},
                  year = 2022,
                 month = jan,
                volume = {509},
                number = {1},
                 pages = {619-633},
                   doi = {10.1093/mnras/stab2950},
         archivePrefix = {arXiv},
                eprint = {2110.15651},
          primaryClass = {astro-ph.HE},
                adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..619I},
               adsnote = {Provided by the SAO/NASA Astrophysics Data System}
               }`;
            const el = document.createElement('textarea');
            el.value = bibtex;
            document.body.appendChild(el);
            el.select();
            document.execCommand('copy');
            document.body.removeChild(el);       
        }
    </script>

    <ul>
        <li>Ingram et al., 2022 MNRAS, 509, 619
            [<a href="https://doi.org/10.1093/mnras/stab2950">DOI</a>]
            [<a href="https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..619I">ADS</a>]
            [<a onclick="copyRTdistBib()">Copy BibTeX to clipboard</a>]</li>
    </ul>

Finally, if you use the double lamp post flavour, we ask that you also cite:

.. raw:: html

    <script type="text/javascript">
        function copydoubleLPBib() {
            var bibtex = `@ARTICLE{2023ApJ...951...19L,
                author = {{Lucchini}, Matteo and {Mastroserio}, Guglielmo and {Wang}, Jingyi and {Kara}, Erin and {Ingram}, Adam and {Garcia}, Javier and {Dauser}, Thomas and {van der Klis}, Michiel and {K{\"o}nig}, Ole and {Lewin}, Collin and {Nathan}, Edward and {Panagiotou}, Christos},
                 title = "{Investigating the Impact of Vertically Extended Coronae on X-Ray Reverberation Mapping}",
               journal = {\apj},
              keywords = {Accretion, Black hole physics, Reverberation mapping, 14, 159, 2019, Astrophysics - High Energy Astrophysical Phenomena},
                  year = 2023,
                 month = jul,
                volume = {951},
                number = {1},
                   eid = {19},
                 pages = {19},
                   doi = {10.3847/1538-4357/acd24f},
         archivePrefix = {arXiv},
                eprint = {2305.05039},
          primaryClass = {astro-ph.HE},
                adsurl = {https://ui.adsabs.harvard.edu/abs/2023ApJ...951...19L},
               adsnote = {Provided by the SAO/NASA Astrophysics Data System}
               }`;
            const el = document.createElement('textarea');
            el.value = bibtex;
            document.body.appendChild(el);
            el.select();
            document.execCommand('copy');
            document.body.removeChild(el);       
        }
    </script>
    
    <ul>
        <li>Lucchini et al., 2023 ApJ, 951, 19
            [<a href="https://doi.org/10.3847/1538-4357/acd24f">DOI</a>]
            [<a href="https://ui.adsabs.harvard.edu/abs/2023ApJ...951...19L">ADS</a>]
            [<a onclick="copydoubleLPBib()">Copy BibTeX to clipboard</a>]</li>
    </ul>


.. |Docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://reltransdocs.readthedocs.io/en/
.. |License| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT   
.. |GitHub release| image:: https://img.shields.io/github/v/release/reltrans/reltrans
   :target: https://github.com/reltrans/reltrans/releases/latest   
   
.. _issues: https://github.com/reltrans/reltrans/issues
