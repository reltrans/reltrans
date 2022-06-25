***************Installing reltrans***************

Once you have downloaded the reltrans source code, go into the top directory (reverberation_code) and unpack the fftw tar file"
> tar -xvf fftw-3.3.9.tar.gz
Then run the configure script:
> ./config.sh
You may need ot give the required permissions with a chmod command.
Then run the compile script. For macs this is:
> ./compile_xspec_mac.sh
For linux this is instead:
> ./compile_xspec_linux.sh
This should create a file called libreltrans.dylib.

From now on, to use reltrans, you just need to open xspec:
> xspec
go to the top reltrans directory (reverberation_code) and type:
XSPEC12>load libreltrans.dylib

Typing "mo ?" into the xspec prompt will then show you the models at your desposal. The models listed with a * are user defined models loaded for this session. You can see that the reltrans models are:
reltrans*            - reltrans with cut off power law continuum
reltransCp*        - reltrans with nthComp continuum
reltransD*          - reltrans with cut off power law continuum, density is a free parameter
reltransDCp*      - reltrans with nthComp continuum, density is a free parameter
reltransx*           - reltrans with cut off power law continuum, uses reflionx instead of xillver
rtdist*                - same a reltransDCp except distance is a parameter instead of ionisation parameter
rtdistx*              - like rtdist but uses reflionx instead of xillver
simrtdist*          - simulates lag-energy spectrum
xillver*               - restframe reflection spectrum for cut off power law continuum
xillverCp*           - restframe reflection spectrum fornthComp continuum
xillverD*            - cut off power law continuum, free density
xillverDCp*        - nthComp continuum, free density.

The most advanced models are reltransDCp and rtdist. The simplest is reltrans.


***************Using reltrans***************

Go to the top reltrans directory (reverberation_code) and type:
XSPEC12>load libreltrans.dylib

You can call the basic reltrans model with the command:
XSPEC12>mo reltrans

You will be prompted for parameters. Important ones to dicuss here are:
14:reltrans:fmin>
15:reltrans:fmax>
These are the minimum and maximum frequencies of the range you wish to model. If you are simply modelling the time-averaged spectrum (DC component), set both of these parameters to 0. If you are modelling a Fourier product, these must be non-zero.
16:reltrans:ReIm>
This parameter governs which Fourier product is output (only if fmin and fmax are both >0).
ReIm = 1: Real part of the cross spectrum.
ReIm = 2: Imaginary part of the cross spectrum.
ReIm = 3: Modulus of the cross spectrum.
ReIm = 4: Lag-energy spectrum.
ReIm = 5: Modulus of the folded cross spectrum.
ReIm = 6: Lag-energy spectrum calculated from the folded cross spectrum.
For each of these options, you will be asked for a response matrix and an energy range for the reference band.
These are required to calculate the phase of the reference band (i.e. the zero point of the lag spectrum).
To fit to data, either use ReIm=1 and ReIm=3 or ReIm=5 and/or ReIm=6. Do not fit to data using ReIm=3 or ReIm=4.
When plotting, it makes sense to use e.g. "iplot eemodel" for ReIm=1,2,3. For ReIm=4,5, instead use "iplot model". You will also need to un-log the y-axis.
You can alternatively use negative numbers here; e.g. ReIm=-4. In this case, the reference band phase is calculated assuming the instrument response is diagonal and no response matrix at all is needed.
20:reltrans:RESP>
This enables the user to simultaneously fit lags measured by two instruments. If you're only using one instrument, just fix this parameter to 1.
21:reltrans:norm>
XSPEC always adds a norm parameter. If you are using ReIm=1,2,3,5 then norm must be a free parameter. If you are using ReIm=4,6, then norm *must* be fixed to 1. For the rtdist model, norm must *always* be fixed to 1.


***************Using simrtdist***************

simrtdist is a xspec local model that simulates a lag-energy
spectrum. You call it like a noral local model but there are a few nuiances:

XSPEC12>mo simrtdist

You will then be prompted for model parameters. Some key parameters:

8:simrtdist:Dkpc>100000.0
Distance in kpc. This is important to get right, as it goes into the
calculation of the ionisation parameter.

20:simrtdist:coh2>0.8
              0         -1(         1)     -6.283     -6.283
              6.283      6.283

This is the squared coherence of the frequency band; \gamma^2 =
|cross-spectrum|^2 / (  Pr * Ps )/

21:simrtdist:phiA>
              0         -1(         1)     -6.283     -6.283   6.283      6.283

Always swet this to zero unless you *really* know what you're doing!

22:simrtdist:phiAB>
              0       0.01(      0.01)          0          0        0.5        0.5
23:simrtdist:g>
          7e-05      1e-07(     7e-07)      1e-12      1e-12      1e+10      1e+10

Continuum lag parameters (see Mastroserio et al 2021).

24:simrtdist:Anorm>2.2e-4
         130000         -1(      1300)      1e-10      1e-10      1e+10      1e+10

The normalisation of the continuum spectrum.

25:simrtdist:Texp>260000.0
           0.01         -1(    0.0001)      1e-10      1e-10      1e+10      1e+10

Exposure time in seconds.

26:simrtdist:pow>25.
              1         -1(      0.01)          1          1          2          2

Power (squared fractional rms / Hz) in the frequency
band. Technically, it is the variability power of the normalisation of
the continuum. It is vital to get this right. It has no bearing on the
model itself, only on the errors. See e.g. Fig 2 of Vaughan et al
(2011) for a typical AGN power spectrum in the correct units. Pow is
much larger for an AGN than for an XRB.

27:simrtdist:RESP>1
              1       0.01(      0.01)          0          0      1e+20      1e+24

Just set to unity unless you really know what you're doing.

28:simrtdist:norm>1

The xspec norm parameter *must* be fixed to unity. It is only here
because it's impossible to stop XSPEC from including it!

You will then be prompted with the following things (with example
answers written below):

 Enter name of the response file (with full path)
PN.rmf
 Enter name of the anciliary (arf) response file (with full path)
PN.arf
 Enter lower energy in reference band
0.5
 Enter upper energy in reference band
10.0
 Enter name of the background file (with full path)
PNbackground_spectrum.fits
 Enter BACKSCAL factor (enter 1 if you dont know what this is)
1.

This will simulate an XMM-Newton EPIC-pn observation with the
reference band being 0.5-10 keV. You must enter the requested
files. For example, typing "none" for the background will result in an
error.

You will then be prompted to enter a root name for the files that
contain the simulated products:

 Enter root name of simulation products
adameg
 -----------------------------------------------
 Outputs: adameg.dat, xadameg.dat
 command: flx2xsp xadameg.dat xadameg.pha xadameg.rsp
 -----------------------------------------------

You will see that files such as adameg.dat have been created. You're
not done yet though, because these files will be nonesense! You first
need to come up with a sensible broad energy grid to simulate on! You
can do this the usual way with:

XSPEC12>dummy 0.5 10.0 20

You will again be prompted to enter a root name for the products:

 Enter root name of simulation products
adameg
 -----------------------------------------------
 Outputs: adameg.dat, xadameg.dat
 command: flx2xsp xadameg.dat xadameg.pha xadameg.rsp
 -----------------------------------------------

In a different terminal, you can take a look at the *.dat file that is
produced. This is in qdp format. Data group 1 is the simulated
lag-energy spectrum and data group 2 is the model. You can use this to
decide if you want to change anything.

If you do want to change something, you can just do this in the xspec
session. For example, say you want to change a parameter. You do this
the usual way ("newpar x y"), and the model prompts you to enter the
root name and writes out new files.

Once you are happy with the output, you can create an xspec-readable
pha file with the command:

flx2xsp xadameg.dat xadameg.pha xadameg.rsp

(properly adjusted to account for your chosen root name). The exact
command you need will be printed to screen.

This will create a pha file: adameg.pha. You can fit this synthetic
data in the usual way:

XSPEC12>model rtdist

XSPEC12>data 1:1 adameg.pha

