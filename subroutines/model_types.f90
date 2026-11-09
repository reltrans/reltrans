module model_types
    !> Model parameters and configuration shared by the model routines.
    use rtconstants, only: wp
    use kerrz, only: krz_KerrMetric
    implicit none
    private
    public :: t_config, t_model_arguments

    ! This type is instantiated per-call, so is thread safe.
    ! TODO: This should be moved so that all of the wrappers pass the following
    ! instead of a params array
    type :: t_model_arguments
        ! These are hardcoded with a fixed size of 2, since nlp <= 2. by not
        ! having them be pointers or dynamically allocated, it makes reasoning
        ! about the code a little easier.
        real(wp) :: h(2)
        real(wp) :: DelAB(2), g(2)
        ! Black hole spin (dimensionless):
        real(wp) :: a
        ! Observer inclination in degrees:
        real(wp) :: inc
        ! Inner and outer radii of the accretion disc (rg):
        real(wp) :: rin, rout
        ! Scale height of the accretion disc
        real(wp) :: honr
        ! Cosmological redshift:
        real(wp) :: zcos
        ! Photon index
        real(wp) :: Gamma
        real(wp) :: logxi, Afe, lognep, Cutoff_obs, Cutoff_s
        real(wp) :: Dkpc, Anorm, beta_p
        real(wp) :: Nh, boost, Mass, floHz, fhiHz, DelA
        integer :: nlp, ReIm, resp_matr, Cp

        ! eta: Fourier frequency dependent normalisation ratio C1(vc) / C2(vc)
        ! See Lucchini et al. 2023 for details.
        real(wp) :: eta

        ! eta_0: Time-averaged normalisation ratio C1 / C2 between the continuum
        ! of the two lampposts. It sets the continuum cutoff and disc
        ! ionisation.
        ! See Lucchini et al. 2023 for details.
        real(wp) :: eta_0

        ! Asymmetry parameter of angular emissivity function:
        real(wp) :: qboost

        ! Linear coefficient of angular emissivity function (b1) and quadratic
        ! coefficient of angular emissivity function (b2):
        real(wp) :: b1, b2

        ! The below are computed from the above
        ! The cosine angle
        real(wp) :: muobs

        ! The kerrz metric used during computations.
        type(krz_KerrMetric) :: metric
    end type t_model_arguments

    ! This type is shared between calls so is not thread safe.
    type :: t_config
        integer :: verbose = 0
        ! firstcall: is this the first time the model has been called?
        logical :: firstcall = .true., needtrans = .true., needconv = .true.
        ! me: Number of mue bins
        ! xe: Number of logr bins: bins 1:xe-1 are logarithmically spaced, bin
        ! xe is everything else
        ! ionvar: ionvar=1 means include ionization variations, ionvar=0 means
        ! don't.
        integer :: me, xe, nex, m, ionvar, refvar

        ! TODO: are these really constants, or are they hidden variables?
        ! Constants

        ! Relativistic grid size, i.e. pixel count splitting up the image plane
        ! in polar coordinates.
        integer :: nphi = 200, nro = 200
        ! Non-relativistic grid size:
        integer :: nron = 100, nphin = 100

        ! The number of time bins used for calculating the impulse response
        ! function:
        integer :: nt = 2**9

        ! Emin, Emax: min and max of the internal energy grid
        ! (different from output grid)
        real(wp) :: Emin = 1e-2_wp, Emax = 3e3_wp, dyn = 0.0_wp

        ! rnmax: max radius for which GR ray-tracing is used
        ! dlogf: a resolution parameter in base 10
        real(wp) :: rnmax = 300.0_wp, dlogf = 0.09_wp

        real(wp) :: DeltaGamma = 0.01_wp

        ! Use ring-like coronal model. Should future models get added, this
        ! could be promoted to an enumeration of some description.
        logical :: ring_like = .false.

        ! Toggle whether to calculate the impulse response or not. Since the
        ! impulse response is not directly used in calculations, and must be
        ! read out in another way, this defaults to false.
        logical :: calculate_impulse_response = .false.

        ! internal frequency grid
        ! Number of frequency bins
        integer :: nf
        real(wp) :: f, fac
        real(wp) :: fc, flo, fhi
        ! internal frequency grid, for when we do lag/frequency spectra
        integer :: fbinx
        real(wp), allocatable :: fix(:)

        ! relativistic parameters and limit on rin and h
        real(wp) :: rmin, rh
        real(wp), allocatable :: height(:)

        ! internal energy grid (nex) and output/xspec (ne) energy grid
        ! dloge: logarithmic resolution of the internal energy grid
        real(wp) :: E, dE, dloge

        ! Radial and angle profile
        integer :: mubin, rbin, ibin

        ! variable for non linear effects
        ! DC == 1 means calculate the DC component, else opposite.
        integer :: DC, ionvariation
        real(wp) :: dlogxi1, dlogxi2

    end type t_config

end module model_types
