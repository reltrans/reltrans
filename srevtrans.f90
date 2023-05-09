! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for \ell * g^{2+Gamma}
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) Ionisation profile
! 5) The dependence on incident angle (mimicked by tweaking the ionization)

!CURRENT BRANCH
! This branch is adding the multiple flavours to the reltrans model
  
  !  Cp : chooses xillver model
  !      -1 reltrans      1e15 density and powerlaw illumination  
  !       1 reltransD     high density and powerlaw illumination
  !      -2 reltransCp    1e15 density and nthcomp  illumination
  !       2 reltransDCp   high density and nthcomp  illumination
  
  
include 'subroutines/header.h'


!-----------------------------------------------------------------------
subroutine tdreltrans(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(20), photar(ne), par(32)
! Settings
  Cp   = 1   !|Cp|=1 means cut-off pl, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !logxi
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = 300.0            !Ecut_obs
  par(12) = param(11)        !Nh
  par(13) = param(12)        !boost
  par(14) = 1.0              !qboost
  par(15) = param(13)        !Mass
  par(16) = 0.0              !honr
  par(17) = 0.0              !b1
  par(18) = 0.0              !b2
  par(19) = param(14)        !floHz
  par(20) = param(15)        !fhiHz
  par(21) = param(16)        !ReIm
  par(22) = param(17)        !DelA
  par(23) = param(18)        !DelAB
  par(24) = param(19)        !g
  par(25) = 1.0              !Anorm
  par(26) = param(20)        !telescope response
  par(27) = 0.0              !First lc lowEn band
  par(28) = 0.0              !First lc highEn band
  par(29) = 0.0              !Second lc lowEn band
  par(30) = 0.0              !Second lc highEn band
  par(31) = 0.0              !AB slope 
  par(32) = 0.0              !g slope
! Call general code
  call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  return
end subroutine tdreltrans
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine tdreltransCp(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !logxi
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = param(13)        !boost
  par(14) = 1.0              !qboost
  par(15) = param(14)        !Mass
  par(16) = 0.0              !honr
  par(17) = 0.0              !b1
  par(18) = 0.0              !b2
  par(19) = param(15)        !floHz
  par(20) = param(16)        !fhiHz
  par(21) = param(17)        !ReIm
  par(22) = param(18)        !DelA
  par(23) = param(19)        !DelAB
  par(24) = param(20)        !g
  par(25) = 1.0              !Anorm
  par(26) = param(21)        !telescope response
  par(27) = 0.0              !First lc lowEn band
  par(28) = 0.0              !First lc highEn band
  par(29) = 0.0              !Second lc lowEn band
  par(30) = 0.0              !Second lc highEn band
  par(31) = 0.0              !AB slope 
  par(32) = 0.0              !g slope
! Call general code
  call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  return
end subroutine tdreltransCp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransx(ear,ne,param,ifl,photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !logxi
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = param(13)        !boost
  par(14) = 1.0              !qboost
  par(15) = param(14)        !Mass
  par(16) = 0.0              !honr
  par(17) = 0.0              !b1
  par(18) = 0.0              !b2
  par(19) = param(15)        !floHz
  par(20) = param(16)        !fhiHz
  par(21) = param(17)        !ReIm
  par(22) = param(18)        !DelA
  par(23) = param(19)        !DelAB
  par(24) = param(20)        !g
  par(25) = 1.0              !Anorm
  par(26) = param(21)        !telescope response
  par(27) = 0.0              !First lc lowEn band
  par(28) = 0.0              !First lc highEn band
  par(29) = 0.0              !Second lc lowEn band
  par(30) = 0.0              !Second lc highEn band
  par(31) = 0.0              !AB slope 
  par(32) = 0.0              !g slope
! Call general code
  call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  return
end subroutine tdreltransx
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine tdrtdist(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !Dkpc
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = 1.0              !boost
  par(14) = param(13)        !qboost
  par(15) = param(14)        !Mass
  par(16) = param(15)        !honr
  par(17) = param(16)        !b1
  par(18) = param(17)        !b2
  par(19) = param(18)        !floHz
  par(20) = param(19)        !fhiHz
  par(21) = param(20)        !ReIm
  par(22) = param(21)        !DelA
  par(23) = param(22)        !DelAB
  par(24) = param(23)        !g
  par(25) = param(24)        !Anorm
  par(26) = param(25)        !telescope response
  par(27) = 0.0              !First lc lowEn band
  par(28) = 0.0              !First lc highEn band
  par(29) = 0.0              !Second lc lowEn band
  par(30) = 0.0              !Second lc highEn band
  par(31) = 0.0              !AB slope 
  par(32) = 0.0              !g slope
! Check that we're not looking at the side of the disc
  honr = par(16)
  pi   = acos(-1.d0)
  cosi = cos( par(3) * pi / 180.d0 )
  cos0 = honr / sqrt( honr**2 + 1.d0  )  
! Call general code
  if( cos0 .ge. cosi )then
     photar = 0.0   !XSPEC *hates* this. Just do it with limits, and flag here.
     write(*,*)"Warning! Disc thickness is too high for this inclinaiton!"
     write(*,*)"Model output set to zero -- XSPEC *hates* this and may get lost"
     write(*,*)"leading to crash and seg fault. Better to set hard max on inc, incmax"
     write(*,*)"and set honr_max to cos(incmax)/sqrt(1-cos^2(incmax))."
  else
     call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  end if
  
  return
end subroutine tdrtdist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine tdrtdistX(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !Dkpc
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = 1.0              !boost
  par(14) = param(13)        !qboost
  par(15) = param(14)        !Mass
  par(16) = param(15)        !honr
  par(17) = param(16)        !b1
  par(18) = param(17)        !b2
  par(19) = param(18)        !floHz
  par(20) = param(19)        !fhiHz
  par(21) = param(20)        !ReIm
  par(22) = param(21)        !DelA
  par(23) = param(22)        !DelAB
  par(24) = param(23)        !g
  par(25) = param(24)        !Anorm
  par(26) = param(25)        !telescope response
! Check that we're not looking at the side of the disc
  honr = par(16)
  pi   = acos(-1.d0)
  cosi = cos( par(3) * pi / 180.d0 )
  cos0 = honr / sqrt( honr**2 + 1.d0  )  
! Call general code
  if( cos0 .ge. cosi )then
     photar = 0.0   !XSPEC *hates* this. Just do it with limits, and flag here.
     write(*,*)"Warning! Disc thickness is too high for this inclinaiton!"
     write(*,*)"Model output set to zero -- XSPEC *hates* this and may get lost"
     write(*,*)"leading to crash and seg fault. Better to set hard max on inc, incmax"
     write(*,*)"and set honr_max to cos(incmax)/sqrt(1-cos^2(incmax))."
  else
     call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  end if
  
  return
end subroutine tdrtdistX
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine simrtdist(ear, ne, param, ifl, photar)
  use telematrix
  use env_variables
  implicit none
  integer :: ne, ifl, Cp, dset, i
  real    :: ear(0:ne), param(27), photar(ne), par(32)
  real    :: gammac2, Texp, E, dE, getcountrate
  real    :: rephotar(ne), imphotar(ne)
  real, parameter :: Emin = 1e-1, Emax = 300.0
  integer, parameter :: nex=2**12
  real :: earx(0:nex),photarx(nex),pow
  real :: Pr,rephotarx(nex),imphotarx(nex),mur,mus
  real :: dlag(ne),G2,ReG,ImG,Psnoise,Prnoise,br,bs(ne)
  real :: flo,fhi,fc,lag(ne),gasdev,lagsim(ne)
  real, parameter :: pi = acos(-1.0)
  integer  unit,xunit,status,j
  real E1,E2,frac
  character (len=200) command,flxlagfile,phalagfile,rsplagfile,lagfile,root
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !Dkpc
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = 1.0              !boost
  par(14) = param(13)        !qboost
  par(15) = param(14)        !Mass
  par(16) = param(15)        !honr
  par(17) = param(16)        !b1
  par(18) = param(17)        !b2
  par(19) = param(18)        !floHz
  par(20) = param(19)        !fhiHz
  gammac2 = param(20)        !squared coherence
  par(22) = param(21)        !DelA
  par(23) = param(22)        !DelAB
  par(24) = param(23)        !g
  par(25) = param(24)        !Anorm
  Texp    = param(25)        !Texp (s)
  pow     = param(26)        !power in [rms/mean]^2/Hz units (alpha(nu))
  par(26) = param(27)        !telescope response
  par(27) = 0.0              !First lc lowEn band
  par(28) = 0.0              !First lc highEn band
  par(29) = 0.0              !Second lc lowEn band
  par(30) = 0.0              !Second lc highEn band
  par(31) = 0.0              !AB slope 
  par(32) = 0.0              !g slope
  
  flo = param(18)
  fhi = param(19)
  fc  = 0.5 * ( fhi + flo )
  
! Get `folded' lags
  par(21) = 6.0   !ReIm
  call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  do i = 1,ne
     lag(i) = photar(i) / ( ear(i) - ear(i-1) )
  end do
  
! Set internal energy grid
  do i = 0, nex
     earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
  end do
  
! Calculate real and imaginary parts on the fine energy grid
  par(21) = 1.0   !ReIm
  call genreltrans(Cp, dset, earx, nex, par, ifl, rephotarx)
  par(21) = 2.0   !ReIm
  call genreltrans(Cp, dset, earx, nex, par, ifl, imphotarx)
! Get DC component
  par(19) = 0.0   !floHz
  par(20) = 0.0   !fhiHz
  par(21) = 1.0   !ReIm
  call genreltrans(Cp, dset, earx, nex, par, ifl, photarx)

! Read in background array
  if( needbkg )then
     allocate(bkgcounts(1:numchn))
     allocate(bkgrate(1:numchn))
     call readinbkg
     needbkg = .false.
  end if

! Calculate background in reference band
  br = 0.0
  do i = ilo,ihi
     br = br + bkgrate(i)
  end do

! Calculate background in subject
  do j = 1,ne
     bs(j) = 0.0
     do i = 1,numchn
        if( ECHN(i) .gt. ear(j-1) .and. ECHN(i-1) .le. ear(j) )then
           E1 = max( ear(j-1) , ECHN(i-1) )
           E2 = min( ear(j)   , ECHN(i)   )
           frac = ( E2-E1 ) / ( ECHN(i)-ECHN(i-1) )
           bs(j) = bs(j) + frac * bkgrate(i)
        end if
     end do
     E = 0.5 * ( ear(j) + ear(j-1) )
  end do
  
! Calculate reference band power (in units of *absolute rms^2*)
  Pr = pow * getcountrate(Elo,Ehi,nex,earx,rephotarx)
! Calculate reference band Poisson noise (in *absolutem rms^2)
  mur = getcountrate(Elo,Ehi,nex,earx,photarx)
  Prnoise = 2.0 * ( br + mur )
  write(*,*)"br,mur=",br,mur
  write(*,*)"Pr (fractional rms)^2/Hz",Pr/mur**2
  
! Open file to write the lag simulation to
  
  ! write(*,*)"Enter root name of simulation products"
  ! read(*,'(a)')root
  write(*,*)"The simulation products have the root name reltrans_sim"
  write(*,*)"(change the name of the products if you do not want to overwrite them with the next simulation)"
  root = 'reltrans_sim'
  lagfile = trim(root) // '.dat'
  flxlagfile = 'x' // trim(root) // '.dat'
  phalagfile = 'x' // trim(root) // '.pha'
  rsplagfile = 'x' // trim(root) // '.rsp'
  status = 0
  call ftgiou(xunit,status)
  open(xunit,file=flxlagfile)
  call ftgiou(unit,status)
  open(unit,file=lagfile)
  
! Loop through energy bins
  write(unit,*)"skip on"
  write(unit,*)"read serr 1 2"
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     mus = getcountrate(ear(i-1),ear(i),nex,earx,photarx)
     Psnoise = 2.0 * ( mus + bs(i) )
     ReG = getcountrate(ear(i-1),ear(i),nex,earx,rephotarx)
     ImG = getcountrate(ear(i-1),ear(i),nex,earx,imphotarx)
     G2  = pow**2 * ( ReG**2 + ImG**2 )
     ! Can finally calculate error     
     dlag(i) = 1.0 + Prnoise/Pr
     dlag(i) = dlag(i) * ( G2*(1.0-gammac2) + Psnoise*Pr )
     dlag(i) = dlag(i) / ( gammac2*G2 )
     dlag(i) = dlag(i) / ( 2.0 * Texp * (fhi-flo) )
     dlag(i) = sqrt( dlag(i) )
     dlag(i) = dlag(i) / ( 2.0 * pi * fc )
     !Now generate simulated data
     lagsim(i) = lag(i) + gasdev(idum) * dlag(i)
     !Write out
     write(unit,*)E,0.5*dE,lagsim(i),dlag(i),lag(i)
     write(xunit,*)ear(i-1),ear(i),dE*lagsim(i),dE*dlag(i)
  end do
  close(unit)
  call ftfiou(unit,status)
  close(xunit)
  call ftfiou(xunit,status)
 
  command = 'flx2xsp ' // trim(flxlagfile) // ' ' // trim(phalagfile)
  command = trim(command) // ' ' // trim(rsplagfile)
  write(*,*)"-----------------------------------------------"
  write(*,*)"Outputs: ",trim(lagfile),', ',trim(flxlagfile)
  write(*,*)"command: ",trim(command)
  write(*,*)"-----------------------------------------------"
 
  return
end subroutine simrtdist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine tdreltransCpF(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(22), photar(ne), par(32)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h
  par(2)  = param(2)         !a
  par(3)  = param(3)         !inc
  par(4)  = param(4)         !rin
  par(5)  = param(5)         !rout
  par(6)  = param(6)         !zcos
  par(7)  = param(7)         !Gamma
  par(8)  = param(8)         !logxi
  par(9)  = param(9)         !Afe
  par(10) = param(10)        !lognep
  par(11) = param(11)        !kTe
  par(12) = param(12)        !Nh
  par(13) = param(13)        !boost
  par(14) = 1.0              !qboost
  par(15) = param(14)        !Mass
  par(16) = 0.0              !honr
  par(17) = 0.0              !b1
  par(18) = 0.0              !b2
  par(19) = 0.0              !floHz
  par(20) = 0.0              !fhiHz
  par(21) = 0.0              !ReIm
  par(22) = 0.0              !DelA
  par(23) = param(19)        !DelAB
  par(24) = param(21)        !g
  par(25) = 1.0              !Anorm
  par(26) = param(21)        !telescope response
  par(27) = param(15)        !First lc lowEn band
  par(28) = param(16)        !First lc highEn band
  par(29) = param(17)        !Second lc lowEn band
  par(30) = param(18)        !Second lc highEn band
  par(31) = param(20)        !AB slope 
  par(32) = param(22)        !g slope
! Call general code
  call genreltrans(Cp, dset, ear, ne, par, ifl, photar)
  return
end subroutine tdreltransCpF
!-----------------------------------------------------------------------
    
!-----------------------------------------------------------------------
subroutine genreltrans(Cp, dset, ear, ne, param, ifl, photar)
!
! All reltrans flavours are calculated in this subroutine.
! Cp and dset are the settings:
! |Cp|=1 means use cut-off power-law, |Cp|=2 means use nthcomp
! Cp>1 means there is a density parameter, Cp<1 means density is hardwired  
! dset=0 means ionisation is a parameter, dset=1 means ionization is calculated
! from distance. What to do about ION_ZONES=1 in the distance model?
!
! The parameter array has 26 parameters. No one model actually has 26
! parameters. In each model, some of these parameters are hardwired, but
! the parameters must be sorted into the param(1:26) array for this subroutine.
  
!    Arg:
! 

  !  Internal variables:
  !      constants:
  !         pi: greek pi
  !         rnmax: maximum radius to consider GR effects
  !         nphi, rno: resolution variables, number of pixels on the observer's camera(b and phib)
  !         Emax, Emin: minimum and maximum range of the internal energy grid which is different than the xspec one
  !         dlogf: resolution parameter of the frequency grid
  !         dyn:   limit to check the saved values
  !         ionvar: sets the ionisation variation (1 = w/ ion var; 0 = w/o ion var)
  
  
  use dyn_gr
  use conv_mod
  implicit none
!Constants
  integer         , parameter :: nphi = 200, nro = 200, ionvar = 1 
  real            , parameter :: Emin = 1e-2, Emax = 3e3, dyn = 1e-7
  double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, &
       dlogf = 0.09 !This is a resolution parameter (base 10)
       
!Args:
  integer, intent(inout) :: ifl
  integer, intent(in)    :: Cp, dset, ne
  real   , intent(inout) :: param(32)
  real   , intent(out)   :: photar(ne)
  
!Variables of the subroutine
!initializer
  integer          :: verbose, me, xe
  logical          :: firstcall, needtrans, needconv
  double precision :: d
!Parameters of the model:
  double precision :: h, a, inc, rin, rout, zcos, Gamma, honr, muobs
  real             :: logxi, Afe, lognep, Ecut_obs, Ecut_s, Dkpc, Anorm
  real             :: Nh, boost, Mass, floHz, fhiHz, DelA, DelAB, g
  integer          :: ReIm, resp_matr
  double precision :: qboost,b1,b2
! >>>>>>> distance
!internal frequency grid
  integer          :: nf 
  real             :: f, fac
  double precision :: fc, flo, fhi 
! internal energy grid and xspec energy grid
  real             :: E, dE, dloge
  real             :: earx(0:nex)   
  real             :: ear(0:ne)
!relativistic parameters and limit on rin and h
  real             :: gso, gsd
  double precision :: rmin, rh, lens

!TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
  complex, dimension(:,:,:,:), allocatable :: transe, transea
  real   , dimension(:,:)    , allocatable :: ReW0, ImW0, ReW1, ImW1,&
       ReW2, ImW2, ReW3, ImW3, ReSraw, ImSraw, ReSrawa, ImSrawa,&
       ReGrawa, ImGrawa, ReG, ImG
  double precision :: frobs, frrel  !reflection fraction variables (verbose)
!Radial and angle profile 
  integer                       :: mubin, rbin
  double precision, allocatable :: logxir(:),gsdr(:), logner(:)
!Reflection + total corss spectrum
  real    :: mue, logxi0, reline(nex), imline(nex), photarx(nex), photerx(nex)
  real    :: absorbx(nex),contx(nex), ImGbar(nex), ReGbar(nex)
  real    :: ReGx(nex),ImGx(nex),ReS(ne),ImS(ne)
!variable for non linear effects
  integer ::  DC, ionvariation
  real    :: photarx_1(nex), photarx_2(nex), photarx_delta(nex), &
       reline_a(nex),imline_a(nex),photarx_dlogxi(nex), &
       dlogxi1, dlogxi2, Gamma1, Gamma2, DeltaGamma  
!SAVE 
  integer          :: nfsave, Cpsave
  real             :: paramsave(32)
  double precision :: fhisave, flosave
!Functions
  integer          :: i, j, myenv
  double precision :: disco, dgsofac

! New  
  double precision :: fcons,get_fcons,ell13pt6,lacc,get_lacc
  real             :: Gamma0,logne,Ecut0,thetae,logxiin
  integer          :: Cp_cont

!lag frequency
  real                :: Ea1keV, Ea2keV, Eb1keV, Eb2keV
  real                :: df, gslope, ABslope
  integer             :: fbinx, Ea1, Ea2, Eb1, Eb2
  real, allocatable   :: ReSrawL(:), ImSrawL(:), ReSL(:), ImSL(:), fix(:)
  integer :: resp_matr_a, resp_matr_b

  
  data firstcall /.true./
  data Cpsave/2/
  data nfsave /-1/  
!Save the first call variables
  save firstcall, dloge, earx, me, xe, d, verbose
  save paramsave, fhisave, flosave, nfsave
  save frobs, frrel, lens, Cpsave, needtrans
  save transe, transea, logxir, gsdr, logner
  save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,contx
  save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG
  
  ifl = 1
  !call FNINIT
  
  ! Initialise some parameters 
  call initialiser(firstcall, Emin, Emax, dloge, earx, rnmax, d, needtrans, me, xe, verbose)
  
!Allocate dynamically the array to calculate the trasfer function 
  if (.not. allocated(re1)) allocate(re1(nphi,nro))
  if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
  if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

  call set_param(dset, param, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
     lognep, Ecut_obs, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, Anorm, resp_matr, Ea1keV, Ea2keV, Eb1keV, Eb2keV, &
     ABslope, gslope)
     
  muobs = cos( inc * pi / 180.d0 )

      

  if (Ea1keV .ne. 0.0 .and.  Ea2keV .ne. 0.0) then
     fc = 0.5d0 * ( floHz + fhiHz )
     nf = 1000
     if (.not. allocated(fix)) allocate(fix(0:nf))
     floHz = 0.1
     fhiHz = 500.
! Set frequency array
     do fbinx = 0, nf
        fix(fbinx) = floHz * (fhiHz / floHz)**( real(fbinx) / real(nf) )
     end do
     DC     = 0
     boost  = abs(boost)

  else
  
     !Work out how many frequencies to average over
     fc = 0.5d0 * ( floHz + fhiHz )
     nf = ceiling( log10(fhiHz/floHz) / dlogf )
     if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
        fhiHz = 0.d0
        floHz = 0.d0
        nf    = 1
     end if

!Decide if this is the DC component or not
     if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
        DC     = 1
! ionvar = 0 This is not necessary because in rawS there is the cond 
        g      = 0.0
        DelAB  = 0.0
        DelA   = 0.0
        ReIm   = 1
     else
        DC     = 0
        boost  = abs(boost)
     end if

endif


!Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
  fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass !4.916d-6 * Mass
  flo   = dble(floHz) * 4.92695275718945d-06 * Mass !4.916d-6 * Mass



!this could go into a subroutine 
!Set minimum r (ISCO) and convert rin and h to rg
  if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
  rmin   = disco( a )
  if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
  rh     = 1.d0+sqrt(1.d0-a**2)
  if( h .lt. 0.d0 ) h = abs(h) * rh
  if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
  if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h
  if( rin .lt. rmin )then
     write(*,*)"Warning! rin<ISCO! Set to ISCO"
     rin = rmin
  end if
  if( h .lt. 1.5d0*rh )then
     write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
     h = 1.5d0 * rh
  end if

!Calculate source to observer g-factor and source frame Ecut
  gso    = real( dgsofac(a,h) )
  Ecut_s = real(1.d0+zcos) * Ecut_obs / gso
  if( verbose .gt. 0 )then
     if( abs(Cp) .eq. 1 )then
        write(*,*)"Ecut in source restframe (keV)=",Ecut_s
     else
        write(*,*)"kTe in source restframe (keV)=", Ecut_s
     end if
  end if
  
!Determine if I need to calculate the kernel
  call need_check(Cp,Cpsave,param,paramsave,fhi,flo,fhisave,flosave,nf,nfsave,needtrans,needconv)
   
! Allocate arrays that depend on frequency
  if( nf .ne. nfsave )then
        if( allocated(transe ) ) deallocate(transe )
        if( allocated(transea) ) deallocate(transea)
        allocate(  transe(nex,nf,me,xe) )
        allocate( transea(nex,nf,me,xe) )
        if( allocated(ReW0) ) deallocate(ReW0)
        if( allocated(ImW0) ) deallocate(ImW0)
        if( allocated(ReW1) ) deallocate(ReW1)
        if( allocated(ImW1) ) deallocate(ImW1)
        if( allocated(ReW2) ) deallocate(ReW2)
        if( allocated(ImW2) ) deallocate(ImW2)
        if( allocated(ReW3) ) deallocate(ReW3)
        if( allocated(ImW3) ) deallocate(ImW3)
        allocate( ReW0(nex,nf) )
        allocate( ImW0(nex,nf) )
        allocate( ReW1(nex,nf) )
        allocate( ImW1(nex,nf) )
        allocate( ReW2(nex,nf) )
        allocate( ImW2(nex,nf) )
        allocate( ReW3(nex,nf) )
        allocate( ImW3(nex,nf) )
        if( allocated(ReSraw) ) deallocate(ReSraw)
        if( allocated(ImSraw) ) deallocate(ImSraw)
        allocate( ReSraw(nex,nf) )
        allocate( ImSraw(nex,nf) )
        if( allocated(ReSrawa) ) deallocate(ReSrawa)
        if( allocated(ImSrawa) ) deallocate(ImSrawa)
        allocate( ReSrawa(nex,nf) )
        allocate( ImSrawa(nex,nf) )
        if( allocated(ReGrawa) ) deallocate(ReGrawa)
        if( allocated(ImGrawa) ) deallocate(ImGrawa)
        allocate( ReGrawa(nex,nf) )
        allocate( ImGrawa(nex,nf) )
        if( allocated(ReG) ) deallocate(ReG)
        if( allocated(ImG) ) deallocate(ImG)
        allocate( ReG(nex,nf) )
        allocate( ImG(nex,nf) )

        if( allocated(ReSrawL) ) deallocate(ReSrawL)
        if( allocated(ImSrawL) ) deallocate(ImSrawL)
        allocate( ReSrawL(nf) )
        allocate( ImSrawL(nf) )
        if( allocated(ReSL) ) deallocate(ReSL)
        if( allocated(ImSL) ) deallocate(ImSL)
        allocate( ReSL(ne) )
        allocate( ImSL(ne) )

  end if

! Calculate continuum - so fast there is no need for an if statement
  Cp_cont = Cp
  if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx
  call ad_getcont(nex, earx, Gamma, Afe, Ecut_obs, lognep, logxi, Cp_cont, contx)
     
  if( dset .eq. 1 )then
     fcons = get_fcons(h,a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)
     if( verbose .gt. 0 )then
        lacc = get_lacc(h,a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)
        write(*,*)"Lacc/Ledd=",lacc
        ell13pt6 = fcons * Mass * 1.73152e-28
        write(*,*)"13.6eV-13.6keV luminosity of single source=",ell13pt6
     end if
  else
     fcons = 0.0 !this never gets used for dset=0
     if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
  end if
  
  if( needtrans )then
     !Allocate arrays for kernels     
     if( .not. allocated(logxir) ) allocate( logxir(xe) )
     if( .not. allocated(gsdr)   ) allocate( gsdr  (xe) )
     if( .not. allocated(logner) ) allocate( logner(xe) )
     !Calculate the Kernel for the given parameters
     status_re_tau = .true.
     call rtrans(a, h, muobs, Gamma, rin, rout, honr, d, rnmax, zcos,&
          b1,b2,qboost,dset,fcons,nro, nphi, nex, dloge, nf, fhi, flo, me, xe, &
          logxi, lognep, transe, transea, frobs, frrel, &
          lens, logxir, gsdr, logner)
  end if
  
  if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",boost*frobs
  if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

  
  if( needconv )then
     !needtrans = .false.
     !Initialize arrays for transfer functions
     ReW0 = 0.0
     ImW0 = 0.0
     ReW1 = 0.0
     ImW1 = 0.0
     ReW2 = 0.0
     ImW2 = 0.0
     ReW3 = 0.0
     ImW3 = 0.0
     DeltaGamma = 0.01
     Gamma1 = real(Gamma) - 0.5*DeltaGamma
     Gamma2 = real(Gamma) + 0.5*DeltaGamma
     
     !Get logxi values corresponding to Gamma1 and Gamma2
     call xilimits(nex,earx,contx,DeltaGamma,gso,real(zcos),dlogxi1,dlogxi2)

!Set the ion-variation to 1, there is an if inside the radial loop to check if either the ionvar is 0 or the logxi is 0 to set ionvariation to 0
! it is important that ionvariation is different than ionvar because ionvar is used also later in rawS routine to calculate the cross-spectrum
     ionvariation = 1

     !Loop over radius, emission angle and frequency
     do rbin = 1, xe  !Loop over radial zones
        !Set parameters with radial dependence
        Gamma0 = real(Gamma)
        logne  = logner(rbin)
        Ecut0  = real( gsdr(rbin) ) * Ecut_s
        logxi0 = real( logxir(rbin) )
        if( xe .eq. 1 )then
           Ecut0  = Ecut_s
           logne  = lognep
           logxi0 = logxi
        end if
        !Avoid negative values of the ionisation parameter 
        if (logxi0 .eq. 0.0 .or. ionvar .eq. 0) then
           ionvariation = 0.0
        end if
        do mubin = 1, me      !loop over emission angle zones
           !Calculate input emission angle
           mue    = ( real(mubin) - 0.5 ) / real(me)
           thetae = acos( mue ) * 180.0 / real(pi)
           if( me .eq. 1 ) thetae = real(inc)
           !Call restframe reflection model
           call myreflect(earx,nex,Gamma0,Afe,logne,Ecut0,logxi0,thetae,Cp,photarx)
           !NON LINEAR EFFECTS
           if (DC .eq. 0) then 
              !Gamma variations
              logxiin = logxi0 + ionvariation * dlogxi1
              call myreflect(earx,nex,Gamma1,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_1)
              logxiin = logxi0 + ionvariation * dlogxi2
              call myreflect(earx,nex,Gamma2,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_2)
              photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
              !xi variations
              logxiin = logxi0 + ionvariation * dlogxi1
              call myreflect(earx,nex,Gamma0,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_1)
              logxiin = logxi0 + ionvariation * dlogxi2
              call myreflect(earx,nex,Gamma0,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_2)
              photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10           
           end if


           !Loop through frequencies
           do j = 1,nf
                 do i = 1,nex
                    reline(i)   = real(  transe(i,j,mubin,rbin) )
                    imline(i)   = aimag( transe(i,j,mubin,rbin) )
                    reline_a(i) = real(  transea(i,j,mubin,rbin) )
                    imline_a(i) = aimag( transea(i,j,mubin,rbin) )
                 end do
                 call conv_all_FFTw(dyn, photarx, photarx_delta, reline, imline, reline_a , imline_a,&
                         photarx_dlogxi, ReW0(:,j), ImW0(:,j), ReW1(:,j), ImW1(:,j), &
                         ReW2(:,j), ImW2(:,j), ReW3(:,j), ImW3(:,j), DC)
           end do
        end do
     end do

  end if
  
  !call absorption model 
  call tbabs(earx,nex,nh,Ifl,absorbx,photerx)

  ! write(*,*) nex, Ea1kev, Ea2kev, Eb1kev, Eb2kev
  !Calculate the lag between the two bands
  if (Ea1keV .ne. 0.0 .and.  Ea2keV .ne. 0.0) then
     ! !Figure out the 2 energy bands in terms of array index
     ! Ea1 = ceiling( real(nex) * log10(Ea1kev / Emin) / log10(Emax / Emin))
     ! Ea2 = ceiling( real(nex) * log10(Ea2kev / Emin) / log10(Emax / Emin))
     ! Eb1 = ceiling( real(nex) * log10(Eb1kev / Emin) / log10(Emax / Emin))
     ! Eb2 = ceiling( real(nex) * log10(Eb2kev / Emin) / log10(Emax / Emin))
     ! ! write(*,*) Ea1, Ea2, Eb1, Eb2
     ! ! write(*,*) earx(Ea1), earx(Ea2), earx(Eb1), earx(Eb2)
     ! call lag_freq(nex, earx, Ea1, Ea2, Eb1, Eb2, contx,&
     !      ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3, ImW3, &
     !      absorbx, g, gslope, DelAB, ABslope, nf, fix, boost, real(zcos), gso, &
     !      real(lens), real(Gamma), ionvar, ReSrawL, ImSrawL)

     resp_matr_a = -1
     resp_matr_b = -1
     call lag_freq_resp(nex,earx,Emin,dloge,Ea1keV,Ea2keV,Eb1keV,Eb2keV,resp_matr_a,resp_matr_b,&
     contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
     absorbx, g, gslope,  DelAB, ABslope, nfx, fix, boost, z, gso, lens,&
     Gamma, ionvar, ReGraw, ImGraw)
     
     !Rebin for xspec
     call rebinE(fix, ReSrawL, nf, ear, ReSL, ne)
     call rebinE(fix, ImSrawL, nf, ear, ImSL, ne)

     do j = 1, ne
        df = ear(j) - ear(j - 1)
        photar(j) = atan2( ImSL(j) , ReSL(j) ) / ( pi * (ear(j) + ear(j - 1)) ) * df! times 2.0 and 0.5 they cancel out
     end do

  else 
     
! Calculate raw FT of the full spectrum without absorption
     call rawS(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,real(zcos),&
          gso,real(lens),real(Gamma),ionvar,DC,ReSraw,ImSraw)

! Multiply the raw FT by the absorption 
     do j = 1, nf
        do i = 1, nex
           ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
           ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
        end do
     end do

     if( DC .eq. 1 )then
        do i = 1, nex
           ReGbar(i) = Anorm * ReSrawa(i,1)
           !Norm is applied internally for DC component of dset=1
           !        ImGbar(i) = ImSrawa(i,1)  !No need for the immaginary part in DC
        end do
     else
        ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
        if (ReIm .gt. 0.0) then
           call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa, resp_matr)
        else
           call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
        endif

        ! Apply phase correction parameter to the cross-spectral model (for bad calibration)
        do j = 1,nf
           do i = 1,nex
              ReG(i,j) = cos(DelA) * ReGrawa(i,j) - sin(DelA) * ImGrawa(i,j)
              ImG(i,j) = cos(DelA) * ImGrawa(i,j) + sin(DelA) * ReGrawa(i,j)
           end do
        end do

        ReGbar = 0.0
        ImGbar = 0.0
        fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
        do j = 1,nf
           f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
           do i = 1,nex
              ReGbar(i) = ReGbar(i) + ReG(i,j) / f
              ImGbar(i) = ImGbar(i) + ImG(i,j) / f
           end do
        end do
        ReGbar = ReGbar * fac * Anorm**2  !This means that norm for the AC
        ImGbar = ImGbar * fac * Anorm**2  !components in the dset=1 model
        !is power in squared fractional rms format
     end if

     ! Write output depending on ReIm parameter
     !  if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) ) ReIm = 1
     if( abs(ReIm) .le. 4 )then
        call crebin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is in photar form
        if( abs(ReIm) .eq. 1 )then        !Real part
           photar = ReS
        else if( abs(ReIm) .eq. 2 )then   !Imaginary part
           photar = ImS
        else if( abs(ReIm) .eq. 3 )then   !Modulus
           photar = sqrt( ReS**2 + ImS**2 )
           write(*,*) "Warning ReIm=3 should not be used for fitting!"
        else if( abs(ReIm) .eq. 4 )then   !Time lag (s)
           do i = 1,ne
              dE = ear(i) - ear(i-1)
              photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
           end do
           write(*,*)"Warning ReIm=4 should not be used for fitting!"
        end if
     else
        call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS,resp_matr) !S is count rate in energy channel (photar)
        if( abs(ReIm) .eq. 5 )then        !Modulus
           do i = 1, ne
              dE = ear(i) - ear(i-1)
              photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 )
           end do
        else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
           do i = 1, ne
              dE = ear(i) - ear(i-1)
              photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
           end do
        end if
     end if

     if (verbose .gt. 1 .and. abs(ReIm) .gt. 0) then
        !this writes the individual components to file
        call write_ComponentS(ne,ear,nex,earx,nf,contx,absorbx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,floHz,fhiHz,&
             ReIm,DelA,DelAB,g,boost,real(zcos),gso,real(lens),real(Gamma),ionvar,resp_matr)
        !this writes the full model as returned to Xspec 
        !note that xspec gets output in e.g. lags*dE, and we want just the lags, so a factor dE needs to be included
        open (unit = 14, file = 'Output/Total.dat', status='replace', action = 'write')
        do i = 1,ne 
           dE = ear(i) - ear(i-1)
           write (14,*) (ear(i)+ear(i-1))/2., photar(i)/dE
        end do
        close(14)  
     endif

     fhisave   = fhi
     flosave   = flo
     nfsave    = nf
     paramsave = param
     Cpsave    = Cp

  endif
  
end subroutine genreltrans
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine xilimits(nex,earx,contx,DeltaGamma,gso,z,dlogxi1,dlogxi2)
! Inputs: nex,earx,contx,DeltaGamma,gso
! Outputs: dlogxi1,dlogxi2
! logxi1 = logxi0 + dlogxi1
  implicit none
  integer nex,i
  real earx(0:nex),contx(nex),DeltaGamma,gso,z,dlogxi1,dlogxi2
  real num1,num2,den,logS1,logS2,gsoz,E
  gsoz = gso / (1.0+z)
  num1 = 0.0
  num2 = 0.0
  den = 0.0
  do i = 1,nex
     E   = 0.5 * ( earx(i) + earx(i-1) )
     num2 = num2 + E**(1.0-0.5*DeltaGamma) * contx(i)
     num1 = num1 + E**(1.0+0.5*DeltaGamma) * contx(i)
     den  = den  + E * contx(i)
  end do
  logS1 = log10(num1/den)
  logS2 = log10(num2/den)
  dlogxi1 = -0.5*DeltaGamma * log10(gsoz) + logS1
  dlogxi2 =  0.5*DeltaGamma * log10(gsoz) + logS2
  return
end subroutine xilimits
!-----------------------------------------------------------------------  
      
