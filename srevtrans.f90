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
subroutine tdreltrans(ear,ne,param,ifl,photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(20), photar(ne), par(27)
! Settings
  Cp   = -1  !|Cp|=1 means cut-off pl, Cp<1 means no density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10)  = param(9)        !Afe
  par(11) = 15.0             !lognep
  par(12) = param(10)        !Ecut_obs
  par(13) = param(11)        !Nh
  par(14) = param(12)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(13)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(14)        !floHz
  par(21) = param(15)        !fhiHz
  par(22) = param(16)        !ReIm
  par(23) = param(17)        !DelA
  par(24) = param(18)        !DelAB
  par(25) = param(19)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(20)        !telescope response
  ! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)
  return
end subroutine tdreltrans
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine tdreltransD(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(20), photar(ne), par(27)
! Settings
  Cp   = 1   !|Cp|=1 means cut-off pl, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = 300.0            !Ecut_obs
  par(13) = param(11)        !Nh
  par(14) = param(12)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(13)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(14)        !floHz
  par(21) = param(15)        !fhiHz
  par(22) = param(16)        !ReIm
  par(23) = param(17)        !DelA
  par(24) = param(18)        !DelAB
  par(25) = param(19)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(20)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)
  return
end subroutine tdreltransD
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransCp(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(20), photar(ne), par(27)
! Settings
  Cp   = -2  !|Cp|=2 means nthcomp, Cp<1 means no density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = 15.0             !lognep
  par(12) = param(10)        !kTe
  par(13) = param(11)        !Nh
  par(14) = param(12)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(13)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(14)        !floHz
  par(21) = param(15)        !fhiHz
  par(22) = param(16)        !ReIm
  par(23) = param(17)        !DelA
  par(24) = param(18)        !DelAB
  par(25) = param(19)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(20)        !telescope response
  ! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)     
  return
end subroutine tdreltransCp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine tdreltransDCp(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(21), photar(ne), par(27)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = param(12)        !Nh
  par(14) = param(13)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(14)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(15)        !floHz
  par(21) = param(16)        !fhiHz
  par(22) = param(17)        !ReIm
  par(23) = param(18)        !DelA
  par(24) = param(19)        !DelAB
  par(25) = param(20)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(21)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransDCp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransx(ear,ne,param,ifl,photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(21), photar(ne), par(27)
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = param(12)        !Nh
  par(14) = param(13)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(14)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(15)        !floHz
  par(21) = param(16)        !fhiHz
  par(22) = param(17)        !ReIm
  par(23) = param(18)        !DelA
  par(24) = param(19)        !DelAB
  par(25) = param(20)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(21)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransx
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine tdreltransDbl(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(22), photar(ne), par(27)
!Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
  nlp =2     !use a double lamp post 
! Transfer to general parameter array REDO THIS PROPERLY -- EDIT GENRELTRANS AND SET_PARAM FIRST
  par(1)  = param(1)         !h1
  par(2)  = param(2)         !h2
  par(3)  = param(3)         !a
  par(4)  = param(4)         !inc
  par(5)  = param(5)         !rin
  par(6)  = param(6)         !rout
  par(7)  = param(7)         !zcos
  par(8)  = param(8)         !Gamma
  par(9)  = param(9)         !logxi
  par(10) = param(10)        !Afe
  par(11) = param(11)        !lognep
  par(12) = param(12)        !kTe
  par(13) = param(13)        !Nh
  par(14) = param(14)        !boost
  par(15) = 1.0              !qboost
  par(16) = param(15)        !Mass
  par(17) = 0.0              !honr
  par(18) = 0.0              !b1
  par(19) = 0.0              !b2
  par(20) = param(16)        !floHz
  par(21) = param(17)        !fhiHz
  par(22) = param(18)        !ReIm
  par(23) = param(19)        !DelA
  par(24) = param(20)        !DelAB
  par(25) = param(21)        !g
  par(26) = 1.0              !Anorm
  par(27) = param(22)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransDbl
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine tdrtdist(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(25), photar(ne), par(27), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !Dkpc
  par(10)  = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = param(12)        !Nh
  par(14) = 1.0              !boost
  par(15) = param(13)        !qboost
  par(16) = param(14)        !Mass
  par(17) = param(15)        !honr
  par(18) = param(16)        !b1
  par(19) = param(17)        !b2
  par(20) = param(18)        !floHz
  par(21) = param(19)        !fhiHz
  par(22) = param(20)        !ReIm
  par(23) = param(21)        !DelA
  par(24) = param(22)        !DelAB
  par(25) = param(23)        !g
  par(26) = param(24)        !Anorm
  par(27) = param(25)        !telescope response
! Check that we're not looking at the side of the disc
  honr = par(17)
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
     call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  end if
  
  return
end subroutine tdrtdist
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine tdrtdistX(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(25), photar(ne), par(27), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density 
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !Dkpc
  par(10)  = param(9)        !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = param(12)        !Nh
  par(14) = 1.0              !boost
  par(15) = param(13)        !qboost
  par(16) = param(14)        !Mass
  par(17) = param(15)        !honr
  par(18) = param(16)        !b1
  par(19) = param(17)        !b2
  par(20) = param(18)        !floHz
  par(21) = param(19)        !fhiHz
  par(22) = param(20)        !ReIm
  par(23) = param(21)        !DelA
  par(24) = param(22)        !DelAB
  par(25) = param(23)        !g
  par(26) = param(24)        !Anorm
  par(27) = param(25)        !telescope response
! Check that we're not looking at the side of the disc
  honr = par(17)
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
     call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  end if
  
  return
end subroutine tdrtdistX
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine simrtdist(ear, ne, param, ifl, photar)
  use telematrix
  implicit none
  integer :: ne, ifl, Cp, dset, nlp, i
  real    :: ear(0:ne), param(27), photar(ne), par(27)
  real    :: gammac2, Texp, E, dE, getcountrate
  real    :: rephotar(ne), imphotar(ne)
  real, parameter :: Emin = 1e-1, Emax = 300.0
  integer, parameter :: nex=2**12
  real :: earx(0:nex),photarx(nex),pow
  real :: Pr,rephotarx(nex),imphotarx(nex),mur,mus
  real :: dlag(ne),G2,ReG,ImG,Psnoise,Prnoise,br,bs(ne)
  real :: flo,fhi,fc,lag(ne),gasdev,lagsim(ne)
  real, parameter :: pi = acos(-1.0)
  integer idum, unit,xunit,status,j
  real E1,E2,frac
  data idum/-2851043/
  save idum
  character (len=200) command,flxlagfile,phalagfile,rsplagfile,lagfile,root
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
  nlp = 1    !use a single lamp post
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !Dkpc
  par(10)  = param(9)        !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = param(12)        !Nh
  par(14) = 1.0              !boost
  par(15) = param(13)        !qboost
  par(16) = param(14)        !Mass
  par(17) = param(15)        !honr
  par(18) = param(16)        !b1
  par(19) = param(17)        !b2
  par(20) = param(18)        !floHz
  par(21) = param(19)        !fhiHz
  gammac2 = param(20)        !squared coherence
  par(23) = param(21)        !DelA
  par(24) = param(22)        !DelAB
  par(25) = param(23)        !g
  par(26) = param(24)        !Anorm
  Texp    = param(25)        !Texp (s)
  pow     = param(26)        !power in [rms/mean]^2/Hz units (alpha(nu))
  par(27) = param(27)        !telescope response
  
  flo = param(19)
  fhi = param(20)
  fc  = 0.5 * ( fhi + flo )
  
! Get `folded' lags
  par(22) = 6.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  do i = 1,ne
     lag(i) = photar(i) / ( ear(i) - ear(i-1) )
  end do
  
! Set internal energy grid
  do i = 0, nex
     earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
  end do
  
! Calculate real and imaginary parts on the fine energy grid
  par(22) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  par(22) = 2.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
! Get DC component
  par(20) = 0.0   !floHz
  par(21) = 0.0   !fhiHz
  par(22) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  

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
  
  write(*,*)"Enter root name of simulation products"
  read(*,'(a)')root
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
subroutine genreltrans(Cp, dset, nlp, ear, ne, param, ifl, photar)
! All reltrans flavours are calculated in this subroutine.
! Cp and dset are the settings:
! |Cp|=1 means use cut-off power-law, |Cp|=2 means use nthcomp
! Cp>1 means there is a density parameter, Cp<1 means density is hardwired  
! dset=0 means ionisation is a parameter, dset=1 means ionization is calculated
! from distance. What to do about ION_ZONES=1 in the distance model?

! The parameter array has 27 parameters. No one model actually has 27
! parameters. In each model, some of these parameters are hardwired, but
! the parameters must be sorted into the param(1:27) array for this subroutine.
  
!    Arg:
! 
!  Internal variables:
!         constants:
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
    integer, intent(in)    :: Cp, dset, ne, nlp
    real   , intent(inout) :: param(27)
    real   , intent(out)   :: photar(ne)  
    !Variables of the subroutine
    !initializer
    integer          :: verbose, me, xe, m
    logical          :: firstcall, needtrans, needconv
    double precision :: d
    !Parameters of the model:
    double precision :: h(nlp), a, inc, rin, rout, zcos, Gamma, honr, muobs 
    real             :: logxi, Afe, lognep, Ecut_obs, Ecut_s, Dkpc, Anorm
    real             :: Nh, boost, Mass, floHz, fhiHz, DelA, DelAB, g
    integer          :: ReIm, resp_matr
    double precision :: qboost,b1,b2
    double precision :: lumratio
    !internal frequency grid
    integer          :: nf 
    real             :: f, fac
    double precision :: fc, flo, fhi
    ! internal energy grid (nex) and output/xspec (ne) energy grid
    real             :: E, dE, dloge
    real             :: earx(0:nex)   
    real             :: ear(0:ne)
    !relativistic parameters and limit on rin and h
    double precision :: rmin, rh 
    double precision :: gso(nlp),tauso(nlp),cosdelta_obs(nlp),height(nlp),contx_int(nlp)
    !note: lens needs to be allocatable to save it. Unclear why only this quantity is saved though....
    double precision, allocatable :: lens(:)
    !TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
    complex, dimension(:,:,:,:), allocatable :: transe, transea
    real   , dimension(:,:)    , allocatable :: ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                                                ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG                                                
    double precision :: frobs, frrel  !reflection fraction variables (verbose)
    !Radial and angle profile 
    integer                       :: mubin, rbin, ibin
    double precision, allocatable :: logxir(:),gsdr(:), logner(:)
    real    :: contx(nex,nlp)
    real    :: mue, logxi0, reline(nex), imline(nex), photarx(nex), photerx(nex)
    real    :: absorbx(nex), ImGbar(nex), ReGbar(nex)
    real    :: ReGx(nex),ImGx(nex),ReS(ne),ImS(ne)
    !variable for non linear effects
    integer ::  DC, ionvariation
    real    :: photarx_1(nex), photarx_2(nex), photarx_delta(nex), &
               reline_a(nex),imline_a(nex),photarx_dlogxi(nex), &
               dlogxi1, dlogxi2, Gamma1, Gamma2, DeltaGamma  
    !SAVE 
    integer          :: nfsave, Cpsave
    real             :: paramsave(27)
    double precision :: fhisave, flosave
    !Functions
    integer          :: i, j, myenv
    double precision :: disco, dgsofac
    real             :: Eintegrate
    ! New  
    double precision :: fcons,get_fcons,ell13pt6,lacc,get_lacc,contx_temp
    real             :: Gamma0,logne,Ecut0,thetae,logxiin
    integer          :: Cp_cont
 
    data firstcall /.true./
    data Cpsave/2/
    data nfsave /-1/  
    !Save the first call variables
    save firstcall, dloge, earx, me, xe, d, verbose
    save paramsave, fhisave, flosave, nfsave
    save frobs, frrel, Cpsave, needtrans, lens
    save transe, transea, logxir, gsdr, logner
    save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3
    save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG

    ifl = 1
    call FNINIT

    ! Initialise some parameters 
    call initialiser(firstcall, Emin, Emax, dloge, earx, rnmax, d, needtrans, me, xe, verbose)

    !Allocate dynamically the array to calculate the trasfer function 
    if (.not. allocated(re1)) allocate(re1(nphi,nro))
    if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
    if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

    !Note: the two different calls are because for the double lP we set the temperature from the coronal frame(s), but for the single
    !LP we use the temperature in the observer frame
    if (nlp .eq. 1) then
        call set_param(dset, param, nlp, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
        lognep, Ecut_obs, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
        DelA, DelAB, g, Anorm, resp_matr)
    else 
        call set_param(dset, param, nlp, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
        lognep, Ecut_s, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
        DelA, DelAB, g, Anorm, resp_matr) 
    end if 
    !placeholder for different luminosities in each LPs DO BETTER
    lumratio = 1.

    muobs = cos( inc * pi / 180.d0 )

    !Work out how many frequencies to average over
    fc = 0.5d0 * ( floHz + fhiHz )
    nf = ceiling( log10(fhiHz/floHz) / dlogf )
    if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
        fhiHz = 0.d0
        floHz = 0.d0
        nf    = 1
    end if
 
    !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
    fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass !4.916d-6 * Mass
    flo   = dble(floHz) * 4.92695275718945d-06 * Mass !4.916d-6 * Mass

    !Decide if this is the DC component or not
    if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
        DC     = 1
        g      = 0.0
        DelAB  = 0.0
        DelA   = 0.0
        ReIm   = 1
    else
        DC     = 0
        boost  = abs(boost)
    end if
  
    !this could go into a subroutine 
    !Set minimum r (ISCO) and convert rin and h to rg
    if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
    rmin   = disco( a )
    if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
    rh     = 1.d0+sqrt(1.d0-a**2)
    if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
    if( rin .lt. rmin )then
        write(*,*)"Warning! rin<ISCO! Set to ISCO"
        rin = rmin
    end if
    do m=1,nlp 
        if( h(m) .lt. 0.d0 ) h(m) = abs(h(m)) * rh
        if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h(m)
        if( h(m) .lt. 1.5d0*rh )then
            write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
            h(m) = 1.5d0 * rh
        end if 
    end do

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
    end if
  
    !allocate lensing array if necessary
    if( needtrans ) then
        if( .not. allocated(lens) ) allocate( lens(nlp) )
    end if
  
    if (nlp .eq. 1) then 
        gso(1) = real( dgsofac(a,h(1)) ) 
        Ecut_s = real(1.d0+zcos) * Ecut_obs / gso(1)
        call getlens(a,h(1),muobs,lens(1),tauso(1),cosdelta_obs(1))
        if( tauso(1) .ne. tauso(1) ) stop "tauso is NaN"
        Cp_cont = Cp
        if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx
        call ad_getcont(nex, earx, Gamma, Afe, Ecut_obs, lognep, logxi, Cp_cont, contx)        
        if( dset .eq. 1 )then
            fcons = get_fcons(h(1),a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)     
        else
            fcons = 0.0
        end if         
        if( verbose .gt. 0 )then
            if( dset .eq. 1 )then    
                lacc = get_lacc(h(1),a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)
                write(*,*)"Lacc/Ledd=",lacc
                ell13pt6 = fcons * Mass * 1.73152e-28
                write(*,*)"13.6eV-13.6keV luminosity of single source=",ell13pt6
            else
                call sourcelum(nex,earx,contx,real(mass),real(gso(1)),real(Gamma))
            end if   
        end if         
        if( abs(Cp) .eq. 1 )then
            write(*,*)"Ecut in source restframe (keV)=",Ecut_s
        else
            write(*,*)"kTe in source restframe (keV)=", Ecut_s
        end if 
        contx_int = 1. !note: for a single LP we don't need to account for this factor in the ionisation profile, so it's defaulted to 1
        contx = lens(1) * (gso(1)/(real(1.d0+zcos))**Gamma) * contx          
    else 
        do m=1,nlp   
            !here the observed cutoffs are set from the temperature in the source frame   
            gso(m) = real( dgsofac(a,h(m)) )
            call getlens(a,h(m),muobs,lens(m),tauso(m),cosdelta_obs(m))
            if( tauso(m) .ne. tauso(m) ) stop "tauso is NaN"
            Ecut_obs = Ecut_s * gso(m) / real(1.d0+zcos)
            Cp_cont = Cp 
            if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx        
            call ad_getcont(nex, earx, Gamma, Afe, Ecut_obs, lognep, logxi, Cp_cont, contx(:,m))
            if (m .gt. 1) contx(:,m) = lumratio*contx(:,m)  
            !TODO fix this section, calculate luminosities better
            if( verbose .gt. 0 )then
                call sourcelum(nex,earx,contx(:,m),real(mass),real(gso(m)),real(Gamma))
                if( abs(Cp) .eq. 1 )then
                    write(*,*)"Ecut observed from source #", m, "is (keV)=" ,Ecut_obs
                else
                    write(*,*)"kTe observed from source #", m, "is (keV)=" ,Ecut_obs
                end if
            end if  
            contx_int(m) = Eintegrate(0.1,1e3,nex,earx,contx(:,m),dlogE)             
            contx(:,m) = lens(m) * (gso(m)/(real(1.d0+zcos))**Gamma) * contx(:,m)   
        end do  
    end if  
 
    if( needtrans )then
        !Allocate arrays for kernels     
        if( .not. allocated(logxir) ) allocate( logxir(xe) )
        if( .not. allocated(gsdr)   ) allocate( gsdr  (xe) )
        if( .not. allocated(logner) ) allocate( logner(xe) )
        !Calculate the Kernel for the given parameters
        status_re_tau = .true.
        call rtrans(verbose,dset,nlp,a,h,muobs,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,lumratio,&
                    fcons,contx_int,tauso,lens,cosdelta_obs,nro,nphi,nex,dloge,nf,fhi,flo,me,xe,logxi,lognep,&
                    transe,transea,frobs,frrel,logxir,gsdr,logner)         
    end if
  
    !do this for each lamp post, then find some sort of weird average?
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
        call xilimits(nex,earx,nlp,contx,DeltaGamma,real(gso),real(lens),real(zcos),dlogxi1,dlogxi2)
        !Set the ion-variation to 1, there is an if inside the radial loop to check if either the ionvar is 0 or the logxi is 0 to
        !set ionvariation to 0  it is important that ionvariation is different than ionvar because ionvar  is used also later in
        !the rawS subroutine to calculate the cross-spectrum
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

    !Calculate raw FT of the full spectrum without absorption
    call rawS(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,real(zcos),&
                real(gso),real(lens),real(Gamma),ionvar,DC,ReSraw,ImSraw)

    ! Calculate absorption and multiply by the raw FT
    call tbabs(earx,nex,nh,Ifl,absorbx,photerx)

    do j = 1, nf
        do i = 1, nex
            ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
            ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
        end do
    end do

    if( DC .eq. 1 )then
        !Norm is applied internally for DC component of dset=1
        !No need for the immaginary part in DC
        do i = 1, nex
            ReGbar(i) = (Anorm/real(1.+lumratio)) * ReSrawa(i,1)
        end do
    else
        !Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
        if (ReIm .gt. 0.0) then
            call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa, resp_matr)
        else
            call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
        endif
        !Apply phase correction parameter to the cross-spectral model (for bad calibration)
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
        !This means that norm for the AC components in the dset=1 model is power in squared fractional rms format
        ReGbar = ReGbar * fac * Anorm**2  
        ImGbar = ImGbar * fac * Anorm**2  
    end if

    !Write output depending on ReIm parameter
    !if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) ) ReIm = 1
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
        call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS,resp_matr) !S is count rate
        if( abs(ReIm) .eq. 5 )then        !Modulus
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
            end do
        else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
            end do
        end if
    end if
  
    !TBD: benchmark how long this loop takes
    if (verbose .gt. 1 .and. abs(ReIm) .gt. 0) then
        !this writes the individual components to file
        call write_components(ne,ear,nex,earx,nf,nlp,contx,absorbx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,floHz,fhiHz,&
                              ReIm,DelA,DelAB,g,boost,real(zcos),real(gso),real(lens),real(Gamma),ionvar,resp_matr)
        !this writes the full model as returned to Xspec 
        !note that xspec gets output in e.g. lags*dE, and we want just the lags, so a factor dE needs to be included
        open (unit = 14, file = 'Output/Total.dat', status='replace', action = 'write')     
        do i = 1,ne 
            dE = ear(i) - ear(i-1)
            write (14,*) (ear(i)+ear(i-1))/2., photar(i)/dE        
        end do 
        close(14)  
        !print continuum for both single and multiple LPs REDO THIS 
        open (unit = 24, file = 'Output/Continuum_spec.dat', status='replace', action = 'write')
        do i=1,nex
            dE = earx(i) - earx(i-1)
            if( nlp .eq. 1 ) then
                contx_temp = contx(i,1)/dE
            else
                contx_temp = 0.
                do m=1,nlp 
                    contx_temp = contx_temp + contx(i,m)
                end do
                contx_temp =  contx_temp/((1.+lumratio)*dE)      
            end if
            write (24,*) (earx(i)+earx(i-1))/2., contx_temp
        end do
        close(24)
    endif 
  
    fhisave   = fhi
    flosave   = flo
    nfsave    = nf
    paramsave = param
    Cpsave    = Cp
  
end subroutine genreltrans
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine xilimits(nex,earx,nlp,contx,DeltaGamma,gso,lens,z,dlogxi1,dlogxi2)
! Inputs: nex,earx,contx,DeltaGamma,gso
! Outputs: dlogxi1,dlogxi2
! logxi1 = logxi0 + dlogxi1
    implicit none
    integer nex,i,m,nlp
    real earx(0:nex),contx(nex,nlp),contx_sum(nex),DeltaGamma,gso(nlp),lens(nlp),z,dlogxi1,dlogxi2
    real num1,num2,den,logS1,logS2,gsoz,E,gso_avg

    !before calculating the differential of the ionisation, set up the arrays properly depending on the number of lampposts
    !note: if we have multiple LPs we have to a) calculate an effective gso factor by averaging over the lensing factor (ie, how
    !luminous each LP appears given their height) and b) get a total continuum flux by summing over each LPs array.
    !for a single lamp posts there is no need to do all this stuff, but we need to read in the contx/gso factors properly to be
    !able to call the following code in the same way
    if(nlp .eq. 1 ) then
        contx_sum = contx(:,1)
        gso_avg = gso(1)
    else
        contx_sum = 0.
        gso_avg = 0.
        do m=1,nlp
            contx_sum = contx_sum + contx(:,m)
            gso_avg = gso_avg + lens(m)*gso(m)
        end do
        gso_avg = gso_avg/sum(lens)
    end if

    gsoz = gso_avg / (1.0+z)
    num1 = 0.0
    num2 = 0.0
    den = 0.0
    do i = 1,nex
        E   = 0.5 * ( earx(i) + earx(i-1) )
        num2 = num2 + E**(1.0-0.5*DeltaGamma) * contx_sum(i)
        num1 = num1 + E**(1.0+0.5*DeltaGamma) * contx_sum(i)
        den  = den  + E * contx_sum(i)
    end do
    logS1 = log10(num1/den)
    logS2 = log10(num2/den)
    dlogxi1 = -0.5*DeltaGamma * log10(gsoz) + logS1
    dlogxi2 =  0.5*DeltaGamma * log10(gsoz) + logS2
    return
end subroutine xilimits
!-----------------------------------------------------------------------  
      
