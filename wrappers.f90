! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for \ell * g^{2+Gamma}
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) Ionisation profile
! 5) The dependence on incident angle (mimicked by tweaking the ionization)

!CURRENT BRANCH
! This branch is adding the multiple flavours to the reltrans model

include 'subroutines/amodules.f90'
include 'subroutines/header.h'
          
!-----------------------------------------------------------------------
subroutine tdreltransDCp(ear, ne, param, ifl, photar)
  implicit none
  integer, parameter :: nlp = 1 !use a single lamp post
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
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
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = param(15)        !floHz
  par(24) = param(16)        !fhiHz
  par(25) = param(17)        !ReIm
  par(26) = param(18)        !DelA
  par(27) = param(19)        !DelAB
  par(28) = param(20)        !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = 1.0              !Anorm
  par(32) = param(21)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransDCp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine tdreltransPL(ear, ne, param, ifl, photar)
  implicit none
  integer, parameter :: nlp = 1 !use a single lamp post
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 1   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
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
  par(12) = param(11)        !Ecut
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = param(15)        !floHz
  par(24) = param(16)        !fhiHz
  par(25) = param(17)        !ReIm
  par(26) = param(18)        !DelA
  par(27) = param(19)        !DelAB
  par(28) = param(20)        !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = 1.0              !Anorm
  par(32) = param(21)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransPL
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransx(ear,ne,param,ifl,photar)
  implicit none
  integer, parameter :: nlp = 1 !use a single lamp post
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density
  dset = 0   !dset=0 means distance is not set, logxi set instead
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
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = param(15)        !floHz
  par(24) = param(16)        !fhiHz
  par(25) = param(17)        !ReIm
  par(26) = param(18)        !DelA
  par(27) = param(19)        !DelAB
  par(28) = param(20)        !g
  par(29) = 0.               !DelAB2   
  par(30) = 0.               !g2
  par(31) = 1.0              !Anorm
  par(32) = param(21)        !telescope response
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransx
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransDbl(ear, ne, param, ifl, photar)
  implicit none
  integer, parameter :: nlp = 2 !use a double lamp post 
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(27), photar(ne), par(32)
!Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
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
  par(13) = param(13)        !eta_0
  par(14) = param(14)        !eta 
  par(15) = param(15)        !beta_p
  par(16) = param(16)        !Nh
  par(17) = param(17)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(18)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = param(19)        !floHz
  par(24) = param(20)        !fhiHz
  par(25) = param(21)        !ReIm
  par(26) = param(22)        !DelA
  par(27) = param(23)        !DelAB1
  par(28) = param(24)        !g1
  par(29) = param(25)        !DelAB2
  par(30) = param(26)        !g2
  par(31) = 1.0              !Anorm
  par(32) = param(27)        !resp
! Call general code
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  return
end subroutine tdreltransDbl
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdrtdist(ear, ne, param, ifl, photar)
  implicit none
  integer, parameter :: nlp = 1 !use a single lamp post
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
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
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = 1.0              !boost
  par(18) = param(13)        !qboost
  par(19) = param(14)        !Mass
  par(20) = param(15)        !honr
  par(21) = param(16)        !b1
  par(22) = param(17)        !b2
  par(23) = param(18)        !floHz
  par(24) = param(19)        !fhiHz
  par(25) = param(20)        !ReIm
  par(26) = param(21)        !DelA
  par(27) = param(22)        !DelAB
  par(28) = param(23)        !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(24)        !Anorm
  par(32) = param(25)        !telescope response
! Check that we're not looking at the side of the disc
  honr = par(20)
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
  integer, parameter :: nlp = 1 !use a single lamp post
  integer :: ne, ifl, Cp, dset
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
  double precision    :: honr,pi,cosi,cos0
! Settings
  Cp   = 0   !Cp=0 means use the reflionx model with nthcomp and free density 
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
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
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = 1.0              !boost
  par(18) = param(13)        !qboost
  par(19) = param(14)        !Mass
  par(20) = param(15)        !honr
  par(21) = param(16)        !b1
  par(22) = param(17)        !b2
  par(23) = param(18)        !floHz
  par(24) = param(19)        !fhiHz
  par(25) = param(20)        !ReIm
  par(26) = param(21)        !DelA
  par(27) = param(22)        !DelAB
  par(28) = param(23)        !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(24)        !Anorm
  par(32) = param(25)        !telescope response
! Check that we're not looking at the side of the disc
  honr = par(20)
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
subroutine simrtdbl(ear, ne, param, ifl, photar)
  use telematrix
  use env_variables
  implicit none
  integer :: ne, ifl, Cp, dset, i
  real    :: ear(0:ne), param(28), photar(ne), par(32)
  real    :: gammac2, Texp, E, dE, getcountrate
  real    :: rephotar(ne), imphotar(ne)
  real, parameter :: Emin = 1e-1, Emax = 300.0
  integer, parameter :: nex=2**12
  integer, parameter :: nlp = 2 !use a double lamp post
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
  dset = 0   !dset=1 means distance is set, logxi is calculated internally
! Transfer to general parameter array
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
  par(13) = param(13)        !eta_0
  par(14) = param(14)        !eta
  par(15) = 0.               !beta_p
  par(16) = param(15)        !Nh
  par(17) = param(16)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(17)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = param(18)        !floHz
  par(24) = param(19)        !fhiHz
  gammac2 = param(20)        !squared coherence
  par(26) = param(21)        !DelA
  par(27) = param(22)        !DelAB
  par(28) = param(23)        !g
  par(20) = param(24)        !DelAB2
  par(30) = param(25)        !g2
  par(31) = 1.0              !Anorm
  Texp    = param(26)        !Texp (s)
  pow     = param(27)        !power in [rms/mean]^2/Hz units (alpha(nu))
  par(32) = param(28)        !telescope response
  
  flo = par(23)
  fhi = par(24)
  fc  = 0.5 * ( fhi + flo )
  
! Get `folded' lags
  par(25) = 6.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  do i = 1,ne
     lag(i) = photar(i) / ( ear(i) - ear(i-1) )
  end do
  
! Set internal energy grid
  do i = 0, nex
     earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
  end do
  
! Calculate real and imaginary parts on the fine energy grid
  par(25) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, rephotarx)  
  par(25) = 2.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, imphotarx)  
! Get DC component
  par(23) = 0.0   !floHz
  par(24) = 0.0   !fhiHz
  par(25) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, photarx)  

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
  
  write(*,*)"Enter file name of simulation products"
  read(*,'(a)')root
  lagfile =  trim(root) // '.dat'
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
end subroutine simrtdbl
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
  integer, parameter :: nlp = 1 !use a single lamp post
  real :: earx(0:nex),photarx(nex),pow
  real :: Pr,rephotarx(nex),imphotarx(nex),mur,mus
  real :: dlag(ne),G2,ReG,ImG,Psnoise,Prnoise,br,bs(ne)
  real :: flo,fhi,fc,lag(ne),gasdev,lagsim(ne)
  real, parameter :: pi = acos(-1.0)
  integer unit,xunit,status,j
  real E1,E2,frac
  character (len=200) command,flxlagfile,phalagfile,rsplagfile,lagfile,root
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 1   !dset=1 means distance is set, logxi is calculated internally
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
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = 1.0              !boost
  par(18) = param(13)        !qboost
  par(19) = param(14)        !Mass
  par(20) = param(15)        !honr
  par(21) = param(16)        !b1
  par(22) = param(17)        !b2
  par(23) = param(18)        !floHz
  par(24) = param(19)        !fhiHz
  gammac2 = param(20)        !squared coherence
  par(26) = param(21)        !DelA
  par(27) = param(22)        !DelAB
  par(28) = param(23)        !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(24)        !Anorm
  Texp    = param(25)        !Texp (s)
  pow     = param(26)        !power in [rms/mean]^2/Hz units (alpha(nu))
  par(32) = param(27)        !telescope response
  
  flo = par(23)
  fhi = par(24)
  fc  = 0.5 * ( fhi + flo )
  
! Get `folded' lags
  par(25) = 6.0   !ReIm
  call genreltrans(Cp, dset, nlp, ear, ne, par, ifl, photar)  
  do i = 1,ne
     lag(i) = photar(i) / ( ear(i) - ear(i-1) )
  end do
  
! Set internal energy grid
  do i = 0, nex
     earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
  end do
  
! Calculate real and imaginary parts on the fine energy grid
  par(25) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, rephotarx)  
  par(25) = 2.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, imphotarx)  
! Get DC component
  par(23) = 0.0   !floHz
  par(24) = 0.0   !fhiHz
  par(25) = 1.0   !ReIm
  call genreltrans(Cp, dset, nlp, earx, nex, par, ifl, photarx)  

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
