! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for \ell * g^{2+Gamma}
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) Ionisation profile
! 5) The dependence on incident angle (mimicked by tweaking the ionization)

!CURRENT BRANCH
! This branch is adding the multiple flavours to the reltrans model
  
include 'subroutines/header.h'
          
!-----------------------------------------------------------------------
subroutine tdreltransDCp(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(21), photar(ne), par(32)
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
subroutine tdreltransD(ear, ne, param, ifl, photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(21), photar(ne), par(32)
! Settings
  Cp   = 1   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
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
end subroutine tdreltransD
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine tdreltransx(ear,ne,param,ifl,photar)
  implicit none
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(21), photar(ne), par(32)
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
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(27), photar(ne), par(32)
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
  integer :: ne, ifl, Cp, dset, nlp
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
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
  honr = par(19)
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
  real    :: ear(0:ne), param(25), photar(ne), par(32), getcountrate
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
  implicit none
  integer :: ne, ifl, Cp, dset, nlp, i
  real    :: ear(0:ne), param(28), photar(ne), par(32)
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
  dset = 0   !dset=1 means distance is set, logxi is calculated internally
  nlp = 2    !use a single lamp post
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
  integer :: ne, ifl, Cp, dset, nlp, i
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
  integer unit,xunit,status,j
  real E1,E2,frac
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
    integer         , parameter :: nphi = 200, nro = 200!, ionvar! = 1 
    real            , parameter :: Emin = 1e-2, Emax = 3e3, dyn = 1e-7
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, &
                                   dlogf = 0.09 !This is a resolution parameter (base 10)       
    !Args:
    integer, intent(inout) :: ifl
    integer, intent(in)    :: Cp, dset, ne, nlp
    real   , intent(inout) :: param(32)
    real   , intent(out)   :: photar(ne)  
    !Variables of the subroutine
    !initializer
    integer          :: verbose, me, xe, m, ionvar, refvar
    logical          :: firstcall, needtrans, needconv
    double precision :: d
    !Parameters of the model:
    double precision :: h(nlp), a, inc, rin, rout, zcos, Gamma, honr, muobs 
    real             :: logxi, Afe, lognep, Ecut_obs, Ecut_s, Dkpc, Anorm, beta_p
    real             :: Nh, boost, Mass, floHz, fhiHz, DelA, DelAB(nlp), g(nlp)
    integer          :: ReIm, resp_matr
    double precision :: qboost,b1,b2, eta, eta_0
    !internal frequency grid
    integer          :: nf 
    real             :: f, fac
    double precision :: fc, flo, fhi
    ! internal energy grid (nex) and output/xspec (ne) energy grid
    real             :: E, dE, dloge
    real             :: earx(0:nex)   
    real             :: ear(0:ne)
    ! internal frequency grid, for when we do lag/frequency spectra
    integer           :: fbinx 
    real, allocatable :: fix(:)
    !relativistic parameters and limit on rin and h
    double precision :: rmin, rh 
    double precision :: gso(nlp),tauso(nlp),cosdelta_obs(nlp),height(nlp),contx_int(nlp)
    !lens needs to be allocatable to save it. 
    double precision, allocatable :: lens(:),frobs(:),frrel(:)
    !TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
   ! complex, dimension(:,:,:,:,:), allocatable :: transe, transea
    complex, dimension(:,:,:,:,:), allocatable :: ker_W0,ker_W1,ker_W2,ker_W3
    real   , dimension(:,:,:)    , allocatable :: ReW0,ImW0,ReW1,ImW1
    real   , dimension(:,:,:)    , allocatable :: ReW2,ImW2,ReW3,ImW3
    real   , dimension(:,:)      , allocatable :: ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG                                                
    !double precision :: frobs(nlp), frrel(nlp)  !reflection fraction variables (verbose)
    !Radial and angle profile 
    integer                       :: mubin, rbin, ibin
    double precision, allocatable :: logxir(:),gsdr(:), logner(:)
    real    :: contx(nex,nlp)
    real    :: mue, logxi0, reline_w0(nlp,nex), imline_w0(nlp,nex), photarx(nex), photerx(nex)
    real    :: absorbx(nex), ImGbar(nex), ReGbar(nex)
    real    :: ReGx(nex),ImGx(nex),ReS(ne),ImS(ne)
    !variable for non linear effects
    integer ::  DC, ionvariation
    real    :: photarx_1(nex), photarx_2(nex), photarx_delta(nex), photarx_dlogxi(nex)
    real    :: reline_w1(nlp,nex),imline_w1(nlp,nex),reline_w2(nlp,nex),imline_w2(nlp,nex)
    real    :: reline_w3(nlp,nex),imline_w3(nlp,nex)
    real    :: dlogxi1, dlogxi2, Gamma1, Gamma2, DeltaGamma  
    !SAVE 
    integer          :: nfsave, Cpsave
    real             :: paramsave(32)
    double precision :: fhisave, flosave
    !Functions
    integer          :: i, j, myenv
    double precision :: disco, dgsofac
    ! New  
    double precision :: fcons,get_fcons,contx_temp!,ell13pt6,lacc,get_lacc,
    real             :: Gamma0,logne,Ecut0,thetae,logxiin
    integer          :: Cp_cont
    real time_start,time_end        !runtime stuff
 
    data firstcall /.true./
    data Cpsave/2/
    data nfsave /-1/  
    !Save the first call variables
    save firstcall, dloge, earx, me, xe, d, verbose
    save paramsave, fhisave, flosave, nfsave, refvar, ionvar
    save frobs, frrel, Cpsave, needtrans, lens
    save ker_W0, ker_W1, ker_W2, ker_W3, logxir, gsdr, logner
    save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3
    save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG

    ifl = 1
    !call FNINIT

    ! Initialise some parameters 
    call initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans,me,xe,refvar,ionvar,verbose)

    !Allocate dynamically the array to calculate the trasfer function 
    if (.not. allocated(re1)) allocate(re1(nphi,nro))
    if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
    if (.not. allocated(pem1)) allocate(pem1(nphi,nro))
    
    !Note: the two different calls are because for the double lP we set the temperature from the coronal frame(s), but for the single
    !LP we use the temperature in the observer frame
    if (nlp .eq. 1) then
        call set_param(dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut_obs,&
                       eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                       g,Anorm,resp_matr,refvar,verbose)        
    else 
        call set_param(dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut_s,&
                       eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                       g,Anorm,resp_matr,refvar,verbose) 
    end if 
    
    muobs = cos( inc * pi / 180.d0 )

    !this needs to go in a subroutine - model_mode or something
    !rework this logic so that low frequencies always result in time independent spectrum, not just for reim<7
    if( ReIm .eq. 7 .and. fhiHz .gt. tiny(fhiHz) .and. floHz .gt. tiny(floHz)) then
        !set up frequency grid in Hz if using lag frequency mode, depending on whether we're looking at AGN or XRBs
        if( Mass .gt. 1000 ) then
            floHz = 1.e-5
            fhiHz = 5.e-2 
        else
            floHz = 0.07
            fhiHz = 700.
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
        !Note that the frequency grid is using a higher resolution since it's what we care about in this mode 
        !TBD: find a way to reduce energy resolution since we don't need it, otherwise the runtime is insanely slow  
        nf = ceiling( log10(fhiHz/floHz) / 0.01 )
        allocate(fix(0:nf))
        do fbinx = 0, nf 
            fix(fbinx) = floHz * (fhiHz / floHz)**( real(fbinx) / real(nf) )
        end do
    else 
        !if doing lag-energy spectra, just work out how many frequencies to average over 
        fc = 0.5d0 * ( floHz + fhiHz )
        nf = ceiling( log10(fhiHz/floHz) / dlogf )
        if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
            fhiHz = 0.d0
            floHz = 0.d0
            nf    = 1
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
    end if   

    !Decide if this is the DC component/time averaged spectrum or not
    if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
        DC     = 1
        g      = 0.0
        DelAB  = 0.0
        DelA   = 0.0
        ReIm   = 1
        eta    = eta_0
        beta_p = 1. !this is an ugly hack for the double LP model, to calculate the time-averaged spectrum
    else
        DC     = 0
        boost  = abs(boost)
    end if
    !this could go into a subroutine -- just put it in set_params?
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

    !Determine if I need to calculate the kernel - double check which parameters cause a new kernel (e.g. new eta)
    call need_check(Cp,Cpsave,param,paramsave,fhi,flo,fhisave,flosave,nf,nfsave,needtrans,needconv)

    ! Allocate arrays that depend on frequency
    if( nf .ne. nfsave )then
        if( allocated(ker_W0 ) ) deallocate(ker_W0 )
        if( allocated(ker_W1) ) deallocate(ker_W1 )
        if( allocated(ker_W2 ) ) deallocate(ker_W2 )
        if( allocated(ker_W3 ) ) deallocate(ker_W3 )
        allocate( ker_W0(nlp,nex,nf,me,xe) )
        allocate( ker_W1(nlp,nex,nf,me,xe) )
        allocate( ker_W2(nlp,nex,nf,me,xe) )
        allocate( ker_W3(nlp,nex,nf,me,xe) )
        if( allocated(ReW0) ) deallocate(ReW0)
        if( allocated(ImW0) ) deallocate(ImW0)
        if( allocated(ReW1) ) deallocate(ReW1)
        if( allocated(ImW1) ) deallocate(ImW1)
        if( allocated(ReW2) ) deallocate(ReW2)
        if( allocated(ImW2) ) deallocate(ImW2)
        if( allocated(ReW3) ) deallocate(ReW3)
        if( allocated(ImW3) ) deallocate(ImW3)
        allocate( ReW0(nlp,nex,nf) )
        allocate( ImW0(nlp,nex,nf) )
        allocate( ReW1(nlp,nex,nf) )
        allocate( ImW1(nlp,nex,nf) )
        allocate( ReW2(nlp,nex,nf) )
        allocate( ImW2(nlp,nex,nf) )
        allocate( ReW3(nlp,nex,nf) )
        allocate( ImW3(nlp,nex,nf) )
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
  
    !allocate lensing/reflection fraction arrays if necessary
    if( needtrans ) then
        if( allocated(lens) ) deallocate( lens )
        allocate (lens(nlp))
        if( allocated(frobs) ) deallocate( frobs )
        allocate (frobs(nlp))
        if( allocated(frrel) ) deallocate( frrel )
        allocate (frrel(nlp))
    end if

    !set up the continuum spectrum plus relative quantities (cutoff energies, lensing/gfactors, luminosity, etc)
    call init_cont(nlp,a,h,zcos,Ecut_s,Ecut_obs,gso,muobs,lens,tauso,cosdelta_obs,Cp_cont,Cp,fcons,Gamma,&
                   Dkpc,Mass,earx,Emin,Emax,contx,dlogE,verbose,dset,Anorm,contx_int,eta)    

    if (verbose .gt. 2) call CPU_TIME (time_start)
    if( needtrans )then
        !Allocate arrays for kernels   
        if (allocated (logxir)) deallocate (logxir)
        allocate (logxir(xe))
        if (allocated (gsdr)) deallocate (gsdr)
        allocate (gsdr(xe))
        if (allocated (logner)) deallocate (logner)
        allocate (logner(xe))
        !Calculate the Kernel for the given parameters
        status_re_tau = .true.
        call rtrans(verbose,dset,nlp,a,h,gso,muobs,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,eta_0,&
                    fcons,contx_int,tauso,lens,cosdelta_obs,nro,nphi,nex,dloge,nf,fhi,flo,me,xe,logxi,lognep,&
                    ker_W0,ker_W1,ker_W2,ker_W3,frobs,frrel,logxir,gsdr,logner)         
    end if
    if( verbose .gt. 2 ) then
        call CPU_TIME (time_end)
        print *, 'Transfer function runtime: ', time_end - time_start, ' seconds'
    end if
  
    !do this for each lamp post, then find some sort of weird average?
    if( verbose .gt. 0) write(*,*)"Observer's reflection fraction for each source:",boost*frobs
    if( verbose .gt. 0) write(*,*)"Relxill reflection fraction for each source:",frrel    
    
    if( verbose .gt. 2) call CPU_TIME (time_start)  
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
                !Loop through frequencies and lamp posts
                do j = 1,nf
                    do i = 1,nex
                        do m=1,nlp
                            reline_w0(m,i) = real( ker_W0(m,i,j,mubin,rbin) )
                            imline_w0(m,i) = aimag( ker_W0(m,i,j,mubin,rbin) )
                            reline_w1(m,i) = real( ker_W1(m,i,j,mubin,rbin) )
                            imline_w1(m,i) = aimag( ker_W1(m,i,j,mubin,rbin) )
                            reline_w2(m,i) = real( ker_W2(m,i,j,mubin,rbin) )
                            imline_w2(m,i) = aimag( ker_W2(m,i,j,mubin,rbin) )
                            reline_w3(m,i) = real( ker_W3(m,i,j,mubin,rbin) )
                            imline_w3(m,i) = aimag( ker_W3(m,i,j,mubin,rbin) )
                        end do  
                    end do
                    !always: convolution for reverberation/DC spectrum
                    !TBD: add flag here to do this convolution if no reflection time, or different convolution with complex
                    !xillver if tref > 0 or something.
                    call conv_one_FFTw(dyn,photarx,reline_w0,imline_w0,ReW0(:,:,j),ImW0(:,:,j),DC,nlp)
                    if(DC .eq. 0 .and. refvar .eq. 1) then
                        call conv_one_FFTw(dyn,photarx,reline_w1,imline_w1,ReW1(:,:,j),ImW1(:,:,j),DC,nlp)
                        call conv_one_FFTw(dyn,photarx_delta,reline_w2,imline_w2,ReW2(:,:,j),ImW2(:,:,j),DC,nlp)
                    end if
                    if(DC .eq. 0 .and. ionvar .eq. 1) then
                        call conv_one_FFTw(dyn,photarx_dlogxi,reline_w3,imline_w3,ReW3(:,:,j),ImW3(:,:,j),DC,nlp)
                    end if
                    !old call: always convolve every single transfer function in one go
                    !call conv_all_FFTw(dyn,photarx,photarx_delta,photarx_dlogxi,reline_w0,imline_w0,reline_w1,imline_w1,&
                    !     reline_w2,imline_w2,reline_w3,imline_w3,ReW0(:,:,j),ImW0(:,:,j),ReW1(:,:,j),ImW1(:,:,j),&
                    !     ReW2(:,:,j),ImW2(:,:,j),ReW3(:,:,j),ImW3(:,:,j),DC,nlp)                        
                end do
            end do
        end do
    end if
    if( verbose .gt. 2 ) then
        call CPU_TIME (time_end)
        print *, 'Convolutions runtime: ', time_end - time_start, ' seconds' 
    endif
  
    ! Calculate absorption 
    call tbabs(earx,nex,nh,Ifl,absorbx,photerx)

    !TBD coherence check - if zero coherence between lamp posts, call a different subroutine 
    if( ReIm .eq. 7 ) then
        !tbd - implement zero cohernece in lag_freq
        if(nlp .gt. 1 .and. beta_p .eq. 0. ) then
            call lag_freq_nocoh(nex,earx,nf,fix,real(flo),real(fhi),Emin,Emax,nlp,contx,absorbx,real(tauso),real(gso),&
                                ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),boost,&
                                g,DelAB,ionvar,ReGbar,ImGbar)
        else
            call lag_freq(nex,earx,nf,fix,real(flo),real(fhi),Emin,Emax,nlp,contx,absorbx,real(tauso),real(gso),&
                          ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),beta_p,&
                          boost,g,DelAB,ionvar,ReGbar,ImGbar)        
        end if
    else if (nlp .gt. 1 .and. beta_p .eq. 0.) then
        call rawG(nex,earx,nf,real(flo),real(fhi),nlp,contx,absorbx,real(tauso),real(gso),ReW0,ImW0,&
                  ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),boost,ReIm,g,DelAB,&
                  ionvar,DC,resp_matr,ReGrawa,ImGrawa)                
    else
        !Calculate raw FT of the full spectrum without absorption
        call rawS(nex,earx,nf,real(flo),real(fhi),nlp,contx,real(tauso),real(gso),ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                  real(h),real(zcos),real(Gamma),real(eta),beta_p,boost,g,DelAB,ionvar,DC,ReSraw,ImSraw)
        !Include absorption in the model
        do j = 1, nf
            do i = 1, nex
                ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
                ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
            end do
        end do        
    end if

    if( DC .eq. 1 )then
        !Norm is applied internally for DC/time averaged spectrum component of dset=1
        !No need for the immaginary part in DC
        do i = 1, nex
            ReGbar(i) = (Anorm/real(1.+eta)) * ReSrawa(i,1)
        end do
    else if (ReIm .eq. 7) then     
        !if calculating the lag-frequency spectrum, just rebin the arrays 
        call rebinE(fix, ReGbar, nf, ear, ReS, ne)
        call rebinE(fix, ImGbar, nf, ear, ImS, ne)     
    else 
        !In this case, calculate the lag-energy spectrum
        !Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
        !note: this must be done by rawG for two incoherent lamp posts, hence the skip below
        if(nlp .eq. 1 .or. beta_p .ne. 0.) then
            if (ReIm .gt. 0.0) then
                call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa, resp_matr)
            else
                call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
            endif
        end if
        !Apply phase correction parameter to the cross-spectral model (for bad calibration)
        !this is where coherence = 0 or = 1 cases merge back 
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
        !note: the factor eta is to have the same normalization as the single LP model, it's 100% arbitrary
        ReGbar = ReGbar * fac * (Anorm/real(1.+eta))**2  
        ImGbar = ImGbar * fac * (Anorm/real(1.+eta))**2  
    end if

    !Write output depending on ReIm parameter
    if( ReIm .eq. 7 ) then
        do i=1,ne 
            dE = ear(i) - ear(i-1) 
            photar(i) = atan2(ImS(i),ReS(i))/(pi*(ear(i) + ear(i-1)))*dE
        end do
    else if( abs(ReIm) .le. 4 )then
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
                photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) 
            end do
        else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
            end do
        end if
    end if
    
    if (verbose .gt. 1 .and. abs(ReIm) .gt. 0 .and. ReIm .lt. 7) then
        if( DC .eq. 0 .and. beta_p .eq. 0) then
            call write_components(ne,ear,nex,earx,nf,real(flo),real(fhi),nlp,contx,absorbx,real(tauso),real(gso),&
                                  ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),&
                                  beta_p,boost,floHz,fhiHz,ReIm,DelA,DelAB,g,ionvar,resp_matr)
        !catch case here for coherence = 0 or 1
        end if                
        !this writes the full model as returned to Xspec 
        !note that xspec gets output in e.g. lags*dE, and we want just the lags, so a factor dE needs to be included
        !add writing of components for lag frequency spectrum 
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
                contx_temp =  contx_temp/((1.+eta)*dE)      
            end if
            write (24,*) (earx(i)+earx(i-1))/2., contx_temp
        end do
        close(24)
    else if (ReIm .eq. 7) then
        open (unit = 14, file = 'Output/Total.dat', status='replace', action = 'write')     
        do i = 1,ne 
            dE = ear(i) - ear(i-1)
            write (14,*)  (ear(i)+ear(i-1))/2., photar(i)/dE        
        end do 
        close(14)
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
      
