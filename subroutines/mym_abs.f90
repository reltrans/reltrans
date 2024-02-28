
! etable file on my system:
!filenm = '/home/njm442/Documents/xstar/xstar2xspec/erg/c_fe/fe_5/xout_etable.fits'

!=======================================================================
subroutine mym(ear,ne,param,ifl,photar)
! Calculates transmission fraction from a user-defined etable that
! has logxi as the only free parameter  
  implicit none
  integer ne,ifl,i
  real ear(0:ne),param(4),photar(ne),photer(ne)
  real tau(ne),Nh,logxi,fcov,z,par(2),Nh0,trans(ne),Emax0
  real ear_rest(0:ne),Emin0
  character (len=255) filenm
  logical firstcall
  save filenm,Emax0,Nh0,firstcall,Emin0
!  save Emax0,Nh0,firstcall,Emin0
  data firstcall/.true./
!  filenm = '/home/njm442/Documents/xstar/xstar2xspec/erg/c_fe/fe_5/xout_etable.fits'
! Read in parameters
  Nh    = param(1)
  logxi = param(2)
  fcov  = param(3)
  z     = param(4)
  
! Read into mtable parameter array
  par(1) = logxi
  par(2) = 0.0  !Call mtable with z=0 because the redshift functionality is unreliable

! Adjust energy array by redshift
  ear_rest = ear * (1.0+z)
  
! Initialize
  if(firstcall)then
     !Get the name of the etable file
     write(*,*)"Enter name of etable (with full path)"
     read(*,'(a)')filenm

     !Read the largest energy in the grid
!     call getEmax(filenm,Emax0)
     !Read the smallest energy in the grid
!     call getEmin(filenm,Emin0)
     emax0=20.0
     emin0=0.1
     !Read column used from the etable and convert to units of 1e22
!     write(*,*)"Ive just set the energy boundary values..."
     call getNh(filenm,Nh0)
!     write(*,*)"Ive just called getNh..."
     Nh0 = Nh0/1e22
!     write(*,*)"Ive just normalised Nh..."
     !Set firstcall
     firstcall = .false.
  end if
 
!  write(*,*)"About to call mtable subroutine..."   
  
!  write(*,*)filenm
! Call mtable routine to interpolate from grid
  call xsmtbl(ear_rest,ne,par,filenm,ifl,tau,photer)
  
!  write(*,*)"I've just called mtable subroutine..." 
  
! Extrapolate above the highest energy in the grid
  call extrap(ear,ne,z,filenm,Emax0,tau)
 
!  write(*,*)"I've just called extrap subroutine..." 

! Extrapolate below the lowest energy in the grid
  call extraplo(ear,ne,z,filenm,Emin0,tau)
 
!  write(*,*)"I've just called extraplo subroutine..." 
 
! Convert tau for the input Nh
  tau = tau * Nh/Nh0
  
!  write(*,*)"I've just coverted tau to the input Nh..." 

! Calculate transmission fraction
  trans = exp( -tau )
  
!  write(*,*)"I've just calculated the transmission fraction..." 
  
! Apply covering fraction
  photar = (1.0-fcov) + fcov*trans
  
!  write(*,*)"I've just applied the covering fraction..."
  
  return
end subroutine mym
!=======================================================================





!-----------------------------------------------------------------------
subroutine extrap(ear,ne,z,filenm,Emax0,tau)
!
! Routine to extrapolate the optical depth above the highest energy
! in the grid. Assume an E^{-3} dependence above that.
!
! Emax0 = highest energy in the grid for z=0
!
  integer ne
  real ear(0:ne),z,Emax0,tau(ne)
  character (len=255) filenm
  real Emax,E
  integer i,ihi
  !Adjust for redshift
  Emax = Emax0 / (1.0+z)
  !Find the energy bin that this corresponds to
  i = 1
  do while( ear(i) .lt. Emax .and. i .le. ne )
     i = i + 1
  end do
  ihi = i - 1
  !Extrapolate
  do i = ihi+1,ne
     E      = 0.5 * ( ear(i) + ear(i-1) )
     tau(i) = tau(ihi) * ( E / ear(ihi) )**(-3)
  end do
  return
end subroutine extrap
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine extraplo(ear,ne,z,filenm,Emin0,tau)
!
! Routine to extrapolate the optical depth below the lowest energy
! in the grid. Assume an E^{-3} dependence below that.
!
! Emin0 = lowest energy in the grid for z=0
!
  integer ne
  real ear(0:ne),z,Emin0,tau(ne)
  character (len=255) filenm
  real Emin,E
  integer i,ilo
  !Adjust for redshift
  Emin = Emin0 / (1.0+z)
  !Find the energy bin that this corresponds to
  i = 1
  do while( ear(i-1) .lt. Emin .and. i .lt. ne )
     i = i + 1
  end do
  ilo = i
  !Extrapolate
  do i = 1,ilo-1
     E      = 0.5 * ( ear(i) + ear(i-1) )
     tau(i) = tau(ilo) * ( E / ear(ilo) )**(-3)
  end do
  return
end subroutine extraplo
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getNh(filenm,Nh)
  implicit none
  real Nh
  character (len=255) filenm
  integer status,unit,readwrite,blocksize,ecol
  logical exact,anynull
  character (len=50) comment

! Open etable file
  status = 0
  call ftgiou(unit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(unit,filenm,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open etable file'

! Shift to  extension "ENERGIES"
  status = 0
  call ftmnhd(unit,2,'PARAMETERS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension PARAMETERS'

! Get column density used
  call ftgkye(unit,'COLUMN',Nh,comment,status)
  if(status .ne. 0) stop 'Cannot determine COLUMN'

! Close file and free up unit
  call ftclos(unit,status)
  call ftfiou(unit,status)
  
  return
end subroutine getNh
!-----------------------------------------------------------------------
  



!-----------------------------------------------------------------------
subroutine getEmax(filenm,Emax)
  implicit none
  real Emax
  character (len=255) filenm
  integer status,unit,readwrite,blocksize,ecol,netab
  logical exact,anynull
  character (len=50) comment

! Open etable file
  status = 0
  call ftgiou(unit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(unit,filenm,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open etable file'

! Shift to  extension "ENERGIES"
  status = 0
  call ftmnhd(unit,2,'ENERGIES',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension ENERGIES'
      
! Get the column number for required inputs
  exact=.false.
  status = 0
  call ftgcno(unit,exact,'ENERG_HI',ecol,status)
  if(status.ne.0) stop 'Cannot determine column number'

! Get number of rows in the table
  call ftgkyj(unit,'NAXIS2',netab,comment,status)
  if(status .ne. 0) stop 'Cannot determine No of rows'

! Get final energy in table
  call ftgcve(unit,ecol,netab,1,1,-1.0,Emax,anynull,status)

! Close file and free up unit
  call ftclos(unit,status)
  call ftfiou(unit,status)
  
  return
end subroutine getEmax
!-----------------------------------------------------------------------
  


!-----------------------------------------------------------------------
subroutine getEmin(filenm,Emin)
  implicit none
  real Emin
  character (len=255) filenm
  integer status,unit,readwrite,blocksize,ecol,netab
  logical exact,anynull
  character (len=50) comment

! Open etable file
  status = 0
  call ftgiou(unit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(unit,filenm,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open etable file'

! Shift to  extension "ENERGIES"
  status = 0
  call ftmnhd(unit,2,'ENERGIES',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension ENERGIES'
      
! Get the column number for required inputs
  exact=.false.
  status = 0
  call ftgcno(unit,exact,'ENERG_HI',ecol,status)
  if(status.ne.0) stop 'Cannot determine column number'

! Get first energy in table
  call ftgcve(unit,ecol,1,1,1,-1.0,Emin,anynull,status)

! Close file and free up unit
  call ftclos(unit,status)
  call ftfiou(unit,status)
  
  return
end subroutine getEmin
!-----------------------------------------------------------------------
