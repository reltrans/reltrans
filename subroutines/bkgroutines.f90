!-----------------------------------------------------------------------
subroutine readinbkg
! Reads in the background spectrum
! ***Must already know numchn***
! ***Must have already initialised bkgcounts and bkgrate***  
  use telematrix
  implicit none
  integer status,U1,readwrite,blocksize,i,colnum,felem
  integer nelem
  real nullval,Texp,bcorr, myenv_real
  logical anynull
  character (len=500) strenv
  character (len=200) comment, bkgenv

! Get the name of the background fits file
  bkgenv = 'BKG_SET'
  bkgname = strenv(bkgenv)
  if( trim(bkgname) .eq. 'none' )then
     write(*,*)"Enter name of the background file (with full path)"
     read(*,'(a)')bkgname
  endif

! Get background scaling factor
  bcorr = myenv_real("BACKSCL",0.0)
  if (bcorr .eq. 0.0) then
     write(*,*)"Enter BACKSCAL factor (enter 1 if you dont know what this is)"
     read(*,*)bcorr
  endif
  
! Open an unused unit
  status = 0
  call ftgiou(U1,status)
  
! Open the background fits file with read-only access
  readwrite = 0
  call ftopen(U1,bkgname,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open background file'

! Shift to SPECTRUM extension
  status = 0
  call FTMNHD(U1, 2, 'SPECTRUM', 0, status)
  if( status .ne. 0 ) stop 'cannot switch to SPECTRUM extension'

! Read in COUNTS column
  do i = 1,numchn
     colnum  = 2
     felem   = 1
     nelem   = 1
     nullval = -1.0
     anynull = .false.
     status = 0
     call FTGCVJ(U1,colnum,i,felem,nelem,nullval,bkgcounts(i),anynull,status)
     if( status .ne. 0 ) stop 'problem reading in counts column'
  end do

! Read in the exposure time
  call ftgkye(U1,'EXPOSURE',Texp,comment,status)

! Convert to count rate
  do i = 1,numchn
     bkgrate(i) = real( bkgcounts(i) ) / Texp
  end do

! Apply background scaling factor
  bkgrate = bkgrate * bcorr
  
! Close unit
  call ftclos(U1,status)
  call ftfiou(U1,status)

  return
end subroutine readinbkg
!-----------------------------------------------------------------------
  
