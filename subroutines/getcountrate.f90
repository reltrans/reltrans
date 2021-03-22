!-----------------------------------------------------------------------
function getcountrate(E1,E2,nex,earx,photarx)
  use telematrix
  implicit none
  integer :: nex
  real :: getcountrate,E1,E2,earx(0:nex),photarx(nex)
  real,    allocatable :: spec(:)
  integer :: I1,I2,i
  
!Read from response file
  if( needresp )then
     call initmatrix
  end if

!Allocate spectrum array
  if( .not. allocated(spec) ) allocate(spec(numchn))

!Convert energy range to channel range
  I1 = 1
  I2 = numchn
  do i = 0, numchn
     if( echn(i) .lt. E1 ) I1 = i
     if( echn(i) .le. E2 ) I2 = i
  end do
  I1 = I1 + 1
  if( I1 .gt. I2 ) I2 = I1
  
! Fold photarx spectrum around response matrix
  call fold(nex, earx, photarx, spec)
  
! Calculate count rate from spec
  getcountrate = 0.0
  do i = I1,I2
     getcountrate = getcountrate + spec(i)
  end do
  
  return
end function getcountrate
!-----------------------------------------------------------------------
