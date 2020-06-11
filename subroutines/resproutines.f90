!-----------------------------------------------------------------------
subroutine cfoldandbin(nex,earx,ReGx,ImGx,ne,ear,ReG,ImG)
! Initmatrix must have alreadt been called
! Input:  {ReGx(nex),ImGx(nex)]: in units of photar; i.e. (dN/dE)*dE
! Output: {ReG(nex) ,ImG(nex) ]: in units of photar; i.e. (dN/dE)*dE
! G is folded around the instrument response and re-binned onto the
! input energy array ear(0:ne)
  use telematrix
  implicit none
  integer nex,ne,i,j,k
  real earx(0:nex),ReGx(nex),ImGx(nex),ear(0:ne),ReG(ne),ImG(ne)
  real ReGi(nenerg),ImGi(nenerg),E,dE,E2ReGx(nex),E2ImGx(nex)
  real ReGtel(numchn),ImGtel(numchn)
  
  !Convert to E^2*dN/dE for better accuracy
  do i = 1,nex
     E  = 0.5 * ( earx(i) + earx(i-1) )
     dE = earx(i) - earx(i-1)
     E2ReGx(i) = E**2 * ReGx(i) / dE
     E2ImGx(i) = E**2 * ImGx(i) / dE
  end do
     
  !Rebin input arrays onto interpal telescope energy grid
  call rebinE(earx,E2ReGx,nex,En,ReGi,nenerg)
  call rebinE(earx,E2ImGx,nex,En,ImGi,nenerg)

  !Convert back to (dN/dE)*dE
  do i = 1,nenerg
     E  = 0.5 * ( En(i) + En(i-1) )
     dE = En(i) - En(i-1)
     ReGi(i) = ReGi(i) / E**2 * dE
     ImGi(i) = ImGi(i) / E**2 * dE
  end do
     
  !Fold around response
  ReGtel = 0.0
  ImGtel = 0.0
  do J = 1,NENERG
     do K = 1,NGRP(J)
        do I = FCHAN(J,K)+1,LCHAN(J,K)
           dE = En(J) - En(J-1)
           ReGtel(I) = ReGtel(I) + ReGi(J) * RESP(I,J)
           ImGtel(I) = ImGtel(I) + ImGi(J) * RESP(I,J)
        end do
     end do
  end do

  !Rebin onto input energy grid
  call rebinE(echn,ReGtel,numchn,ear,ReG,ne)
  call rebinE(echn,ImGtel,numchn,ear,ImG,ne)
  
  return
end subroutine cfoldandbin
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine cfold(nex, earx, ReGx, ImGx, ReGtel, ImGtel)
! Initmatrix must have alreadt been called
! Input (ReGx,ImGx) is in terms of **PHOTAR**; i.e. (dN/dE)*dE
! RGtel, ImGtel is in count rate vs channel number
  use telematrix
  implicit none
  integer nex,i,j,k
  real earx(0:nex),ReGx(nex),ImGx(nex),ReGtel(numchn),ImGtel(numchn)
  real ReGi(nenerg),ImGi(nenerg),E,dE,E2ReGx(nex),E2ImGx(nex)
  
  !Convert to E^2*dN/dE for better accuracy
  do i = 1,nex
     E  = 0.5 * ( earx(i) + earx(i-1) )
     dE = earx(i) - earx(i-1)
     E2ReGx(i) = E**2 * ReGx(i) / dE
     E2ImGx(i) = E**2 * ImGx(i) / dE
  end do
  
  !Rebin input arrays onto interpal telescope energy grid
  call rebinE(earx,E2ReGx,nex,En,ReGi,nenerg)
  call rebinE(earx,E2ImGx,nex,En,ImGi,nenerg)
  
  !Convert back to (dN/dE)*dE
  do i = 1,nenerg
     E  = 0.5 * ( En(i) + En(i-1) )
     dE = En(i) - En(i-1)
     ReGi(i) = ReGi(i) / E**2 * dE
     ImGi(i) = ImGi(i) / E**2 * dE
  end do

  !Fold around response
  ReGtel = 0.0
  ImGtel = 0.0
  do J = 1, NENERG
     do K = 1,NGRP(J)
        do I = FCHAN(J,K) + 1, LCHAN(J,K)
           dE = En(J) - En(J-1) !I think this can be commented out
           ReGtel(I) = ReGtel(I) + ReGi(J) * RESP(I,J)
           ImGtel(I) = ImGtel(I) + ImGi(J) * RESP(I,J)
        end do
     end do
  end do

  
  return
end subroutine cfold
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine readinresp
! Reads in the response matrix
! ***Must already know numchn nd nenerg***
  use telematrix
  implicit none
  integer status,U1,readwrite,blocksize,hdutype,i,colnum,felem
  integer nelem,j,rows,k
  character (len=200) exname,comment
  real nullval,area(10000)
  logical anynull
  status = 0
  call ftgiou(U1,status)
  !Open the response matrix fits file with read-only access
  readwrite = 0
  call ftopen(U1,RESPNAME,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open response file'
  !Shift to extension 1
  call ftmrhd(U1,1,hdutype,status)
  !Get the name of this extension
  call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
  !Read in whatever this extension this is
  if( EXNAME .eq. 'EBOUNDS' )then     
     call energyextension(U1)
  else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
     call matrixextension(U1)
  end if
  !Shift to extension 2
  call ftmrhd(U1,1,hdutype,status)
  !Get the name of this extension
  call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
  if( EXNAME .eq. 'EBOUNDS' )then     
     call energyextension(U1)
  else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
     call matrixextension(U1)
  end if
  !Close unit
  call ftclos(U1,status)
  call ftfiou(U1,status)
  !-----------------------------------------------------------------
  !Read in the arf file if required
  if( ARF )then
     !Open file
     status = 0
     call ftgiou(U1,status)
     call ftopen(U1,ARFNAME,readwrite,blocksize,status)
     if( status .ne. 0 ) stop 'cannot open arf file'
     !Move to the SPECRESP extension
     status = 0
     call ftmrhd(U1,1,hdutype,status)
     if(status .ne. 0) stop 'Cannot move to the SPECRESP extension'
     !Check this has the same number of rows as the rmf file
     call ftgkyj(U1,'NAXIS2',rows,comment,status)
     if(status .ne. 0) stop 'Cannot determine NENERG from arf file'
     if( rows .ne. NENERG ) stop 'rmf and arf not compatible!'
     !Read in rows and re-normalise response matrix
     do J = 1,NENERG
        colnum = 3
        call ftgcve(U1,colnum,J,1,1,nullval,AREA(J),anynull,status)
        if( status .ne. 0 ) stop 'problem reading AREA'
        do K = 1,NGRP(J)
           do I = FCHAN(J,K)+1,LCHAN(J,K)
              resp(I,J) = resp(I,J) * AREA(J)
           end do
        end do
     end do
     !Close unit
     call ftclos(U1,status)
     call ftfiou(U1,status)
  end if
  return
end subroutine readinresp
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
subroutine energyextension(U1)
  use telematrix
  implicit none
  integer status,U1,i,colnum,felem,nelem
  character (len=200) exname,comment
  real nullval
  logical anynull
  !Read in En(numchn)
  do I = 1,NUMCHN
     status = 0
     !Read in E_MIN
     colnum  = 2
     felem   = 1
     nelem   = 1
     nullval = -1.0
     anynull = .false.
     call ftgcve(U1,colnum,I,felem,nelem,nullval,ECHN(I-1),anynull,status)
     !Read in E_MAX
     colnum  = 3
     call ftgcve(U1,colnum,I,felem,nelem,nullval,ECHN(I),anynull,status)
     if( status .ne. 0 ) stop 'problem reading in EBOUNDS'
  end do
  return
end subroutine energyextension  
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine matrixextension(U1)
  use telematrix
  implicit none
  integer status,U1,i,colnum,felem,nelem,j,k,i0,arrayj(5000)
  real nullval,arraye(5000)
  logical anynull
  !Read in resp(numchn,nenerg) and En
  resp = 0.0
  do J = 1,NENERG
     status = 0
     !Read in ENERG_LO
     colnum  = 1
     felem   = 1
     nelem   = 1
     nullval = -1.0
     anynull = .false.
     call ftgcve(U1,colnum,J,felem,nelem,nullval,En(J-1),anynull,status)
     if( status .ne. 0 ) stop 'problem reading ENERG_LO'
     !Read in ENERG_HI
     colnum  = 2
     call ftgcve(U1,colnum,J,felem,nelem,nullval,En(J),anynull,status)
     if( status .ne. 0 ) stop 'problem reading ENERG_HI'
     !Read in NGRP(J)
     colnum  = 3
     call ftgcvj(U1,colnum,J,felem,nelem,nullval,NGRP(J),anynull,status)
     if( status .ne. 0 ) stop 'problem reading NGRP'
     !Read in FCHAN(J,K)
     colnum = 4
     anynull = .false.
     call ftgcvj(U1,colnum,J,1,NGRP(J),nullval,ARRAYJ,anynull,status)
     do K = 1,NGRP(J)
        FCHAN(J,K) = ARRAYJ(K)
     end do
     if( status .ne. 0 ) stop 'problem reading FCHAN'
     !Read in NCHAN(J,K)
     colnum = 5
     call ftgcvj(U1,colnum,J,1,NGRP(J),nullval,ARRAYJ,anynull,status)
     do K = 1,NGRP(J)
        NCHAN(J,K) = ARRAYJ(K)
        LCHAN(J,K) = FCHAN(J,K)+NCHAN(J,K)
     end do
     if( status .ne. 0 ) stop 'problem reading NCHAN'
     !Read in MATRIX - first calculate number of elements per row
     colnum = 6
     nelem  = 0
     do K = 1,NGRP(J)
        nelem = nelem + NCHAN(J,K)
     end do
     !Now can read in and put in a sensible format
     call ftgcve(U1,colnum,J,1,nelem,nullval,ARRAYE,anynull,status)
     if( status .ne. 0 ) stop 'problem reading MATRIX'
     if( anynull ) write(*,*)"Null values in MATRIX"
     I0 = 0
     do K = 1,NGRP(J)
        do I = FCHAN(J,K)+1,LCHAN(J,K)
           I0 = I0 + 1
           resp(I,J) = ARRAYE(I0)
        end do
     end do
  end do
  return
end subroutine matrixextension      
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getdim(respname,nenerg,numchn)
  implicit none
  character (len=500) respname
  integer nenerg,numchn
  integer status,U1,readwrite,blocksize,hdutype
  character (len=500) exname,comment
  !Find a unit number not currently being used
  status = 0
  call ftgiou(U1,status)
  !Open the response matrix fits file with read-only access
  readwrite = 0
  call ftopen(U1,RESPNAME,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open response file'
  !Shift to extension 1
  call ftmrhd(U1,1,hdutype,status)
  !Get the name of this extension
  call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
  !Read in whatever this extension this is
  if( EXNAME .eq. 'EBOUNDS' )then
     !Get number of rows in the table: NUMCHN
     status = 0
     call ftgkyj(U1,'NAXIS2',NUMCHN,comment,status)
     if (status .ne. 0) stop 'Cannot determine NUMCHN'
  else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
     !Get number of rows in the table: NENERG
     status = 0
     call ftgkyj(U1,'NAXIS2',NENERG,comment,status)
     if (status .ne. 0) stop 'Cannot determine NENERG'
  end if
  !Shift to extension 2
  call ftmrhd(U1,1,hdutype,status)
  !Get the name of this extension
  call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
  if( EXNAME .eq. 'EBOUNDS' )then
     !Get number of rows in the table: NUMCHN
     status = 0
     call ftgkyj(U1,'NAXIS2',NUMCHN,comment,status)
     if (status .ne. 0) stop 'Cannot determine NUMCHN'
  else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
     !Get number of rows in the table: NENERG
     status = 0
     call ftgkyj(U1,'NAXIS2',NENERG,comment,status)
     if (status .ne. 0) stop 'Cannot determine NENERG'
  end if
  !Close unit
  call ftclos(U1,status)
  call ftfiou(U1,status)
  return
end subroutine getdim
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine initmatrix
  use telematrix
  implicit none
  character (len=500) strenv

!Get name of response file and arf file
  respname = strenv('RMF_SET')
  arfname  = strenv('ARF_SET')        
!If this is not set, ask for it
  if( trim(respname) .eq. 'none' )then
     write(*,*)"Enter name of the response file (with full path)"
     read(*,'(a)') respname
  end if
!Check if I need the arf
  call arfcheck(respname,arf)
!Look for arf
  if( arf )then
     !If not defined, ask for it
     if( trim(arfname) .eq. 'none' )then
        write(*,*)"Enter name of the anciliary (arf) response file (with full path)"
        read(*,'(a)')arfname
     end if
  end if

!Get the dimensions of the arrays in the matrix
  call getdim(respname,nenerg,numchn)
  
!Allocate the arrays
  allocate( En(0:nenerg) )
  allocate( Echn(0:numchn) )
  allocate( resp(numchn,nenerg) )
  allocate( Ngrp(nenerg) )
  allocate( fchan(nenerg,numchn) )
  allocate( lchan(nenerg,numchn) )
  allocate( nchan(nenerg,numchn) )
  
! Read matrix to fill the arrays
  call readinresp

  needresp = .false.
  needchans = .true.
  return
end subroutine initmatrix
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      subroutine arfcheck(respname,arf)
      implicit none
      character (len=500) respname
      logical arf
      integer lresp
      character (len=3) exten
      arf   = .false.
      lresp = len_trim( respname )
      exten = respname(lresp-2:lresp)
      if( exten .eq. 'rmf' .or. exten .eq. 'RMF' ) arf = .true.
      return
      end subroutine arfcheck
!-----------------------------------------------------------------------
