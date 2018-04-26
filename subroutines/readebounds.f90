!-----------------------------------------------------------------------
      subroutine readEBOUNDS(U1,CHNMAX,NUMCHN,ECHN)
! Routine to read the EBOUNDS extension which has already been opened
! INPUT:
! U1             Open fortran unit number
! CHNMAX         Maximum allowed number of channels
! OUTPUT:
! NUMCHN         Number of energy channels
! ECHN(0:NUMCHN) Energy ranges of channels
      implicit none
      integer U1,CHNMAX,NUMCHN,status,I,colnum,felem,nelem
      real ECHN(0:CHNMAX),nullval
      character (LEN=500) comment
      logical anynull
      status = 0
      !-----------------------------------------------------------------
      !Get number of rows in the table: NUMCHN
      call ftgkyj(U1,'NAXIS2',NUMCHN,comment,status)
      if (status .ne. 0) stop 'Cannot determine NUMCHN'
      if( NUMCHN .gt. CHNMAX ) stop 'CHNMAX is too small!!'
      !-----------------------------------------------------------------
      !Read in the ranges of the energy channels
      do I = 1,NUMCHN
        !---------------------------------------------------------------
        !Read in E_MIN
        colnum  = 2
        felem   = 1
        nelem   = 1
        nullval = -1.0
        anynull = .false.
        call ftgcve(U1,colnum,I,felem,nelem,nullval,ECHN(I-1),anynull,status)
        !---------------------------------------------------------------
        !Read in E_MAX
        colnum  = 3
        call ftgcve(U1,colnum,I,felem,nelem,nullval,ECHN(I),anynull,status)
        if( status .ne. 0 ) stop 'problem reading in EBOUNDS'
! c        if( ECHN(I) .lt. ECHN(I-1) )then
! c          write(*,*)I,ECHN(I-1),ECHN(I)
! c          stop 'EBOUNDS are wrong'
! c        end if
      end do
      RETURN
      END
!-----------------------------------------------------------------------
