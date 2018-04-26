!-----------------------------------------------------------------------
      subroutine readMATRIX(U1,NENMAX,CHNMAX,NENERG,En,RD,NGRP,FCHAN,LCHAN)
! Routine to read MATRIX extension which has already been opened
! INPUT:
! U1                 Open fortran unit number
! NENMAX             Maximum allowed value of NENERG
! CHNMAX             Maximum allowed value of NUMCHN
! OUTPUT
! NENERG             Number of rows in the effective area table
! En(0:NENERG)       Energy ranges in effective area table
! RD(NUMCHN,NENERG)  Response matrix
! NGRP(NENERG):      Number of groups with non-zero entries
! FCHAN(NENERG,NGRP) First energy channel in a group with a non-zero entry
! LCHAN(NENERG,NGRP) Last energy channel in a group with a non-zero entry
      implicit none
      integer U1,NENMAX,CHNMAX,NENERG,NGRP(NENMAX),FCHAN(NENMAX,CHNMAX)
      integer LCHAN(NENMAX,CHNMAX),status,J,colnum,felem,nelem
      integer ARRAYJ(CHNMAX),K,I,I0,NCHAN(NENMAX,CHNMAX)
      real En(0:NENMAX),RD(CHNMAX,NENMAX),nullval,ARRAYE(CHNMAX)
      character (LEN=500) comment
      logical anynull
      status = 0
      !-----------------------------------------------------------------
      !Get number of rows in the table: NENERG
      call ftgkyj(U1,'NAXIS2',NENERG,comment,status)
      if (status .ne. 0) stop 'Cannot determine NENERG'
      if( NENERG .gt. NENMAX ) stop 'NENMAX is too small!!'
      !-----------------------------------------------------------------
      !Read from the binary tables
      RD = 0.0
      do J = 1,NENERG
        !---------------------------------------------------------------
        !Read in ENERG_LO
        colnum  = 1
        felem   = 1
        nelem   = 1
        nullval = -1.0
        anynull = .false.
        call ftgcve(U1,colnum,J,felem,nelem,nullval,En(J-1),anynull,status)
        if( status .ne. 0 ) stop 'problem reading ENERG_LO'
        !---------------------------------------------------------------
        !Read in ENERG_HI
        colnum  = 2
        call ftgcve(U1,colnum,J,felem,nelem,nullval,En(J),anynull,status)
        if( status .ne. 0 ) stop 'problem reading ENERG_HI'
        !---------------------------------------------------------------
        !Read in NGRP(J)
        colnum  = 3
        call ftgcvj(U1,colnum,J,felem,nelem,nullval,NGRP(J),anynull,status)
        if( status .ne. 0 ) stop 'problem reading NGRP'
        !---------------------------------------------------------------
        !Read in FCHAN(J,K)
        colnum = 4
        anynull = .false.
        call ftgcvj(U1,colnum,J,1,NGRP(J),nullval,ARRAYJ,anynull,status)
        do K = 1,NGRP(J)
          FCHAN(J,K) = ARRAYJ(K)
        end do
        if( status .ne. 0 ) stop 'problem reading FCHAN'
        !---------------------------------------------------------------
        !Read in NCHAN(J,K)
        colnum = 5
        call ftgcvj(U1,colnum,J,1,NGRP(J),nullval,ARRAYJ,anynull,status)
        do K = 1,NGRP(J)
          NCHAN(J,K) = ARRAYJ(K)
          LCHAN(J,K) = FCHAN(J,K)+NCHAN(J,K)
        end do
        if( status .ne. 0 ) stop 'problem reading NCHAN'
        !---------------------------------------------------------------
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
            RD(I,J) = ARRAYE(I0)
          end do
        end do
      end do
      RETURN
    END subroutine readMATRIX
!-----------------------------------------------------------------------

