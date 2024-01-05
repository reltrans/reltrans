















  
!-----------------------------------------------------------------------
      subroutine getrsp1(arf,NENMAX,CHNMAX,respname,arfname,NENERG,&
        En,RD,NGRP,FCHAN,LCHAN,NUMCHN,ECHN)
! c Subroutine to read in a response file in standard fits format
! c INPUT:
! c ARF                Logical determining whether or not an arf is required
! c NENMAX             Maximum allowed value of NENERG
! c CHNMAX             Maximum allowed value of NUMCHN
! c respname           Name of response fits file
! c arfname            Name of arf file ('none' if just rsp)
! c OUTPUT
! c NENERG             Number of rows in the effective area table
! c En(0:NENERG)       Energy ranges in effective area table
! c RD(NUMCHN,NENERG)  Response matrix
! c NGRP(NENERG):      Number of groups with non-zero entries
! c FCHAN(NENERG,NGRP) First energy channel in a group with a non-zero entry
! c LCHAN(NENERG,NGRP) Last energy channel in a group with a non-zero entry
! c NUMCHN             Number of energy channels
! c ECHN(0:NUMCHN)     Energy range of each channel
      implicit none
      integer status,U1,readwrite,blocksize,hdutype,NENERG,NUMCHN
      integer NENMAX,colnum,J,CHNMAX,K   !,Lr
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX)
      integer LCHAN(NENMAX,CHNMAX),I,rows
      character (LEN=500) comment,RESPNAME,ARFNAME,EXNAME
!      character (LEN=3) EXTEN
      real En(0:NENMAX),nullval,ECHN(0:CHNMAX)
      real RD(CHNMAX,NENMAX),AREA(NENMAX)
      logical anynull,ARF
      !-----------------------------------------------------------------
      !Find a unit number not currently being used
      status = 0
      call ftgiou(U1,status)
      !-----------------------------------------------------------------
      !Open the response matrix fits file with read-only access
      readwrite = 0
      call ftopen(U1,RESPNAME,readwrite,blocksize,status)
      if( status .ne. 0 ) then
         write(*,*) 'response file name: ', RESPNAME
         stop 'cannot open response file'
      endif
      !-----------------------------------------------------------------
      !Shift to extension 1
      call ftmrhd(U1,1,hdutype,status)
      !Get the name of this extension
      call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
      !Read in whatever this extension this is
      if( EXNAME .eq. 'EBOUNDS' )then
        call readEBOUNDS(U1,CHNMAX,NUMCHN,ECHN)
      else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
        call readMATRIX(U1,NENMAX,CHNMAX,NENERG,En,RD,NGRP,FCHAN,LCHAN)
      else
        stop 'response has unrecognised extension name'
      end if
      !-----------------------------------------------------------------
      !Shift to extension 2
      call ftmrhd(U1,1,hdutype,status)
      !Get the name of this extension
      call ftgkys(U1,'EXTNAME',EXNAME,comment,status)
      !Read in whatever this extension this is
      if( EXNAME .eq. 'EBOUNDS' )then
        call readEBOUNDS(U1,CHNMAX,NUMCHN,ECHN)
      else if(EXNAME.eq.'MATRIX'.or.EXNAME.eq.'SPECRESP MATRIX')then
        call readMATRIX(U1,NENMAX,CHNMAX,NENERG,En,RD,NGRP,FCHAN,LCHAN)
      else
        stop 'response has unrecognised extension name'
      end if
      !Close unit
      call ftclos(U1,status)
      call ftfiou(U1,status)
      !-----------------------------------------------------------------
      !Read in the arf file if required
      if( ARF )then
        !Open file
        call ftgiou(U1,status)
        call ftopen(U1,ARFNAME,readwrite,blocksize,status)
        if( status .ne. 0 ) stop 'cannot open arf file'
        !Move to the SPECRESP extension
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
              RD(I,J) = RD(I,J) * AREA(J)
            end do
          end do
        end do
      end if
      !-----------------------------------------------------------------
      !Close unit
      call ftclos(U1,status)
      call ftfiou(U1,status)
    END subroutine getrsp1
!-----------------------------------------------------------------------
