!-----------------------------------------------------------------------
      subroutine phaseA(nex,earx,contx,reconv,imconv,gso,zcos,Gamma,afac,lens,phiA)
! Self-consistently calculates the parameter phiA if phiset=1
      implicit none
      integer nex,phiset,myenv
      real earx(0:nex),contx(nex),reconv(nex),imconv(nex),afac,phiA,gso
      double precision zcos,Gamma,lens
      logical needresp
      integer NENMAX,CHNMAX,NENERG,NUMCHN,Ilo,Ihi,I
      parameter (NENMAX=5000,CHNMAX=5000)
      real EN(0:NENMAX),RESP(CHNMAX,NENMAX),ECHN(0:CHNMAX)      
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX),LCHAN(NENMAX,CHNMAX)
      real conti(nenmax),Rei(nenmax),Imi(nenmax),fI(chnmax),ReWI(chnmax)
      real ImWI(chnmax),Im,Re
      data needresp/.true./
      save needresp,nenerg,En,resp,numchn,echn,ngrp,fchan,lchan,Ilo,Ihi
      
      !Get environment variable phiset
      phiset  = myenv("PHI_SET",0)      !phiA is a parameter (0) or calculated (1)

      !If it is set to 1, then calculate phiA
      if( phiset .eq. 1 )then
        !Read in instrument response matrix
        if( needresp )then
          call fetchresp(nenmax,chnmax,nenerg,En,resp,ngrp,fchan,lchan,numchn,echn,Ilo,Ihi)
          !Turn off needresp
          needresp = .false.
        end if
        !Re-bin onto internal telescope grid
        call rebinE(earx,contx,nex,En,conti,nenerg)
        call rebinE(earx,reconv,nex,En,Rei,nenerg)
        call rebinE(earx,imconv,nex,En,Imi,nenerg)
        !Then fold around instrument response
        call fold(NENMAX,CHNMAX,NENERG,conti,En,RESP,NGRP,FCHAN,LCHAN,fI)
        call fold(NENMAX,CHNMAX,NENERG,Rei,En,RESP,NGRP,FCHAN,LCHAN,ReWI)
        call fold(NENMAX,CHNMAX,NENERG,Imi,En,RESP,NGRP,FCHAN,LCHAN,ImWI)
        Im = 0.0
        Re = 0.0
        do I = Ilo,Ihi
          Im = Im - ImWI(I)
          Re = Re + lens * ( gso / (1.0+zcos) )**(2+Gamma) * fI(I)
          Re = Re + afac * ReWI(I)
        end do
        Im   = afac * Im
        phiA = atan2( Im , Re )
      end if
     
      return
      end subroutine phaseA
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine fetchresp(nenmax,chnmax,nenerg,En,resp,ngrp,fchan,lchan,numchn,echn,Ilo,Ihi)
! Fetches response matrix to use for calculating phiA
      implicit none
      integer nenmax,chnmax,nenerg,numchn
      real EN(0:NENMAX),RESP(CHNMAX,NENMAX),ECHN(0:CHNMAX)
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX),LCHAN(NENMAX,CHNMAX)
      character (len=500) respname,arfname
      character (len=3) exten
      integer lresp,Ilo,Ihi
      logical arf
      !First get the required files
      write(*,*)"Enter name of the response file (with full path)"
      read(*,'(a)')respname
      arf   = .false.
      lresp = len_trim( respname )
      exten = respname(lresp-2:lresp)
      if( exten .eq. 'rmf' .or. exten .eq. 'RMF' ) arf = .true.
      if( arf )then
        write(*,*)"Enter name of the anciliary (arf) response file (with full path)"
        read(*,'(a)')arfname
      end if
      !write(*,*)trim(respname)
      !write(*,*)trim(arfname)
      !Now get the response matrix
      call getrsp1(arf,NENMAX,CHNMAX,respname,arfname,NENERG,&
           En,resp,NGRP,FCHAN,LCHAN,NUMCHN,ECHN)
      !Now ask for energy channel range of reference band
      write(*,*)"Enter lowest channel number of the reference band"
      write(*,*)"(Convention: 1st channel is 1 not zero)"
      read(*,*)Ilo
      write(*,*)"Enter highest channel number of the reference band"
      read(*,*)Ihi
      if( Ihi .gt. numchn )then
        write(*,*)"Warning! Upper channel set too high! Changing to highest allowed."
        Ihi = numchn
      end if
      if( Ilo .lt. 1 )then
        write(*,*)"Warning! Lower channel set too low! Changing to channel 1"
        Ilo = 1
      end if
      return
      end subroutine fetchresp
!-----------------------------------------------------------------------
