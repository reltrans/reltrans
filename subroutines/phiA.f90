!-----------------------------------------------------------------------
      subroutine phaseA(nex,earx,contx,reconv,imconv,gso,zcos,Gamma,afac,lens,phiA)
! Self-consistently calculates the parameter phiA if phiset=1
      implicit none
      integer nex,phiset,myenv
      real earx(0:nex),contx(nex),reconv(nex),imconv(nex),afac,phiA,gso
      double precision zcos,Gamma,lens
      character (len=500) strenv,respname,arfname
!      character (len=3) exten
!      character (len=7) name
      logical needresp,needchans
      integer NENMAX,CHNMAX,NENERG,NUMCHN,Ilo,Ihi,I 
      parameter (NENMAX=5000,CHNMAX=5000)
      real EN(0:NENMAX),RESP(CHNMAX,NENMAX),ECHN(0:CHNMAX),E
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX),LCHAN(NENMAX,CHNMAX)
      real conti(nenmax),Rei(nenmax),Imi(nenmax),fI(chnmax),ReWI(chnmax)
      real ImWI(chnmax),Im,Re,Elo,Ehi
      logical arf
      data needresp/.true./
      data needchans/.true./
      common /need/ needresp,nenerg,En,resp,numchn,echn,ngrp,fchan,lchan,respname,arfname
      !save needresp,nenerg,En,resp,numchn,echn,ngrp,fchan,lchan,Ilo,Ihi,needchans
      save Ilo,Ihi,needchans,Elo,Ehi
      
      !Get environment variable phiset
      phiset  = myenv("PHI_SET",0)      !phiA is a parameter (0) or calculated (1)
      
      !If it is set to 1, then calculate phiA
      if( phiset .eq. 1 )then         
        !Read in instrument response matrix
        if( needresp )then
          !Get name of response file and arf file
          respname = strenv('RMF_SET')
          arfname  = strenv('ARF_SET')
          !If this is not set, ask for it
          if( trim(respname) .eq. 'none' )then
            write(*,*)"Enter name of the response file (with full path)"
            read(*,'(a)')respname
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
          !Now read in the files
          if( trim(respname) .ne. 'none' )then            
            call getrsp1(arf,NENMAX,CHNMAX,respname,arfname,NENERG,&
                 En,resp,NGRP,FCHAN,LCHAN,NUMCHN,ECHN)
          end if
          needresp = .false.
        end if
        !Now calculate phiA
        if( trim(respname) .ne. 'none' )then
          !Read in min and max channels
          if( needchans )then
            call fetchchans(numchn,Ilo,Ihi)
            needchans = .false.
          end if
          !Re-bin onto internal telescope grid
          call rebinE(earx,contx,nex,En,conti,nenerg)
          call rebinE(earx,reconv,nex,En,Rei,nenerg)
          call rebinE(earx,imconv,nex,En,Imi,nenerg)
          !Then fold around instrument response
          call fold(NENMAX,CHNMAX,NENERG,conti,RESP,NGRP,FCHAN,LCHAN,fI)
          call fold(NENMAX,CHNMAX,NENERG,Rei,RESP,NGRP,FCHAN,LCHAN,ReWI)
          call fold(NENMAX,CHNMAX,NENERG,Imi,RESP,NGRP,FCHAN,LCHAN,ImWI)
          !Sum up
          Im = 0.0
          Re = 0.0
          do I = Ilo,Ihi
            Im = Im - ImWI(I)
            Re = Re + lens * ( gso / (1.0+zcos) )**(2+Gamma) * fI(I)
            Re = Re + afac * ReWI(I)
          end do
          Im   = afac * Im
          phiA = atan2( Im , Re )
        else
          !Read in min and max energy range
          if( needchans )then
            write(*,*)"Enter lower energy in reference band"
            read(*,*)Elo
            write(*,*)"Enter upper energy in reference band"
            read(*,*)Ehi
            needchans = .false.
          end if
          !Sum up
          Im = 0.0
          Re = 0.0
          do i = 1,nex
            E = 0.5 * ( earx(i) + earx(i-1) )
            if( E .ge. Elo .and. E .le. Ehi )then
              Im = Im - imconv(i)
              Re = Re + lens * ( gso / (1.0+zcos) )**(2+Gamma) * contx(i)
              Re = Re + afac * reconv(i)
            end if
          end do
          Im   = afac * Im
          phiA = atan2( Im , Re )
        end if
      end if
      return
      end subroutine phaseA
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
      
      
!-----------------------------------------------------------------------
      subroutine fetchchans(numchn,Ilo,Ihi)
      implicit none
      integer numchn,Ilo,Ihi
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
      end subroutine fetchchans
!-----------------------------------------------------------------------
      

!-----------------------------------------------------------------------
      subroutine folder(nex,earx,ReGx,ImGx,ne,ear,ReG,ImG)
! Folds real and imaginary parts of the cross-spectrum around the
! telescope response
      implicit none
      integer nex,ne
      real earx(0:nex),ear(0:ne),ReGx(nex),ImGx(nex),ReG(ne),ImG(ne)
      character (len=500) strenv,respname,arfname
      logical needresp,arf
      integer NENMAX,CHNMAX,NENERG,NUMCHN
      parameter (NENMAX=5000,CHNMAX=5000)
      real EN(0:NENMAX),RESP(CHNMAX,NENMAX),ECHN(0:CHNMAX)      
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX),LCHAN(NENMAX,CHNMAX)
      real ReGi(nenmax),ImGi(nenmax),ReGtel(chnmax),ImGtel(chnmax)
!      real E,dE
!      character (len=3) exten
!      integer I,lresp
      data needresp/.true./
      common /need/ needresp,nenerg,En,resp,numchn,echn,ngrp,fchan,lchan,respname,arfname
!      save needresp,nenerg,En,resp,numchn,echn,ngrp,fchan,lchan

      !Read in response file
      if( needresp )then
        !Get name of response file and arf file
        respname = strenv('RMF_SET')
        arfname  = strenv('ARF_SET')        
        !If this is not set, ask for it
        if( trim(respname) .eq. 'none' )then
          write(*,*)"Enter name of the response file (with full path)"
          read(*,'(a)')respname
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
        !Now read in the files
        if( trim(respname) .ne. 'none' )then            
          call getrsp1(arf,NENMAX,CHNMAX,respname,arfname,NENERG,&
               En,resp,NGRP,FCHAN,LCHAN,NUMCHN,ECHN)
        end if
        needresp = .false.
      end if

      !Fold around that response
      if( respname .eq. 'none' )then
        stop "Error! ReIm parameter set to fold around response, but none loaded!"
      else
        !Re-bin onto internal telescope grid
        call rebinE(earx,ReGx,nex,En,ReGi,nenerg)
        call rebinE(earx,ImGx,nex,En,ImGi,nenerg)
        !Then fold around instrument response
        call fold(NENMAX,CHNMAX,NENERG,ReGi,RESP,NGRP,FCHAN,LCHAN,ReGtel)
        call fold(NENMAX,CHNMAX,NENERG,ImGi,RESP,NGRP,FCHAN,LCHAN,ImGtel)
        !Then re-bin back onto input XSPEC grid
        call rebinE(echn,ReGtel,numchn,ear,ReG,ne)
        call rebinE(echn,ImGtel,numchn,ear,ImG,ne)
      end if
     
      return
      end subroutine folder
!-----------------------------------------------------------------------


   
!-----------------------------------------------------------------------
      subroutine justresp(nenmax,chnmax,nenerg,En,resp,ngrp,fchan,lchan,numchn,echn)
! Fetches response matrix to use for calculating phiA
      implicit none
      integer nenmax,chnmax,nenerg,numchn
      real EN(0:NENMAX),RESP(CHNMAX,NENMAX),ECHN(0:CHNMAX)
      integer NGRP(NENMAX),FCHAN(NENMAX,CHNMAX),LCHAN(NENMAX,CHNMAX)
      character (len=500) respname,arfname
      character (len=3) exten
      integer lresp
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
      !Now get the response matrix
      call getrsp1(arf,NENMAX,CHNMAX,respname,arfname,NENERG,&
           En,resp,NGRP,FCHAN,LCHAN,NUMCHN,ECHN)
      return
      end subroutine justresp
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
