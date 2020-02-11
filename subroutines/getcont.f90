!-----------------------------------------------------------------------
      subroutine getcontD(nex, earx, Gamma, Afe, logdens, logxi, Cp, contx, xillparD)
! Calculates continuum spectrum
      implicit none
      integer nex,ifl,Cp
      real earx(0:nex), Afe, logxi, contx(nex), xillparD(7)
      double precision Gamma
        
      !First continuum
      xillparD(1) = real( Gamma )
      xillparD(2) = Afe
      xillparD(3) = dens
      xillparD(4) = logxi
      xillparD(5) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillparD(6) = 30.0       !inclination angle (doesn't matter for continuum)
      xillparD(7) = 0.0       !reflection fraction of 0

      ifl        = 1

      call myxill_hD(earx, nex, xillparD, ifl, Cp, contx)
      return
    end subroutine getcontD
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine getcont(nex,earx,Gamma,Afe,Ecut_obs,logxi,Cp,contx,xillpar)
! Calculates continuum spectrum
      implicit none
      integer nex,ifl,Cp
      real earx(0:nex),Afe,Ecut_obs,logxi,contx(nex),xillpar(7)
      double precision Gamma
        
      !First continuum
      xillpar(1) = real( Gamma )
      xillpar(2) = Afe
      xillpar(3) = Ecut_obs
      xillpar(4) = logxi
      xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
      xillpar(7) = 0.0       !reflection fraction of 0
      ifl        = 1

      call myxill_hD(earx,nex,xillparD,ifl,Cp,contx)
      return
      end subroutine getcont
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine getcont_kTbb(nex,earx,Gamma,Afe,kTbb,Ecut_obs,logxi,contx,xillpar)
! Calculates continuum spectrum
      implicit none
      integer nex,ifl,i   !,Cp
      real earx(0:nex),Afe,Ecut_obs,logxi,contx(nex),xillpar(7),kTbb
      real photer(nex),comp_par(5),E,dE,norm6
      double precision Gamma
        
      !First continuum
      xillpar(1) = real( Gamma )
      xillpar(2) = Afe
      xillpar(3) = Ecut_obs
      xillpar(4) = logxi
      xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
      xillpar(7) = 0.0       !reflection fraction of 0
      ifl        = 1
      ! Cp = 1
      ! call myxill(earx,nex,xillpar,ifl,Cp,contx)
      comp_par(1) = real( Gamma )
      comp_par(2) = Ecut_obs
      comp_par(3) = kTbb
      comp_par(4) = 1.0            !diskbb type of seed spectrum
      comp_par(5) = 0.0            !cosmological redshift is accounted for by the transfer function     
      call nthcomp(earx,nex,comp_par,Ifl,contx,photer)

      ! do i = 1,nex
      !    E  = 0.5 * ( earx(i) + earx(i-1) )
      !    dE = earx(i) - earx(i-1)
      !    write(83,*)E,E**2*contx(i)/dE

      !    if( E .lt. 6.0 ) norm6 = E**2*contx(i)/dE
         
      ! end do

      ! write(*,*)"***norm6 (new)=",norm6
      
      return
      end subroutine getcont_kTbb
!-----------------------------------------------------------------------
