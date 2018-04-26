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
      call myxill(earx,nex,xillpar,ifl,Cp,contx)
      return
      end subroutine getcont
!-----------------------------------------------------------------------
