!-----------------------------------------------------------------------
subroutine rest_frame(ear,ne,Gamma,Afe,logne,Ecut,logxi,thetae,Cp,photar)
!
!  Cp : chooses reflection model
!      -1 xillver      1e15 density and powerlaw illumination  
!       1 xillverD     high density and powerlaw illumination 
!       2 xillverDCp   high density and nthcomp  illumination
!       0 reflionxDCp  reflionx high density and nthcomp  illumination
!
!       The first 2 have the same number of parameters xillpar(6), 
!       in the first one the Ecut is a parameter and the density is fixed to 10^15
!       in the second one the density is a parameter and Ecut is fixed to 300 keV 
!       The Cp = 2 has one more parameter both Ecut and density are parameters

!       Last change: Gullo - 2022 Oct

   implicit none
   integer, intent(in) :: ne, Cp
   real, intent(in)    :: ear(0:ne), Gamma, Afe, logne, Ecut, logxi, thetae
   real, intent(out)   :: photar(ne)
   real                :: xillpar(6), xillparDCp(7)
   integer :: i 
   
   if( Cp .ne. 0 )then
      !The model is a xillver model
   !   !Set density limits
   !   lognex = min(logne,22.0)
      !Fill parameter arrays
      xillpar(1) = Gamma     !Power law index
      xillpar(2) = Afe       !Iron abundance
      xillpar(3) = logxi     !ionization par
      xillpar(4) = Ecut      !Ecut or kTe
      if( Cp .eq. 1 )then
         xillpar(4) = logne !logne
      end if
      xillpar(5) = thetae    !emission angle
      xillpar(6) = 0.0       !redshift
      xillparDCp(1) = Gamma  !photon index
      xillparDCp(2) = Afe    !Afe
      xillparDCp(3) = logxi  !ionization par
      xillparDCp(4) = Ecut   !kTe
      xillparDCp(5) = logne !logne
      xillparDCp(6) = thetae !emission angle
      xillparDCp(7) = 0.0    !redshift
      write(*,*) 'logxi in rest frame ', logxi, xillparDCp(3)
      write(*,*) 'logne in rest frame ', logne, xillparDCp(5)
      call get_xillver(ear, ne, xillpar, xillparDCp, Cp, photar)
      ! photar = photar / 10**(logxi + logne - 15) !this factor is needed to match the normalisation with the first versions of reltrans
      ! write(*,*) 'xillver normalisation factor 10**(logxi + logne - 15)', 10**(logxi + logne - 15)
   else
      !The model is reflionx
      !Set density limits
   !   lognex = min(logne,22.0)
      call normreflionx(ear,ne,Gamma,Afe,logne,Ecut,logxi,thetae,photar)
   end if
   
   ! do i = 1, ne
   !    ! write(33,*) (ear(i) + ear(i-1))*0.5, &
   !    !      photar(i)/(ear(i) - ear(i-1)) * ((ear(i) + ear(i-1))*0.5)**2
   !    write(33,*) (ear(i) + ear(i-1))*0.5, photar(i)
   ! enddo
   
   return
 end subroutine rest_frame
!-----------------------------------------------------------------------


