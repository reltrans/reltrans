!-----------------------------------------------------------------------
      subroutine myxillDCp(ear, ne, param, ifl, photar)
      implicit none
      integer, intent(in)  :: ne,ifl
      real,    intent(in)  :: ear(0:ne), param(8)
      real,    intent(out) :: photar(ne)
      double precision     :: dxillpar(8), dear(0:ne), dphotar(ne)
      double precision     :: dphoter(ne)
      character (len=100)  :: char

      dxillpar = dble( param )
      dear = dble( ear )      

      call lmodxillverdensnthcompf(dear, ne, dxillpar, ifl, dphotar, dphoter, char)

      photar = sngl( dphotar )      
      return
    end subroutine myxillDCp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    subroutine myxill(ear, ne, param7, param8, ifl, Cp, photar)
!!! This routine calls the correct version of xillver based on the Cp !!!
!!!   Arg:
      !  ear: energy grid
      !  ne: number of grid points
      !  param7: array of xillver parameters for xillver,xillverD,xillverCp
      !  param8: array of xillver parameters for xillverDCp
      !  ifl : parameter to call xspec model
      !  Cp : chooses xillver model
      !      -1 xillver      1e15 density and powerlaw illumination  
      !       1 xillverD     high density and powerlaw illumination
      !      -2 xillverCp    1e15 density and nthcomp  illumination
      !       2 xillverDCp   high density and nthcomp  illumination
      ! photar: (output) xillver energy spectrum
!!!   Internal variables:
      !  

!!! Last change: Gullo - 2020 Jul
      
      implicit none
      integer, intent(in)  ::  ne, ifl, Cp
      real   , intent(in)  ::  ear(0:ne), param7(7), param8(8)
      real   , intent(out) :: photar(ne)

      double precision    :: dxillpar7(7), dxillpar8(8)
      double precision    :: dear(0:ne), dphotar(ne), dphoter(ne)
      character (len=100) char
      dear = dble( ear )      
      if( Cp .eq. -1 )then         !xillver
         dxillpar7 = dble( param7 )
         ! write(*,*) 'my xillver', dxillpar7
         call lmodxillverf(dear, ne, dxillpar7, ifl, dphotar, dphoter, char)
      else if ( Cp .eq. -2 )then   !xillverCp
         dxillpar7 = dble( param7 )
         ! write(*,*) 'my xillverCp', dxillpar7
         call lmodxillvernthcompf(dear, ne, dxillpar7, ifl, dphotar, dphoter, char)
      else if( Cp .eq. 1 )then     !xillverD
         dxillpar7 = dble( param7 )
         ! write(*,*) 'my xillverD', dxillpar7
         call lmodxillverdensf(dear, ne, dxillpar7, ifl, dphotar, dphoter, char)
      else if ( Cp .eq. 2 )then    !xillverDCp
         dxillpar8 = dble( param8 )
         ! write(*,*) 'my xillverDCp', dxillpar8
         call lmodxillverdensnthcompf(dear, ne, dxillpar8, ifl, dphotar, dphoter, char)
      else
         write(*,*) 'No xillver model available for this configuration'
         stop 
      end if
      photar = sngl( dphotar )      
      return
      end subroutine myxill
!-----------------------------------------------------------------------

! !-----------------------------------------------------------------------
!       subroutine myxill_hD(ear, ne, param, ifl, Cp, photar)
!       implicit none
!       integer, intent(in)  :: ne,ifl,Cp
!       real,    intent(in)  :: ear(0:ne), param(7)
!       real,    intent(out) :: photar(ne)
!       double precision     :: dxillpar(7), dear(0:ne), dphotar(ne)
!       double precision     :: dphoter(ne)
!       character (len=100)  :: char

!       dxillpar = dble( param )
!       dear = dble( ear )      

!       call lmodxillverdensf(dear, ne, dxillpar, ifl, dphotar, dphoter, char)

! !This if work only when we add xillver with nthComp for High density model 
! !       if( Cp .eq. 0 )then
! !         call lmodxillverdensf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
! !       else
! !          write(*,*) 'WARNING: there is no nthComp illumination for xillver if you want to use high density'
! !         call lmodxillverdensf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
! !       end if

!       photar = sngl( dphotar )      
!       return
!     end subroutine myxill_hD
! !-----------------------------------------------------------------------

! !-----------------------------------------------------------------------
!       subroutine myxill_T(ear,ne,param,kTbb,ifl,photar)
!       implicit none
!       integer ne,ifl,i
!       real ear(0:ne),param(7),photar(ne),kTbb,old_cont(ne)
!       double precision dxillpar(7),dear(0:ne),dphotar(ne)
!       double precision dphoter(ne),dold_cont(ne)
!       character (len=100) char
!       real comp_par(5),new_cont(ne),photer(ne),E,dE
!       dxillpar = dble( param )
!       dear = dble( ear )
!       !Call xillver with reflection fraction of -1
!       call lmodxillvercpf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
!       photar = sngl( dphotar )
!       !Now call again with reflection fraction = 0
!       dxillpar(7) = 0.d0
!       call lmodxillvercpf(dear,ne,dxillpar,ifl,dold_cont,dphoter,char)
!       old_cont = sngl( dold_cont )
!       !Now call nthcomp to make the new continuum
!       comp_par(1) = param(1) !real( Gamma )
!       comp_par(2) = param(3) !Ecut_obs
!       comp_par(3) = kTbb
!       comp_par(4) = 1.0            !diskbb type of seed spectrum
!       comp_par(5) = 0.0            !cosmological redshift is accounted for by the transfer function     
!       call nthcomp(ear,ne,comp_par,Ifl,new_cont,photer)
!       !Now convert reflection spectrum for the new continuum
!       do i = 1,ne
!          photar(i) = photar(i) * new_cont(i) / old_cont(i)
!          if( photar(i) .ne. photar(i) ) photar(i) = 0.0
!       end do
!       return
!       end subroutine myxill_T
! !-----------------------------------------------------------------------
