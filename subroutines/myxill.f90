!-----------------------------------------------------------------------
      subroutine myxill_hD(ear,ne,param,ifl,Cp,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(7),photar(ne)
      double precision dxillpar(7),dear(0:ne),dphotar(ne)
      double precision dphoter(ne)
      character (len=100) char
      dxillpar = dble( param )
      dear = dble( ear )      
      if( Cp .eq. 0 )then
        call lmodxillverdensf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
      else
         write(*,*) 'WARNING: there is no nthComp illumination for xillver if you want to use high density'
        call lmodxillverdensf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
!         call lmodxillvercpf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
      end if
      photar = sngl( dphotar )      
      return
    end subroutine myxill_hD
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine myxill(ear,ne,param,ifl,Cp,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(7),photar(ne)
      double precision dxillpar(7),dear(0:ne),dphotar(ne)
      double precision dphoter(ne)
      character (len=100) char
      dxillpar = dble( param )
      dear = dble( ear )      
      if( Cp .eq. 0 )then
        call lmodxillverf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
      else
        call lmodxillvercpf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
      end if
      photar = sngl( dphotar )      
      return
      end subroutine myxill
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine myxill_T(ear,ne,param,kTbb,ifl,photar)
      implicit none
      integer ne,ifl,i
      real ear(0:ne),param(7),photar(ne),kTbb,old_cont(ne)
      double precision dxillpar(7),dear(0:ne),dphotar(ne)
      double precision dphoter(ne),dold_cont(ne)
      character (len=100) char
      real comp_par(5),new_cont(ne),photer(ne),E,dE
      dxillpar = dble( param )
      dear = dble( ear )
      !Call xillver with reflection fraction of -1
      call lmodxillvercpf(dear,ne,dxillpar,ifl,dphotar,dphoter,char)
      photar = sngl( dphotar )
      !Now call again with reflection fraction = 0
      dxillpar(7) = 0.d0
      call lmodxillvercpf(dear,ne,dxillpar,ifl,dold_cont,dphoter,char)
      old_cont = sngl( dold_cont )
      !Now call nthcomp to make the new continuum
      comp_par(1) = param(1) !real( Gamma )
      comp_par(2) = param(3) !Ecut_obs
      comp_par(3) = kTbb
      comp_par(4) = 1.0            !diskbb type of seed spectrum
      comp_par(5) = 0.0            !cosmological redshift is accounted for by the transfer function     
      call nthcomp(ear,ne,comp_par,Ifl,new_cont,photer)
      !Now convert reflection spectrum for the new continuum
      do i = 1,ne
         photar(i) = photar(i) * new_cont(i) / old_cont(i)
         if( photar(i) .ne. photar(i) ) photar(i) = 0.0
      end do
      return
      end subroutine myxill_T
!-----------------------------------------------------------------------
