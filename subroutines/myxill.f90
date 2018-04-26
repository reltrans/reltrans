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
