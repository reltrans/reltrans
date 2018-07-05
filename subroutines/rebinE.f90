!-----------------------------------------------------------------------
      subroutine rebinE(earx,px,nex,ear,p,ne)
      !General rebinning scheme, should be nice and robust - BUT IT FUCKING ISN'T
      !i,nex,earx,px = input
      !j,ne,ear,p    = output
      implicit none
      integer i,nex,j,ne
      real earx(0:nex),ear(0:ne),px(nex),p(ne),Ehigh,Elow,upper,lower
      real FRAC,Ej,Ei,pi,Ei2,pi2,grad,cons
      logical interp
      i = 1
      do j = 1,ne
        p(j) = 0.0
        Ehigh = ear(j)
        Elow  = ear(j-1)
        do while( earx(i) .le. Elow .and. i .lt. nex )
          i = i + 1
        end do
        interp = .true.
        do while(earx(i-1).lt.Ehigh.and.i.lt.nex)
          lower = MAX( earx(i-1) , Elow  )
          upper = MIN( earx(i)   , Ehigh )
          FRAC  = (upper-lower) / ( Ehigh - Elow )
          p(j)  = p(j) + px(i) * FRAC
          i = i + 1
          interp = .false.
        end do
        if( interp )then
          !Work out if it's ok to interpolate
          if( Elow  .ge. earx(nex) ) interp = .false.
          if( Ehigh .le. earx(0)   ) interp = .false.
          if( i     .ge. nex-1     ) interp = .false.
        end if
        if( interp )then
          !write(*,*)"need to interpolate!"
          !p(j) is interpolation between px(i) and px(i+1)
          !unless i=nex, in which case p(j) = px(i)
          Ej = 0.5 * ( Ehigh + Elow )        !Centre of newbin
          Ei = 0.5 * ( earx(i+1) + earx(i) ) !Centre of one oldbin
          pi = px(i+1)         !Value at bin centre
          if( Ei .eq. Ej )then
             p(j) = pi
          else
            if( Ei .gt. Ej )then
              Ei2 = 0.5 * (earx(i) + earx(i-1) )
              pi2 = px(i)
            else
              Ei2 = 0.5 * (earx(i+2) + earx(i+1) )
              pi2 = px(+2)
            end if
            grad = ( pi - pi2 ) / ( Ei - Ei2 )
            cons = 0.5 * ( pi + pi2 - grad*(Ei+Ei2) )
            p(j) = grad * Ej + cons
          end if
        end if
        if( i .gt. 2 ) i = i - 1
      end do
      RETURN
      END
!-----------------------------------------------------------------------


