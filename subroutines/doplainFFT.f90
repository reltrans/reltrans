!------------------------------------------------------------------------
      subroutine doplainFFT(n,at,ReA,ImA)
      use rtconstants, only: wp
      implicit none
      integer n,j
      real(wp) at(n),ReA(0:n/2),ImA(0:n/2)
      real(wp) data(2*n)
      do j = 1,n
        data(2*j-1) = at(j)
        data(2*j)   = 0.0_wp
      end do
      call ourfour1(data,n,1)
      do j = 1, n/2
        ReA(j) = data(2*j+1)
        ImA(j) = data(2*j+2)
      end do
      ReA(0) = data(1)
      ImA(0) = 0.0_wp
      return
      end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
      subroutine doplaininvFFT(n,ReA,ImA,at)
      use rtconstants, only: wp
      implicit none
      integer n,j
      real(wp) at(n),ReA(0:n/2),ImA(0:n/2)
      real(wp) data(2*n)
! +ve frequencies
      do j = 1,n/2
        data(2*j+1) = ReA(j)
        data(2*j+2) = ImA(j)
      end do
! -ve frequencies
      do j = 1,n/2-1
        data(2*n-2*j+1) =  ReA(j)
        data(2*n-2*j+2) = -ImA(j)         
      end do
! DC component
      data(1) = ReA(0)
      data(2) = ImA(0)
      call ourfour1(data,n,-1)
      do j = 1,n
        at(j) = data(2*j-1) / real(n, wp)
      end do
      return
      end
!------------------------------------------------------------------------
