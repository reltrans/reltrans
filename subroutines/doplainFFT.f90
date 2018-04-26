!------------------------------------------------------------------------
      subroutine doplainFFT(dt,n,at,ReA,ImA)
      implicit none
      integer n,j
      real at(n),ReA(0:n/2),ImA(0:n/2)
      real data(2*n),dt
      do j = 1,n
        data(2*j-1) = at(j)
        data(2*j)   = 0.0
      end do
      call four1(data,n,1)
      do j = 1, n/2
        ReA(j) = data(2*j+1)
        ImA(j) = data(2*j+2)
      end do
      ReA(0) = data(1)
      ImA(0) = 0.0
      return
      end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
      subroutine doplaininvFFT(dt,n,ReA,ImA,at)
      implicit none
      integer n,j
      real at(n),ReA(0:n/2),ImA(0:n/2)
      real data(2*n),dt
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
      call four1(data,n,-1)
      do j = 1,n
        at(j) = data(2*j-1) / float(n)
      end do
      return
      end
!------------------------------------------------------------------------
