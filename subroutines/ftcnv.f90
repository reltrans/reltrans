!-----------------------------------------------------------------------
      subroutine dFTcnv(nex,line,photarx,conv)
! Takes the arrays line(1-nex/2:nex/2) and photarx(1-nex/2:nex/2)
! and convolves them to get conv(1-nex/2:nex/2)
! Uses FFTs, so nex must be a power of 2.
      implicit none
      integer nex,i
      double precision line(nex),photarx(nex),conv(nex)
      real adata(2*nex),bdata(2*nex),cdata(2*nex)
      complex ac(nex),bc(nex),cc(nex)

! Move arrays into arrays for four1
      adata = 0.0
      bdata = 0.0
      !-ve frequencies
      do i = 1,nex/2-1
        adata(2*i+nex+1) = real( line(i) )
        bdata(2*i+nex+1) = real( photarx(i) )
      end do
      !DC component
      adata(1) = real( line(nex/2) )
      bdata(1) = real( photarx(nex/2) )
      !+ve frequencies
      do i = nex/2,nex
        adata(2*i-nex+1) = real( line(i) )
        bdata(2*i-nex+1) = real( photarx(i) )
      end do
      
! Now do the inverse FFT
      call ourfour1(adata,nex,-1)
      call ourfour1(bdata,nex,-1)

! Now put into complex arrays
      do i = 1,nex
        ac(i) = complex( adata(2*i-1) , adata(2*i) )
        bc(i) = complex( bdata(2*i-1) , bdata(2*i) )
      end do

! Multiply complex numbers together
      cc = ac * bc / float(nex)

! Put back into four1 style arrays
      do i = 1,nex
        cdata(2*i-1) =  real( cc(i) )
        cdata(2*i  ) = aimag( cc(i) )
      end do
      
! Then transform back
      call ourfour1(cdata,nex,1)

! Move arrays back into original format
      !-ve frequencies
      do i = 1,nex/2-1
        conv(i) = dble( cdata(2*i+nex+1) )
      end do
      !DC component
      conv(nex/2) = dble( cdata(1) )
      !+ve frequencies
      do i = nex/2,nex
        conv(i) = dble( cdata(2*i-nex+1) )
      end do
      
      return
      end
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine FTcnv(nex,line,photarx,conv)
! Takes the arrays line(1-nex/2:nex/2) and photarx(1-nex/2:nex/2)
! and convolves them to get conv(1-nex/2:nex/2)
! Uses FFTs, so nex must be a power of 2.
      implicit none
      integer nex,i
      real line(nex),photarx(nex),conv(nex)
      real adata(2*nex),bdata(2*nex),cdata(2*nex)
      complex ac(nex),bc(nex),cc(nex)

! Move arrays into arrays for four1
      adata = 0.0
      bdata = 0.0
      !-ve frequencies
      do i = 1,nex/2-1
        adata(2*i+nex+1) = line(i)
        bdata(2*i+nex+1) = photarx(i)
      end do
      !DC component
      adata(1) = line(nex/2)
      bdata(1) = photarx(nex/2)
      !+ve frequencies
      do i = nex/2,nex
        adata(2*i-nex+1) = line(i)
        bdata(2*i-nex+1) = photarx(i)
      end do
      
! Now do the inverse FFT
      call ourfour1(adata,nex,-1)
      call ourfour1(bdata,nex,-1)
      
! Now put into complex arrays
      do i = 1,nex
        ac(i) = complex( adata(2*i-1) , adata(2*i) )
        bc(i) = complex( bdata(2*i-1) , bdata(2*i) )
      end do

! Multiply complex numbers together
      cc = ac * bc / float(nex)

! Put back into four1 style arrays
      do i = 1,nex
        cdata(2*i-1) =  real( cc(i) )
        cdata(2*i  ) = aimag( cc(i) )
      end do
      
! Then transform back
      call ourfour1(cdata,nex,1)
      
! Move arrays back into original format
      !-ve frequencies
      do i = 1,nex/2-1
        conv(i) = cdata(2*i+nex+1)
      end do
      !DC component
      conv(nex/2) = cdata(1)
      !+ve frequencies
      do i = nex/2,nex
        conv(i) = cdata(2*i-nex+1)
      end do
      
      return
      end
!-----------------------------------------------------------------------
