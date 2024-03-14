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
      call ourfour(adata,nex,-1)
      call ourfour(bdata,nex,-1)

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
      call ourfour(cdata,nex,1)

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
    end subroutine dFTcnv
!-----------------------------------------------------------------------


! !-----------------------------------------------------------------------
!       subroutine FTcnv(nex,line,photarx,conv)
! ! Takes the arrays line(1-nex/2:nex/2) and photarx(1-nex/2:nex/2)
! ! and convolves them to get conv(1-nex/2:nex/2)
! ! Uses FFTs, so nex must be a power of 2.
!       implicit none
!       integer nex,i
!       real line(nex),photarx(nex),conv(nex)
!       real adata(2*nex),bdata(2*nex),cdata(2*nex)
!       complex ac(nex),bc(nex),cc(nex)

!       do i = 1, nex
!          write(10,*) i,  line(i)   
!       enddo
!       write(10, *) 'no no'
!       do i = 1, nex
!          write(11,*) i,  photarx(i)   
!       enddo
!       write(11, *) 'no no'

! ! Move arrays into arrays for four1
!       adata = 0.0
!       bdata = 0.0
!       !-ve frequencies
!       do i = 1,nex/2-1
!         adata(2*i+nex+1) = line(i)
!         bdata(2*i+nex+1) = photarx(i)
!       end do
!       !DC component
!       adata(1) = line(nex/2)
!       bdata(1) = photarx(nex/2)
!       !+ve frequencies
!       do i = nex/2,nex
!         adata(2*i-nex+1) = line(i)
!         bdata(2*i-nex+1) = photarx(i)
!       end do
      
! ! Now do the inverse FFT
!       call ourfour(adata,nex,-1)
!       call ourfour(bdata,nex,-1)
      
! ! Now put into complex arrays
!       do i = 1,nex
!         ac(i) = complex( adata(2*i-1) , adata(2*i) )
!         bc(i) = complex( bdata(2*i-1) , bdata(2*i) )
!       end do
!       do i = 1, nex
!          write(100,*) i,  real(ac(i))   
!       enddo
!          write(100,*) 'no no' 
!       do i = 1, nex
!          write(100,*) i,  aimag(ac(i))   
!       enddo
!          write(100,*) 'no no' 
!       do i = 1, nex
!          write(101,*) i,  real(bc(i))   
!       enddo
!          write(101,*) 'no no' 
!       do i = 1, nex
!          write(101,*) i,  aimag(bc(i))   
!       enddo
!          write(101,*) 'no no' 

! ! Multiply complex numbers together
!       cc = ac * bc / float(nex)

! ! Put back into four1 style arrays
!       do i = 1,nex
!         cdata(2*i-1) =  real( cc(i) )
!         cdata(2*i  ) = aimag( cc(i) )
!       end do
      
! ! Then transform back
!       call ourfour(cdata,nex,1)
      
! ! Move arrays back into original format
!       !-ve frequencies
!       do i = 1,nex/2-1
!         conv(i) = cdata(2*i+nex+1)
!       end do
!       !DC component
!       conv(nex/2) = cdata(1)
!       !+ve frequencies
!       do i = nex/2,nex
!         conv(i) = cdata(2*i-nex+1)
!       end do

      
!       do i = 1, nex
!          write(12,*) i,  conv(i)             
!       enddo
!       write(12, *) 'no no'
      

!       return
!     end subroutine FTcnv
! !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
      subroutine FTcnv(nex,line,photarx,conv)
! Takes the arrays line(1-nex/2:nex/2) and photarx(1-nex/2:nex/2)
! and convolves them to get conv(1-nex/2:nex/2)
! Uses FFTs, so nex must be a power of 2.
      implicit none
      integer nex,i
      real line(nex),photarx(nex),conv(nex)
      real adata(8*nex),bdata(8*nex),cdata(8*nex)
      complex ac(4*nex),bc(4*nex),cc(4*nex)

      do i = 1, nex
         write(10,*) i,  line(i)   
      enddo
      write(10, *) 'no no'
      do i = 1, nex
         write(11,*) i,  photarx(i)   
      enddo
      write(11, *) 'no no'

! Move arrays into arrays for four1
      adata = 0.0
      bdata = 0.0
      do i = 1,nex
        adata(2*i - 1) = line(i)
        bdata(2*i - 1) = photarx(i)
      end do
      
! Now do the inverse FFT
      call ourfour(adata,4*nex, 1)
      call ourfour(bdata,4*nex, 1)
      
! Now put into complex arrays
      do i = 1 , 4*nex
        ac(i) = complex( adata(2*i - 1) , adata(2*i) )
        bc(i) = complex( bdata(2*i - 1) , bdata(2*i) )
      end do

      do i = 1, nex
         write(100,*) i,  real(ac(i))   
      enddo
         write(100,*) 'no no' 
      do i = 1, nex
         write(100,*) i,  aimag(ac(i))   
      enddo
         write(100,*) 'no no' 
      do i = 1, nex
         write(101,*) i,  real(bc(i))   
      enddo
         write(101,*) 'no no' 
      do i = 1, nex
         write(101,*) i,  aimag(bc(i))   
      enddo
         write(101,*) 'no no' 

! Multiply complex numbers together
      cc = ac * bc / (4 * float(nex))

! Put back into four1 style arrays
      do i = 1, 4 * nex
        cdata(2*i-1) =  real( cc(i) )
        cdata(2*i  ) = aimag( cc(i) )
      end do
      
! Then transform back
      call ourfour(cdata,nex,-1)
      
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

      
      do i = 1, nex
         write(12,*) i,  conv(i)             
      enddo
      write(12, *) 'no no'
      

      return
    end subroutine FTcnv
!-----------------------------------------------------------------------






           ! do j = 1,nf
           !    do i = 1,nex
           !       reline(i)   = real(  transe(i,j,mubin,rbin) )
           !       imline(i)   = aimag( transe(i,j,mubin,rbin) )
           !       reline_a(i) = real(  transea(i,j,mubin,rbin) )
           !       imline_a(i) = aimag( transea(i,j,mubin,rbin) )
           !    end do              
           !    !Convolve with line profile
           !    !First FFTs
           !    call pad4FFT(nex,photarx,FTphotarx)
           !    call pad4FFT(nex,photarx_delta,FTphotarx_delta)
           !    call pad4FFT(nex,reline,FTreline)
           !    call pad4FFT(nex,imline,FTimline)
           !    call pad4FFT(nex,reline_a,FTreline_a)
           !    call pad4FFT(nex,imline_a,FTimline_a)
           !    call pad4FFT(nex,photarx_dlogxi,FTphotarx_dlogxi)
           !    !Then the multiplications and inverse FFTs
           !    FTreconv = FTreline * FTphotarx
           !    FTimconv = FTimline * FTphotarx
           !    call pad4invFFT(dyn,nex,FTreconv,reconvmu)
           !    call pad4invFFT(dyn,nex,FTimconv,imconvmu) 
           !    do i = 1,nex
           !       ReW0(i,j) = ReW0(i,j) + reconvmu(i)
           !       ImW0(i,j) = ImW0(i,j) + imconvmu(i)
           !    end do              
           !    FTreconv = FTreline_a * FTphotarx
           !    FTimconv = FTimline_a * FTphotarx
           !    call pad4invFFT(dyn,nex,FTreconv,reconvmu)
           !    call pad4invFFT(dyn,nex,FTimconv,imconvmu)
           !    do i = 1,nex
           !       ReW1(i,j) = ReW1(i,j) + reconvmu(i)
           !       ImW1(i,j) = ImW1(i,j) + imconvmu(i)
           !    end do
           !    FTreconv = FTreline * FTphotarx_delta
           !    FTimconv = FTimline * FTphotarx_delta
           !    call pad4invFFT(dyn,nex,FTreconv,reconvmu)
           !    call pad4invFFT(dyn,nex,FTimconv,imconvmu)
           !    do i = 1,nex
           !       ReW2(i,j) = ReW2(i,j) + reconvmu(i)
           !       ImW2(i,j) = ImW2(i,j) + imconvmu(i)
           !    end do
           !    FTreconv = FTreline * FTphotarx_dlogxi
           !    FTimconv = FTimline * FTphotarx_dlogxi
           !    call pad4invFFT(dyn,nex,FTreconv,reconvmu)
           !    call pad4invFFT(dyn,nex,FTimconv,imconvmu)
           !    do i = 1,nex
           !       ReW3(i,j) = ReW3(i,j) + reconvmu(i)
           !       ImW3(i,j) = ImW3(i,j) + imconvmu(i)
           !    end do
