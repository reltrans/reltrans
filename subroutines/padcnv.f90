!-----------------------------------------------------------------------
subroutine pad4FFT(ne,photar,padFT)
! Takes spectrum photar(1:ne), pads out with zeros to make it a length
! of 4*ne, and Fourier transforms to padFT(1:4*ne), which is a function
! of 1/E
  implicit none
  integer ne,i
  real photar(ne),padphot(4*ne)
  complex padFT(4*ne)
  
! Pad out the array
  padphot = 0.0
  do i = 1,ne
     padphot(i+2*ne) = photar(i)
  end do

! Call the energy Fourier transform code
  call E_FT(4*ne,padphot,padFT)

  return
end subroutine pad4FFT      
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine pad4invFFT(ne,padFT,conv)
! Takes padFT(1:4*ne), and zero-padded function of 1/E and inverse
! Fourier transforms to conv(ne), which is a non-zero padded spectrum
  implicit none
  integer ne,i
  real conv(ne),padconv(4*ne)
  complex padFT(4*ne)

! Inverse Fourier transform padded FT
  call E_invFT(4*ne,padFT,padconv)

! Populate output array
  do i = 1,ne
    conv(i) = padconv(i+5*ne/2)
  end do

  return
end subroutine pad4invFFT
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine E_invFT(nex,cc,conv)
! Takes the complex array cc(1:nex), which is a function of 1/E
! and inverse Fourier transforms to get back a real spectrum as a
! function of E, conv(1:nex)
  implicit none
  integer nex,i
  real conv(nex),cdata(2*nex)
  complex cc(nex)

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
end subroutine E_invFT
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine E_FT(nex,photarx,bc)
! Takes the real array photarx(1:nex), which is a spectrum as a
! function of photon energy E and Fourier transforms to bc(1:nex),
! which is complex and a function of 1/E.
! Uses FFTs, so nex must be a power of 2.
! Uses the inverse transform of four1.
  implicit none
  integer nex,i
  real photarx(nex)
  real bdata(2*nex)
  complex bc(nex)

! Move arrays into arrays for four1
  bdata = 0.0
  !-ve frequencies
  do i = 1,nex/2-1
    bdata(2*i+nex+1) = photarx(i)
  end do 
 !DC component
  bdata(1) = photarx(nex/2)
  !+ve frequencies
  do i = nex/2,nex
    bdata(2*i-nex+1) = photarx(i)
  end do

! #ifdef DEBUG
    ! do i = 1, nex
    !    write(61, *)  i , bdata(2*i - 1)
    ! end do
! #endif

! Now do the inverse FFT
  call ourfour1(bdata,nex,-1)

! Now put into complex arrays
  do i = 1,nex
    bc(i) = complex( bdata(2*i-1) , bdata(2*i) ) / sqrt(float(nex))
  end do

  ! do i = 1, nex
  !      write(64, *)  i , bc(i)
  ! end do

  return
end subroutine E_FT
!-----------------------------------------------------------------------
