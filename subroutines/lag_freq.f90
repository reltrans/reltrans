subroutine lag_freq(nex, earx, Ea1, Ea2, Eb1, Eb2, contx, &
     ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3,ImW3,&
     absorbx, g, gslope,  DelAB, ABslope, nfx, fix, boost, z, gso, lens,&
     Gamma, ionvar, ReGraw, ImGraw)

  implicit none
  integer, intent(in) :: nex, nfx, ionvar
  integer, intent(in) :: Ea1, Ea2, Eb1, Eb2
  real   , intent(in) :: g, gslope, DelAB, ABslope, boost, z, gso, lens, Gamma
  real   , intent(in) :: earx(0:nex), contx(nex), absorbx(nex)
  real   , intent(in) :: fix(0:nfx)
  real   , intent(in) :: ReW0(nex,nfx), ImW0(nex,nfx), &
       ReW1(nex,nfx), ImW1(nex,nfx), ReW2(nex,nfx), ImW2(nex,nfx), &
       ReW3(nex,nfx), ImW3(nex,nfx)
       
                         
  real,    intent(out):: ReGraw(nfx), ImGraw(nfx)
  real                :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
  real                :: corr, sinD, cosD, E, fac, ReW0s, ImW0s,&
       ReWbs, ImWbs, ReW3s, ImW3s, gsoz
  real                :: f, DelAB_nu, g_nu
  integer             :: i, j



  gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
  corr = lens * gsoz**Gamma  !Correction factor for direct component
  !Now calculate the cross-spectrum (/complex covariance)
  do j = 1, nfx
     f = (fix(j) + fix(j - 1) ) * 0.5
     ! f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
     ! DelAB_nu = DelAB * (floHz / f)
     DelAB_nu = DelAB * ( (fix(1) + fix(0)) * 0.5  / f)**ABslope
     g_nu     = g * ( (fix(1) + fix(0)) * 0.5  / f)**gslope
     sinD = sin(DelAB_nu)
     cosD = cos(DelAB_nu)
     ReGrawEa = 0.0
     ImGrawEa = 0.0
     ReGrawEb = 0.0
     ImGrawEb = 0.0     
     do i = Ea1, Ea2
        E   = 0.5 * ( earx(i) + earx(i-1) )
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        ReGrawEa = ReGrawEa + cosD * ( fac * corr * contx(i) + ReWbs )
        ReGrawEa = ReGrawEa - sinD * ImWbs
        ReGrawEa = ReGrawEa * g_nu
        ReGrawEa = ReGrawEa + corr * contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGrawEa = ImGrawEa + sinD * ( fac * corr * contx(i) + ReWbs )
        ImGrawEa = ImGrawEa + cosD * ImWbs
        ImGrawEa = ImGrawEa * g_nu
        ImGrawEa = ImGrawEa + ImW0s + ImW3s
        !Account for absorption
        ReGrawEa = ReGrawEa * absorbx(i)
        ImGrawEa = ImGrawEa * absorbx(i)
     end do
     
     do i = Eb1, Eb2
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        ReGrawEb = ReGrawEb + cosD * ( fac * corr * contx(i) + ReWbs )
        ReGrawEb = ReGrawEb - sinD * ImWbs
        ReGrawEb = ReGrawEb * g_nu
        ReGrawEb = ReGrawEb + corr * contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGrawEb = ImGrawEb + sinD * ( fac * corr * contx(i) + ReWbs )
        ImGrawEb = ImGrawEb + cosD * ImWbs
        ImGrawEb = ImGrawEb * g_nu
        ImGrawEb = ImGrawEb + ImW0s + ImW3s
        !Account for absorption
        ReGrawEb = ReGrawEb * absorbx(i)
        ImGrawEb = ImGrawEb * absorbx(i)
     end do

     ! write(21, *) f, ReGrawEa
     ! write(22, *) f, ReGrawEb
     ! write(23, *) f, ImGrawEa
     ! write(24, *) f, ImGrawEb

!Now cross-spectrum between the two energy bands
     ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
     ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
     
  end do

end subroutine lag_freq



!-----------------------------------------------------------------------
subroutine lag_freq_resp(nex,earx,Emin,dloge,Ea1,Ea2,Eb1,Eb2,resp_matr_a,resp_matr_b,&
     contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
     absorbx, g, gslope,  DelAB, ABslope, nfx, fix, boost, z, gso, lens,&
     Gamma, ionvar, ReGraw, ImGraw)
! Calculates cross-spectrum between two bands (A and B) as a function of frequency.
! Ea1, Ea2, Eb1 and Eb2 are in ****keV**** and define the two band passes  
  implicit none
  integer, intent(in) :: nex, nfx, ionvar
  integer, intent(in) :: resp_matr_a,resp_matr_b
  real   , intent(in) :: Ea1, Ea2, Eb1, Eb2
  real   , intent(in) :: g, gslope, DelAB, ABslope, boost, z, gso, lens, Gamma
  real   , intent(in) :: earx(0:nex), Emin, dloge, contx(nex), absorbx(nex)
  real   , intent(in) :: fix(0:nfx)
  real   , intent(in) :: ReW0(nex,nfx), ImW0(nex,nfx), &
       ReW1(nex,nfx), ImW1(nex,nfx), ReW2(nex,nfx), ImW2(nex,nfx), &
       ReW3(nex,nfx), ImW3(nex,nfx)
       
                         
  real,    intent(out):: ReGraw(nfx), ImGraw(nfx)
  real                :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
  real                :: corr, sinD, cosD, E, fac, ReW0s, ImW0s,&
       ReWbs, ImWbs, ReW3s, ImW3s, gsoz
  real                :: f, DelAB_nu, g_nu
  integer             :: i, j
  real    :: thresh,Wx_a(nex),Wx_b(nex)
  integer :: ilo_a,ihi_a,ilo_b,ihi_b
  real    :: dReGrawEa,dImGrawEa,dReGrawEb,dImGrawEb

! Get weighting function for each band
  thresh = 1e-3
  call weightings(resp_matr_a,resp_matr_b,nex,earx,Emin,dloge,thresh,&
       Ea1,Ea2,Eb1,Eb2,ilo_a,ihi_a,Wx_a,ilo_b,ihi_b,Wx_b)
  
! Calculate FT of spectrum and do weighted sum over energies
  gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
  corr = lens * gsoz**Gamma  !Correction factor for direct component
  !Now calculate the cross-spectrum (/complex covariance)
  do j = 1, nfx
     f = (fix(j) + fix(j - 1) ) * 0.5
     ! f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
     ! DelAB_nu = DelAB * (floHz / f)
     DelAB_nu = DelAB * ( (fix(1) + fix(0)) * 0.5  / f)**ABslope
     g_nu     = g * ( (fix(1) + fix(0)) * 0.5  / f)**gslope
     sinD = sin(DelAB_nu)
     cosD = cos(DelAB_nu)
     ReGrawEa = 0.0
     ImGrawEa = 0.0
     ReGrawEb = 0.0
     ImGrawEb = 0.0
     
     do i = ilo_a, ihi_a
        E   = 0.5 * ( earx(i) + earx(i-1) )
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        dReGrawEa = cosD * ( fac * corr * contx(i) + ReWbs )
        dReGrawEa = dReGrawEa - sinD * ImWbs
        dReGrawEa = dReGrawEa * g_nu
        dReGrawEa = dReGrawEa + corr * contx(i) + ReW0s + ReW3s
        dReGrawEa = dReGrawEa * absorbx(i)
        dReGrawEa = dReGrawEa * Wx_a(i)
        ReGrawEa  = ReGrawEa + dReGrawEa
        !Imaginary part
        dImGrawEa = sinD * ( fac * corr * contx(i) + ReWbs )
        dImGrawEa = dImGrawEa + cosD * ImWbs
        dImGrawEa = dImGrawEa * g_nu
        dImGrawEa = dImGrawEa + ImW0s + ImW3s
        dImGrawEa = dImGrawEa * absorbx(i)
        dImGrawEa = dImGrawEa * Wx_a(i)
        ImGrawEa  = ImGrawEa + dImGrawEa
     end do
     
     do i = ilo_b, ihi_b
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        dReGrawEb = cosD * ( fac * corr * contx(i) + ReWbs )
        dReGrawEb = dReGrawEb - sinD * ImWbs
        dReGrawEb = dReGrawEb * g_nu
        dReGrawEb = dReGrawEb + corr * contx(i) + ReW0s + ReW3s
        dReGrawEb = dReGrawEb * absorbx(i)
        dReGrawEb = dReGrawEb * Wx_b(i)
        !Imaginary part
        dImGrawEb = sinD * ( fac * corr * contx(i) + ReWbs )
        dImGrawEb = dImGrawEb + cosD * ImWbs
        dImGrawEb = dImGrawEb * g_nu
        dImGrawEb = dImGrawEb + ImW0s + ImW3s
        dImGrawEb = dImGrawEb * absorbx(i)
        dImGrawEb = dImGrawEb * Wx_b(i)
     end do

!Now cross-spectrum between the two energy bands
     ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
     ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
     
  end do

end subroutine lag_freq_resp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine weightings(resp_matr_a,resp_matr_b,nex,earx,Emin,dloge,thresh,&
       Ea1,Ea2,Eb1,Eb2,ilo_a,ihi_a,Wx_a,ilo_b,ihi_b,Wx_b)
! Calculates weighting functions for bands A and B based on specified
! response matrices. For resp_matr < 0, no matrix is used (or rather a
! diagonal one).
  implicit none
  integer resp_matr_a,resp_matr_b,nex,ilo_a,ihi_a,ilo_b,ihi_b
  real earx(0:nex),Emin,dloge,thresh,Ea1,Ea2,Eb1,Eb2,Wx_a(nex),Wx_b(nex)
  if( resp_matr_a .gt. 0 )then
     call fill_matrix_a(resp_matr_a)
     call get_weight_a(nex,earx,Emin,dloge,thresh,Ea1,Ea2,ilo_a,ihi_a,Wx_a)
  else
     call get_weight_dummy(nex,Emin,dloge,Ea1,Ea2,ilo_a,ihi_a,Wx_a)     
  end if

  if( resp_matr_b .gt. 0 )then
     call fill_matrix_b(resp_matr_b)
     call get_weight_b(nex,earx,Emin,dloge,thresh,Eb1,Eb2,ilo_b,ihi_b,Wx_b)
  else
     call get_weight_dummy(nex,Emin,dloge,Eb1,Eb2,ilo_b,ihi_b,Wx_b)     
  end if
  return
end subroutine weightings
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine get_weight_dummy(nex,Emin,dloge,E1,E2,ilo,ihi,Wx)
! Calculates top hat function weighting for cases with no response
! function loaded
! IN
! nex,Emin,dloge                 ...internal logarithmic energy grid
! E1,E2                          ...band pass
!
! OUT
! ilo,ihi                        ...Wx(ilo:ihi) > 1, zero otherwise
! Wx(1:nex)                      ...weighting function defined on internal energy grid
  implicit none
  integer nex,ilo,ihi
  real Emin,dloge,E1,E2,Wx(nex)
  Wx  = 0.0
  ilo = ceiling( log10(E1/Emin) / dloge ) + 1
  ihi = floor(   log10(E2/Emin) / dloge )
  Wx(ilo:ihi) = 1.0
  return
end subroutine get_weight_dummy  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine get_weight_a(nex,earx,Emin,dloge,thresh,Ea1,Ea2,ilo_a,ihi_a,Wx_a)
  use telematrix_ab
! IN
! nex,earx(0:nex),Emin,dloge     ...internal logarithmic energy grid
! thresh                         ...threshold: if Wx < thresh*Wx_max, set to zero
! Ea1,Ea2                        ...band passes for band a
!
! OUT
! ilo_a,ihi_a                    ...Wx_a(ilo_a:ihi_a) > 0, zero otherwise
! Wx_a(1:nex)                    ...weighting function defined on internal energy grid
  implicit none
  integer nex,ilo_a,ihi_a
  real earx(0:nex),Emin,dloge,thresh,Ea1,Ea2,Wx_a(nex)
  integer Imin,Imax,J,I,Jmin_a,Jmax_a,nenerg_t
  real Wamax,Et_min,Et_max
  real, allocatable :: En_t(:),Wtel_t(:)
  
! Calculate a weighting factor (effective area of a given range of energy channels)
  call energy_limits(numchn_a,Echn_a,Ea1,Ea2,Imin,Imax)
  if( allocated(Wtel_a) ) deallocate(Wtel_a)
  allocate(Wtel_a(nenerg_a))
  Wamax = 0.0
  do J = 1,nenerg_a
     Wtel_a(J) = 0.0
     do I = Imin,Imax
        Wtel_a(J) = Wtel_a(J) + resp_a(I,J)
     end do
     Wamax = max( Wtel_a(J) , Wamax )
  end do

! Truncate the array
  J = 1
  do while( Wtel_a(J)/Wamax .lt. thresh .and. J .lt. nenerg_a )
     J = J + 1
  end do
  Jmin_a = J !First bin in which Wtel_a > thresh
  J = nenerg_a
  do while( Wtel_a(J)/Wamax .lt. thresh .and. J .gt. 1 )
     J = J - 1
  end do
  Jmax_a = J !Last bin in which Wtel_a > thresh
  Jmax_a = max( Jmax_a , Jmin_a )
  nenerg_t = Jmax_a - Jmin_a + 1
  if( allocated(Wtel_t) ) deallocate( Wtel_t )
  allocate( Wtel_t(nenerg_t) )
  if( allocated(En_t) ) deallocate( En_t )
  allocate( En_t(0:nenerg_t) )
  En_t(0) = En_a(Jmin_a-1)
  do J = Jmin_a,Jmax_a
     Wtel_t(J-Jmin_a+1) = Wtel_a(J)
     En_t(J-Jmin_a+1)   = En_a(J)
  end do
  
! Finally need to re-bin a weighting factor onto earx() grid
  call myinterp(nenerg_t,En_t,Wtel_t,nex,earx,Wx_a)

! Calculate bounds for the re-binned function
  Et_min = 0.5 * ( En_t(1) + En_t(0) )
  ilo_a = ceiling( log10(Et_min/Emin) / dloge ) + 1
  Et_max = 0.5 * ( En_t(nenerg_t) + En_t(nenerg_t-1) )
  ihi_a = floor(   log10(Et_max/Emin) / dloge )

  return
end subroutine get_weight_a
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine get_weight_b(nex,earx,Emin,dloge,thresh,Eb1,Eb2,ilo_b,ihi_b,Wx_b)
  use telematrix_ab
! IN
! nex,earx(0:nex),Emin,dloge     ...internal logarithmic energy grid
! thresh                         ...threshold: if Wx < thresh*Wx_max, set to zero
! Eb1,Eb2                        ...band passes for band b
!
! OUT
! ilo_b,ihi_b                    ...Wx_b(ilo_b:ihi_b) > 0, zero otherwise
! Wx_b(1:nex)                    ...weighting function defined on internal energy grid
  implicit none
  integer nex,ilo_b,ihi_b
  real earx(0:nex),Emin,dloge,thresh,Eb1,Eb2,Wx_b(nex)
  integer Imin,Imax,J,I,nenerg_t,Jmin_b,Jmax_b
  real Et_min,Et_max,Wbmax
  real, allocatable :: En_t(:),Wtel_t(:)

! Calculate b weighting factor (effective area of a given range of energy channels)
  call energy_limits(numchn_b,Echn_b,Eb1,Eb2,Imin,Imax)
  if( allocated(Wtel_b) ) deallocate(Wtel_b)
  allocate(Wtel_b(nenerg_b))
  Wbmax = 0.0
  do J = 1,nenerg_b
     Wtel_b(J) = 0.0
     do I = Imin,Imax
        Wtel_b(J) = Wtel_b(J) + resp_b(I,J)
     end do
     Wbmax = max( Wtel_b(J) , Wbmax )
  end do

! Truncate the array
  J = 1
  do while( Wtel_b(J)/Wbmax .lt. thresh .and. J .lt. nenerg_b )
     J = J + 1
  end do
  Jmin_b = J !First bin in which Wtel_b > thresh
  J = nenerg_b
  do while( Wtel_b(J)/Wbmax .lt. thresh .and. J .gt. 1 )
     J = J - 1
  end do
  Jmax_b = J !Last bin in which Wtel_a > thresh
  Jmax_b = max( Jmax_b , Jmin_b )
  nenerg_t = Jmax_b - Jmin_b + 1
  if( allocated(Wtel_t) ) deallocate( Wtel_t ) 
  allocate( Wtel_t(nenerg_t) )
  if( allocated(En_t) ) deallocate( En_t )
  allocate( En_t(0:nenerg_t) )
  En_t(0) = En_b(Jmin_b-1)
  do J = Jmin_b,Jmax_b
     Wtel_t(J-Jmin_b+1) = Wtel_b(J)
     En_t(J-Jmin_b+1)   = En_b(J)
  end do
  
! Finally need to re-bin b weighting factor onto earx() grid
  call myinterp(nenerg_t,En_t,Wtel_t,nex,earx,Wx_b)

! Calculate bounds for the re-binned function
  Et_min = 0.5 * ( En_t(1) + En_t(0) )
  ilo_b = ceiling( log10(Et_min/Emin) / dloge ) + 1
  Et_max = 0.5 * ( En_t(nenerg_t) + En_t(nenerg_t-1) )
  ihi_b = floor(   log10(Et_max/Emin) / dloge )

  return
end subroutine get_weight_b
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine myinterp(nfx,farx,Gfx,nf,far,Gf)
! Interpolates the function Gfx from the grid farx(0:nfx) to the
! function Gf on the grid far(0:nf)
  implicit none
  integer nfx,nf
  real farx(0:nfx),Gfx(nfx),far(0:nf),Gf(nf)
  integer ix,j
  real fx(nfx),f,fxhi,Gxhi,fxlo,Gxlo
! Define grid of central input frequencies
  do ix = 1,nfx
     fx(ix) = 0.5 * ( farx(ix) + farx(ix-1) )
  end do
! Run through grid of central output frequencies
  ix = 1
  do j = 1,nf
     !Find the input grid frequencies either side of the current
     !output grid frequency
     f = 0.5 * ( far(j) + far(j-1) )
     do while( fx(ix) .lt. f .and. ix .lt. nfx )
        ix = ix + 1
     end do
     ix = max( 2 , ix )
     fxhi = fx(ix)
     Gxhi = Gfx(ix)
     ix = ix - 1
     fxlo = fx(ix)
     Gxlo = Gfx(ix)
     !Interpolate
     Gf(j) = Gxlo + ( Gxhi - Gxlo ) * ( f - fxlo ) / ( fxhi - fxlo )
  end do
  return
end subroutine myinterp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine energy_limits(numchn,Echn,Emin,Emax,Imin,Imax)
  implicit none
  integer numchn,Imin,Imax,I
  real Echn(0:numchn),Emin,Emax
  I = 1
  do while( Echn(I-1) .lt. Emin .and. I .lt. numchn )
     I = I + 1
  end do
  Imin = I
  do while( Echn(I) .le. Emax .and. I .lt. numchn )
     I = I + 1
  end do
  Imax = I - 1
  return
end subroutine energy_limits
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine fill_matrix_a(resp_matr_a)
  use telematrix
  use telematrix2
  use telematrix3
  use telematrix_ab
  implicit none
  integer resp_matr_a
  
! Read in response files
  call get_response(resp_matr_a)
  
! Deallocate arrays (will only need in the reltrans code)
  if( allocated(En_a  ) ) deallocate(En_a  )
  if( allocated(Echn_a) ) deallocate(Echn_a)
  if( allocated(resp_a) ) deallocate(resp_a)

! Allocate arrays
  if( resp_matr_a .eq. 1 )then
     nenerg_a = nenerg
     numchn_a = numchn
  else if( resp_matr_a .eq. 2 )then
     nenerg_a = nenerg2
     numchn_a = numchn2
  else if( resp_matr_a .eq. 3 )then     
     nenerg_a = nenerg3
     numchn_a = numchn3
  else
     nenerg_a = nenerg
     numchn_a = numchn
  end if
  allocate( En_a(0:nenerg_a) )
  allocate( Echn_a(0:numchn_a) )
  allocate( resp_a(numchn_a,nenerg_a) )

! Fill arrays
  if( resp_matr_a .eq. 1 )then
     En_a   = En
     Echn_a = Echn
     resp_a = resp
  else if( resp_matr_a .eq. 2 )then
     En_a   = En2
     Echn_a = Echn2
     resp_a = resp2
  else if( resp_matr_a .eq. 3 )then     
     En_a   = En3
     Echn_a = Echn3
     resp_a = resp3
  else
     En_a   = En
     Echn_a = Echn
     resp_a = resp
  end if

  return
end subroutine fill_matrix_a
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine fill_matrix_b(resp_matr_b)
  use telematrix
  use telematrix2
  use telematrix3
  use telematrix_ab
  implicit none
  integer resp_matr_a,resp_matr_b
  
! Read in response files
  call get_response(resp_matr_b)
  
! Deallocate arrays (will only need in the reltrans code)
  if( allocated(En_b  ) ) deallocate(En_b  )
  if( allocated(Echn_b) ) deallocate(Echn_b)
  if( allocated(resp_b) ) deallocate(resp_b)

! Allocate arrays
  if( resp_matr_b .eq. 1 )then
     nenerg_b = nenerg
     numchn_b = numchn
  else if( resp_matr_b .eq. 2 )then
     nenerg_b = nenerg2
     numchn_b = numchn2
  else if( resp_matr_b .eq. 3 )then     
     nenerg_b = nenerg3
     numchn_b = numchn3
  else
     nenerg_b = nenerg
     numchn_b = numchn
  end if
  allocate( En_b(0:nenerg_b) )
  allocate( Echn_b(0:numchn_b) )
  allocate( resp_b(numchn_b,nenerg_b) )

! Fill arrays
  if( resp_matr_b .eq. 1 )then
     En_b   = En
     Echn_b = Echn
     resp_b = resp
  else if( resp_matr_b .eq. 2 )then
     En_b   = En2
     Echn_b = Echn2
     resp_b = resp2
  else if( resp_matr_b .eq. 3 )then     
     En_b   = En3
     Echn_b = Echn3
     resp_b = resp3
  else
     En_b   = En
     Echn_b = Echn
     resp_b = resp
  end if

  return
end subroutine fill_matrix_b
!-----------------------------------------------------------------------
