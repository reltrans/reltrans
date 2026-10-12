     
!-----------------------------------------------------------------------
      FUNCTION gasdev(idum)
!> Returns a normally distributed random number with zero mean and unit
!> variance, using `ran1` as the source of uniform deviates (polar
!> Box-Muller method, Numerical Recipes `gasdev`). Deviates are generated
!> in pairs; the second is saved and returned on the next call.
      use rtconstants, only: wp
      INTEGER idum
      real(wp) gasdev
!CU    USES ran1
      INTEGER iset
      real(wp) fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.0_wp*ran1(idum)-1.0_wp
        v2=2.0_wp*ran1(idum)-1.0_wp
        rsq=v1**2+v2**2
        if(rsq.ge.1.0_wp.or.rsq.eq.0.0_wp)goto 1
        fac=sqrt(-2.0_wp*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
    END FUNCTION gasdev
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      FUNCTION ran1(idum)
!> Returns a uniform random number in (0, 1), Numerical Recipes `ran1`.
!> This is the Park-Miller "minimal standard" generator
!> (idum = 16807 * idum mod (2**31 - 1)) with a 32-entry Bays-Durham
!> shuffle table which removes low-order serial correlations.
!>
!> Call with a negative `idum` to (re)initialise the sequence, and do not
!> alter `idum` between calls. The same seed should give the same sequence on
!> every compiler and platform.
      use rtconstants, only: sp
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      real(sp) ran1,AM,EPS,RNMX
!     NDIV = 1+(IM-1)/NTAB, written out to avoid an integer-division warning
      PARAMETER (IA=16807,IM=2147483647,AM=1.0_sp/real(IM, sp),IQ=127773,&
      IR=2836,NTAB=32,NDIV=67108864,EPS=1.2e-7_sp,RNMX=1.0_sp-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
!       Initialise: discard the first 8 values, then fill the shuffle table
        idum=max(-idum,1)
        do 11 j=NTAB+8,NTAB+1,-1
          k=idum/IQ
!         Schrage's method: computes IA*idum mod IM without integer overflow
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
11      continue
        do 12 j=NTAB,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          iv(j)=idum
12      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
   end function
!-----------------------------------------------------------------------
