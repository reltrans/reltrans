     
!-----------------------------------------------------------------------
      FUNCTION gasdev(idum)
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
      use rtconstants, only: wp
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      real(wp) ran1,AM,EPS,RNMX
!     NDIV = 1+(IM-1)/NTAB, written out to avoid an integer-division warning
      PARAMETER (IA=16807,IM=2147483647,AM=1.0_wp/real(IM, wp),IQ=127773,&
      IR=2836,NTAB=32,NDIV=67108864,EPS=1.2e-7_wp,RNMX=1.0_wp-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,NTAB+1,-1
          k=idum/IQ
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
