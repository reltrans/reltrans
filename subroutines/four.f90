!-----------------------------------------------------------------------
      SUBROUTINE ourfour1(data,nn,isign)
      use rtconstants, only: wp
      INTEGER isign,nn
      real(wp) data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      real(wp) tempi,tempr
      real(wp) theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959_wp/(isign*mmax)
        wpr=-2.0_wp*sin(0.5_wp*theta)**2
        wpi=sin(theta)
        wr=1.0_wp
        wi=0.0_wp
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=wr*data(j)-wi*data(j+1)
            tempi=wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!-----------------------------------------------------------------------
