!-----------------------------------------------------------------------
      FUNCTION drtbis(func,x1,x2,xacc,par)
      use rtconstants, only: wp
      implicit none
      INTEGER JMAX
      real(wp) drtbis,x1,x2,xacc,func,par(*)
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      real(wp) dx,f,fmid,xmid
      fmid=func(x2,par)
      f=func(x1,par)
      if(f*fmid.ge.0.0_wp) write(*,*) 'root must be bracketed in rtbis'
      if(f.lt.0.0_wp)then
        drtbis=x1
        dx=x2-x1
      else
        drtbis=x2
        dx=x1-x2
      endif
      do j=1,JMAX
        dx=dx*0.5_wp
        xmid=drtbis+dx
        fmid=func(xmid,par)
        if(fmid.le.0.0_wp)drtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.0_wp) return
      end do
      write(*,*) 'too many bisections in rtbis'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
!-----------------------------------------------------------------------
