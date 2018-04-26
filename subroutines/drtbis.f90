!-----------------------------------------------------------------------
      FUNCTION drtbis(func,x1,x2,xacc,par)
      implicit none
      INTEGER JMAX
      double precision drtbis,x1,x2,xacc,func,par(*)
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      double precision dx,f,fmid,xmid
      fmid=func(x2,par)
      f=func(x1,par)
      if(f*fmid.ge.0.) write(*,*) 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        drtbis=x1
        dx=x2-x1
      else
        drtbis=x2
        dx=x1-x2
      endif
      do j=1,JMAX
        dx=dx*.5
        xmid=drtbis+dx
        fmid=func(xmid,par)
        if(fmid.le.0.)drtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
      end do
      write(*,*) 'too many bisections in rtbis'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
!-----------------------------------------------------------------------
