!-----------------------------------------------------------------------
!      subroutine GRtrace(nmax,nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d,pem1,taudo1,re1)
      subroutine GRtrace(nmax,nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
! Traces rays in full GR for the camera defined by rn(nro), nro, nphi
! to convert alpha and beta to r and tau_do (don't care about phi)
        use dyn_gr
        use blcoordinate
      implicit none
      integer nmax,nro,nphi,i,j
      double precision rn(nro),mueff,mu0,spin,rmin,rout,mudisk,d
!      double precision pem1(nmax,nmax),taudo1(nmax,nmax),re1(nmax,nmax)
      double precision phin,alpha,beta,cos0,sin0,scal,velocity(3),f1234(4),lambda,q
      double precision pem,re,mucros,phie,taudo,sigmacros
      cos0  = mu0
      sin0  = sqrt(1.0-cos0**2)
      scal     = 1.d0
      velocity = 0.d0
      taudo1   = 0.0
      re1      = 0.0
      do i = 1,nro
        do j = 1,NPHI
          phin  = (j-0.5) * 2.d0 * pi / dble(nphi)
          alpha = rn(i) * sin(phin)
          beta  = -rn(i) * cos(phin) * mueff
          call lambdaq(-alpha,-beta,d,sin0,cos0,spin,scal,velocity,f1234,lambda,q)
          pem = Pemdisk(f1234,lambda,q,sin0,cos0,spin,d,scal,mudisk,rout,rmin)  !Can try rin instead of rmin to save an if statement
          pem1(j,i) = pem
          !pem > 1 means there is a solution
          !pem < 1 means there is no solution
          if( pem .gt. 0.0d0 )then
            call YNOGK(pem,f1234,lambda,q,sin0,cos0,spin,d,scal,re,mucros,phie,taudo,sigmacros)
            taudo1(j,i) = taudo - d
            re1(j,i)    = re
          end if
        end do
      end do
      return
      end subroutine GRtrace
!-----------------------------------------------------------------------
