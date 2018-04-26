!*****************************************************************************************************
      subroutine getdcos(a_spin,h,mudisk,n,rout,npts,r1,dcosdr,tc,cosd1,cosdout)
! INPUTS
! a_spin       Dimensionless spin parameter
! h            Height of on-axis, isotropically emitting source
! mudisk       cos(theta) of disk surface (mu=0 for h/r=0)
! n            Number of values of emission angle delta (see Fig 1 Dauser et al 2013) calculated
! rout         Disk outer radius
        
! OUTPUTS
! npts         Number of points recorded in arrays (leq n, since some trial values will not hit the disk)
! r1(n)        Radius of disk crossing
! dcosdr(n)    Corresponding d\cos\delta/dr
! tc(n)        Corresponding time coordinate
! cosd1(n)     Corresponding \cos\delta
! cosdout      cosd at the outer disk radius
!        
! For n values of the emission angle, delta, the code calculates the r and t coordinates
! for the geodesic for mu=mudisk; i.e. the crossing points of a thin disk.
! Note that mudisk = (h/r) / sqrt( (h/r)**2 + 1 )
      use blcoordinate
      implicit none
      double precision sins,mus,a_spin,h,lambda,q,scal,mudisk
      double precision rhorizon,velocity(3),f1234(4),pp,pr,pt
      double precision costheta,sintheta,cosphi,sinphi
      double precision deltamin,deltamax,rout,cosdout
      integer  j,n,k,counter,npts,nout
      double precision r1(n),d1(n)
      double precision dcosdr(n),tc(n)
      double precision deltas,cosd1(n),r_min,r_max,disco
      double precision rcros,mucros,phicros,tcros,sigmacros,pcros
      scal     = 1.d0   !Meaningless scaling factor
      mus      = 1.d0   !Position of source: mus=0 means on-axis
      sins     = 0.d0   !sin of same angle
      velocity = 0.0D0  !3-velocity of source
      rhorizon = one+sqrt(one-a_spin**2)
      !Calculate smallest delta worth considering
      deltamin = acos( h / sqrt( h**2 + rhorizon**2 ) )
      !Consider arbitrarily large value of delta
      deltamax = pi
      !Set minimum and maximum disk radii
      r_min = disco( a_spin )
      r_max = 1d10
      !Go through n different values of the angle delta_s
      counter = 0
      nout    = 1
      do j = 1,n
        !Run through linear steps in the angle delta (see Fig 1; Dauser et al 2013)
        deltas   = deltamin + (j-1) * (deltamax-deltamin)/float(n-1)
        !Calculate 4-momentum in source rest frame tetrad
        pr = cos(deltas)           !cosdelta
        pp = sqrt( 1.d0 - pr**2 )  !sindelta
        pt= 0.d0
        !Convert to LNRF (locally non-rotating reference frame)
        call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
        !Calculate value of p-coordinate at mu=0
        pcros = Pemdisk(f1234,lambda,q,sins,mus,a_spin,h,scal,mudisk,r_max,r_min)
        !From that, calculate r, phi and t at mu=0
        call YNOGK(pcros,f1234,lambda,q,sins,mus,a_spin,h,scal,rcros,mucros,phicros,tcros,sigmacros)
        if( pcros .gt. 0.0 )then
          !write(88,*)rcros,pr
          counter        = counter + 1
          r1(counter)    = rcros
          cosd1(counter) = pr    !cosdelta
          tc(counter)    = tcros
          if( rout .gt. r1(counter) ) nout = counter
        end if
      end do
      npts = counter
      
      !Calculate cosdout
      if( nout .eq. npts )then
        !Extrapolate assuming Newtonian profile
        cosdout = h/sqrt(h**2+rout**2)-h/sqrt(h**2+r1(npts)**2)+cosd1(npts)
      else
        !Inperpolate
        cosdout = (cosd1(nout+1)-cosd1(nout))*(rout-r1(nout))/(r1(nout+1)-r1(nout))
        cosdout = cosdout + cosd1(nout)
      end if
         
      !write(88,*)"no no"
      !write(88,*)rout,cosdout
      
      !Calculate d\delta/dr on the r-grid
      npts = npts - 1
      do k = 1,npts
        dcosdr(k) = abs( ( cosd1(k+1) - cosd1(k) ) / ( r1(k+1) - r1(k) ) )
      end do
      !Discard the outer points as unreliable
      npts = npts - 7
     !  !Write out
     !  do k = 1,npts
     !    if( r1(k) .gt. rms( a_spin ) )then
     !      write(92,*)r1(k),dcosdr(k),r1(k)*h/(h**2+r1(k)**2)**1.5
     !    end if
     !  end do
     ! write(92,*)"la x r (R\dg\u)"
     ! write(92,*)"la y dcos\gd/dr"
     ! write(92,*)"v .15 .15 .9 .6"
     ! write(92,*)"cs 1.5"
     ! write(92,*)"la file"
     ! write(92,*)"tim off"
     ! write(92,*)"lw 5"
     ! write(92,*)"skip on"
     ! write(92,*)"ls 2 on 2"
     ! write(92,*)"log x"
      return
      end subroutine getdcos
!*****************************************************************************************************
