!*****************************************************************************************************
subroutine getdcos(a_spin,h,mudisk,n,nlp,rout)
    ! INPUTS
    ! a_spin       Dimensionless spin parameter
    ! h            Height of on-axis, isotropically emitting source
    ! mudisk       cos(theta) of disk surface (mu=0 for h/r=0)
    ! n            Number of values of emission angle delta (see Fig 1 Dauser et al 2013) calculated
    ! rout         Disk outer radius

    ! SETTING VALUES IN MODULE dyn_gr 
    ! npts         Number of points recorded in arrays (leq n, since some trial values will not hit the disk)
    ! dcosdr (n,nlp)    Corresponding d\cos\delta/dr        
    ! rlp    (n,nlp)    Radius of disk crossing
    ! tlp    (n,nlp)    Corresponding time coordinate
    ! cosd   (n,nlp)    Corresponding \cos\delta
    ! cosdout(  nlp)    cosd at the outer disk radius 

  ! For n values of the emission angle, delta, the code calculates the r and t coordinates
    ! for the geodesic for mu=mudisk; i.e. the crossing points of a thin disk.
    ! Note that mudisk = (h/r) / sqrt( (h/r)**2 + 1 )
  use blcoordinate
  use dyn_gr
  use isco
    implicit none
    integer         , intent(in) :: nlp
    double precision, intent(in) :: a_spin, h(nlp), mudisk, rout
    integer  m,j,n,k,counter,nout(nlp), t_r1, t_r2
    double precision sins,mus,lambda,q,scal,deltamin,deltamax
    double precision rhorizon,velocity(3),f1234(4),pp,pr,pt
    double precision deltas,r_min,r_max
    double precision rcros,mucros,phicros,tcros,sigmacros,pcros

    !      double precision cosphi,costheta,d1(n),sinphi,sintheta
    scal     = 1.d0   !Meaningless scaling factor
    mus      = 1.d0   !Position of source: mus=0 means on-axis
    sins     = 0.d0   !sin of same angle
    velocity = 0.0D0  !3-velocity of source

    ! open(31, file='delete_getdcos_photons.dat')
    !loop over h here
    do m=1,nlp
        !Calculate smallest delta worth considering
        deltamin = acos( h(m) / sqrt( h(m)**2 + rh**2 ) )
        !Consider arbitrarily large value of delta
        deltamax = pi
        !Set minimum and maximum disk radii
        ! r_min = disco( a_spin )
        r_min = rh
        r_max = 1d10
        !Go through n different values of the angle delta_s
        counter = 0
        nout(m) = 1
        do j = 1,n
        !Run through linear steps in the angle delta (see Fig 1; Dauser et al 2013)
            deltas   = deltamin + (j-1) * (deltamax-deltamin)/float(n-1)
            !Calculate 4-momentum in source rest frame tetrad
            pr = cos(deltas)           !cosdelta
            pp = sqrt( 1.d0 - pr**2 )  !sindelta
            pt= 0.d0
            !Convert to LNRF (locally non-rotating reference frame)
            call initialdirection(pr,pt,pp,sins,mus,a_spin,h(m),velocity,lambda,q,f1234)
            !Calculate value of p-coordinate at mu=0
            pcros = Pemdisk(f1234,lambda,q,sins,mus,a_spin,h(m),scal,mudisk,r_max,r_min)
            !From that, calculate r, phi and t at mu=0
            call YNOGK(pcros,f1234,lambda,q,sins,mus,a_spin,h(m),scal,rcros,mucros,phicros,tcros,sigmacros, t_r1, t_r2)
            ! write(31,*) deltas, pcros, rcros
            if( pcros .gt. 0.0 )then
                counter         = counter + 1
                rlp(counter,m)  = rcros
                cosd(counter,m) = pr    !cosdelta
                tlp(counter,m)  = tcros
                if( rout .gt. rlp(counter,m) ) nout(m) = counter
            end if
        end do 
        npts(m) = counter
    end do 
        
    ! close(31)
    !Calculate cosdout
    do m=1,nlp 
        if( nout(m) .eq. npts(m) )then
        !Extrapolate assuming Newtonian profile
            cosdout(m) = h(m)/sqrt(h(m)**2+rout**2)-h(m)/sqrt(h(m)**2+rlp(npts(m),m)**2)+cosd(npts(m),m)
        else
        !Inperpolate
            cosdout(m) = (cosd(nout(m)+1,m)-cosd(nout(m),m))*(rout-rlp(nout(m),m))/(rlp(nout(m)+1,m)-rlp(nout(m),m))
            cosdout(m) = cosdout(m) + cosd(nout(m),m)
        end if
    end do
    !Calculate d\delta/dr on the r-grid. Note that we need yet another loop over m because of how the counter npts is set up
    !computationally, this costs no time whatsoever
    do m=1,nlp            
        npts(m) = npts(m) -1           
        do k = 1,npts(m)
            dcosdr(k,m) = abs( ( cosd(k+1,m) - cosd(k,m) ) / ( rlp(k+1,m) - rlp(k,m) ) )
        end do
    end do
    !Discard the outer points as unreliable
    npts = npts - 7

    return
end subroutine getdcos
!*****************************************************************************************************
