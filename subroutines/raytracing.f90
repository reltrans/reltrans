module raytracing
!> This module contains the shim subroutines for raytracing and lensing factor 
!> calculations. All the subroutines only serve as interfaces to the actual 
!> implementations but allows for easy substitution without changing the rest of
!> the code.
    use rtconstants, only: pi
    implicit none

contains
    subroutine trace_disk_observer(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,&
                                    d)
    !> CALCULATION SUBROUTINE
    !> Traces rays in full GR for the camera defined by rn(nro), nro, nphi
    !> to convert alpha and beta to r and tau_do (don't care about phi)
    !> Used to be called GRtrace.
    !> Inputs:
    !>     nro: number of radial points in the disk.
    !>     nphi: number of azimuthal points in the disk.
    !>     rn: array of radial points in the disk.
    !>     mueff: effective cosine of the inclination angle of the disk.
    !>     mu0: cosine of the inclination angle of the observer.
    !>     spin: spin of the black hole.
    !>     rmin: minimum radius of the disk.
    !>     rout: maximum radius of the disk.
    !>     mudisk: cosine of the inclination angle of the disk.
    !>     d: distance from the black hole to the observer.
    !> Outputs:
    !>     pem1: array of p-coordinate at the disk for each ray.
    !>     taudo1: array of time coordinate at the disk for each ray.
    !>     re1: array of radial coordinate at the disk for each ray.
        use dyn_gr, only: pem1, re1, taudo1
        use kerrz, only: krz_TraceResult, trace_impact_parameters,             &
            KRZ_STATUS_NONE, kerr_metric
        implicit none
        integer, intent(in) :: nro,nphi
        double precision, intent(in) :: rn(nro),mueff,mu0,spin,rmin,rout
        double precision, intent(in) :: mudisk,d
        double precision :: phin, alpha, beta, cos0
        integer i,j
        type(krz_TraceResult) :: res
        cos0  = mu0
        taudo1   = 0.0
        re1      = 0.0
        !TODO: kerrz optimisation here! 
        do i = 1,nro
            do j = 1,NPHI
                phin  = (j-0.5) * 2.d0 * pi / dble(nphi)
                alpha = rn(i) * sin(phin)
                beta  = -rn(i) * cos(phin) * mueff

                res = trace_impact_parameters(mu0, alpha, beta, d)

                ! Do not include sub-isco contributions for now:
                if (res%status == KRZ_STATUS_NONE .and.            &
                    res%x_final%r > rmin .and. res%x_final%r < rout) then
                    ! TODO: This is for compatability with YNOGK. It is only
                    ! used in `strans` to later decide whether the photon hit
                    ! the disc or not. Any positive value can be used.
                    !
                    ! It is technically redundant, as a negative value could be
                    ! written into `re1` instead.
                    pem1(j,i) = 1.0d0
                    taudo1(j,i) = res%x_final%t - d
                    re1(j,i) = res%x_final%r
                else
                    pem1(j, i) = -1.0d0
                end if
              end do
        end do
        return
    end subroutine trace_disk_observer

    
    subroutine getdcos(a_spin,h,mudisk,n,nlp,rout,npts,r1,dcosdr,tc,cosd1,     &
                        cosdout)
    !> CALCULATION SUBROUTINE
    !> For n values of the emission angle, delta, the code calculates the r and t 
    !> coordinates for the geodesic for mu=mudisk; i.e. the crossing points of a 
    !> thin disk.
    !> Note that mudisk = (h/r) / sqrt( (h/r)**2 + 1 )
    !> INPUTS
    !>    a_spin       Dimensionless spin parameter
    !>    h            Height of on-axis, isotropically emitting source
    !>    mudisk       cos(theta) of disk surface (mu=0 for h/r=0)
    !>    n            Number of values of emission angle delta (see Fig 1 Dauser 
    !>                 et al 2013) calculated
    !>    rout         Disk outer radius
    !> OUTPUTS
    !>    npts         Number of points recorded in arrays (leq n, since some 
    !>                 trial values will not hit the disk)
    !>    r1(n)        Radius of disk crossing
    !>    dcosdr(n)    Corresponding d\cos\delta/dr
    !>    tc(n)        Corresponding time coordinate
    !>    cosd1(n)     Corresponding \cos\delta
    !>    cosdout      cosd at the outer disk radius
    !>
    !> TODO: replace this with a kerrz emissivity call.
        use kerrz, only: kerr_metric, krz_TraceResult, trace_lamppost,         &
            KRZ_STATUS_NONE
        implicit none
        double precision, intent(in )   :: a_spin, h(2), mudisk, rout
        integer         , intent(in )   :: n, nlp
        integer         , intent(inout) :: npts(nlp)
        double precision, intent(inout) :: r1(n,nlp)
        double precision, intent(out)   :: dcosdr(n,nlp), tc(n,nlp), cosd1(n,nlp), cosdout(nlp)
        integer  m,j,k,counter,nout(nlp)
        double precision rhorizon
        double precision deltamin,deltamax, deltas,r_min,r_max
        type(krz_TraceResult) :: res
        rhorizon = kerr_metric%horizon_radius
        ! Set minimum and maximum disk radii
        r_min = kerr_metric%isco
        r_max = 1d10

        ! Loop over each lamppost here:
        do m=1,nlp
            ! Calculate smallest delta worth considering
            deltamin = acos( h(m) / sqrt( h(m)**2 + rhorizon**2 ) )
            ! Consider arbitrarily large value of delta
            deltamax = pi
            ! Go through n different values of the angle delta_s
            counter = 0
            nout(m) = 1
            do j = 1,n
            ! Run through linear steps in the angle delta (see Fig 1; Dauser et
            ! al 2013)
                deltas   = deltamin + (j-1) * (deltamax-deltamin)/float(n-1)
                res = trace_lamppost(h(m), deltas)
                if (res%status == KRZ_STATUS_NONE                  &
                    .and. r_min <= res%x_final%r .and. r_max >= res%x_final%r  &
                ) then
                    counter = counter + 1
                    r1(counter,m) = res%x_final%r
                    cosd1(counter,m) = cos(deltas)
                    tc(counter,m) = res%x_final%t
                    if( rout .gt. r1(counter,m) ) nout(m) = counter
                end if
            end do 
            npts(m) = counter
        end do 
        
        !Calculate cosdout
        do m=1,nlp 
           if( nout(m) .eq. npts(m) )then
            !Extrapolate assuming Newtonian profile
                cosdout(m) = h(m)/sqrt(h(m)**2+rout**2)-h(m)/                  &
                            sqrt(h(m)**2+r1(npts(m),m)**2)+cosd1(npts(m),m)
            else
                !Inperpolate
                cosdout(m) = cosd1(nout(m)+1,m) - cosd1(nout(m),m)
                cosdout(m) = cosdout(m) * ( rout - r1( nout(m),m ) )
                cosdout(m) = cosdout(m) / ( r1(nout(m)+1,m) - r1(nout(m),m) )
                cosdout(m) = cosdout(m) + cosd1(nout(m),m)
            end if
        end do
        !Calculate d\delta/dr on the r-grid. Note that we need yet another loop 
        !over m because of how the counter npts is set up computationally, this 
        !costs no time whatsoever
        do m=1,nlp            
            npts(m) = npts(m) -1           
            do k = 1,npts(m)
                dcosdr(k,m) = abs((cosd1(k+1,m)-cosd1(k,m))/(r1(k+1,m)-r1(k,m)))
            end do
        end do
        !Discard the outer points as unreliable
        npts = npts - 7

        return
    end subroutine getdcos

    subroutine getlens(a_spin,h,muobs,lens,delt,cosdelta1)
    !> CALCULATION SUBROUTINE
    !> Routine to calculate the lensing factor l=d\cos\delta/d\cos(i)
    !> and the source to observer time lag.
    !> Both calculations need us to know the delta value for the geodesic
    !> that ends up at angle i at infinity.
    !> INPUTS
    !>     a_spin       Dimensionless spin parameter
    !>     h            Height of on-axis, isotropically emitting source
    !>     muobs        Cosine of inclination angle
    !
    !> OUTPUTS
    !>     lens         Lensing factor
    !>     delt         Source to observer time lag 
        use kerrz, only: trace_lensing, LamppostContinuum
        implicit none
        double precision, intent(in)    :: a_spin,h, muobs
        double precision, intent(inout) :: cosdelta1
        double precision, intent(out)   :: lens, delt
        double precision :: d
        double precision, parameter :: r_at_inf = 1.0d5
        type(LamppostContinuum) :: continuum
        continuum = trace_lensing(h, r_at_inf, muobs)
        d = continuum%alpha**2 + continuum%beta**2
        d = sqrt( r_at_inf**2 - d  )
        delt = continuum%time - d
        cosdelta1 = continuum%cos_delta
        lens = continuum%lensing_factor
        return
    end subroutine getlens

    
end module raytracing
