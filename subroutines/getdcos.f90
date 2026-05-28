!*****************************************************************************************************
subroutine getdcos(model_args, mudisk, cosdout)
    ! INPUTS
    ! model_args    Derived type containing a (spin), h (heights), nlp, rout
    ! mudisk        cos(theta) of disk surface (mu=0 for h/r=0)
    !
    ! OUTPUTS (via dyn_gr module globals)
    ! npts          Number of points recorded in arrays
    ! rlp(n,nlp)    Radius of disk crossing
    ! dcosdr(n,nlp) Corresponding d\cos\delta/dr
    ! tlp(n,nlp)    Corresponding time coordinate
    ! cosd(n,nlp)   Corresponding \cos\delta
    ! cosdout(nlp)  cosd at the outer disk radius (scalar output argument)
    !
    ! For ndelta values of the emission angle, delta, the code calculates
    ! the r and t coordinates for the geodesic for mu=mudisk; i.e. the
    ! crossing points of a thin disk.
    ! Note that mudisk = (h/r) / sqrt( (h/r)**2 + 1 )
    use common_types
    use dyn_gr
    use blcoordinate
    implicit none
    type(t_model_arguments), intent(in) :: model_args
    double precision, intent(in)        :: mudisk
    double precision, intent(out)       :: cosdout(model_args%nlp)

    double precision :: sins, mus, lambda, q, scal
    double precision :: rhorizon, velocity(3), f1234(4), pp, pr, pt
    double precision :: deltamin, deltamax
    integer  :: m, j, k, counter, nout(model_args%nlp)
    double precision :: deltas, r_min, r_max, disco
    double precision :: rcros, mucros, phicros, tcros, sigmacros, pcros

    scal     = 1.d0   !Meaningless scaling factor
    mus      = 1.d0   !Position of source: mus=0 means on-axis
    sins     = 0.d0   !sin of same angle
    velocity = 0.0D0  !3-velocity of source
    rhorizon = one + sqrt(one - model_args%a**2)

    !loop over h here
    do m = 1, model_args%nlp
        !Calculate smallest delta worth considering
        deltamin = acos(model_args%h(m)                                        &
            / sqrt(model_args%h(m)**2 + rhorizon**2))
        !Consider arbitrarily large value of delta
        deltamax = pi
        !Set minimum and maximum disk radii
        r_min = disco(model_args%a)
        r_max = 1d10
        !Go through ndelta different values of the angle delta_s
        counter = 0
        nout(m) = 1
        do j = 1, ndelta
        !Run through linear steps in the angle delta
            deltas   = deltamin + (j - 1)                                      &
                * (deltamax - deltamin) / float(ndelta - 1)
            !Calculate 4-momentum in source rest frame tetrad
            pr = cos(deltas)           !cosdelta
            pp = sqrt(1.d0 - pr**2)    !sindelta
            pt = 0.d0
            !Convert to LNRF (locally non-rotating reference frame)
            call initialdirection(pr, pt, pp, sins, mus, model_args%a,         &
                model_args%h(m), velocity, lambda, q, f1234)
            !Calculate value of p-coordinate at mu=0
            pcros = Pemdisk(f1234, lambda, q, sins, mus, model_args%a,         &
                model_args%h(m), scal, mudisk, r_max, r_min)
            !From that, calculate r, phi and t at mu=0
            call YNOGK(pcros, f1234, lambda, q, sins, mus, model_args%a,      &
                model_args%h(m), scal, rcros, mucros, phicros, tcros,          &
                sigmacros)
            if(pcros .gt. 0.0) then
                counter          = counter + 1
                rlp(counter, m)  = rcros
                cosd(counter, m) = pr    !cosdelta
                tlp(counter, m)  = tcros
                if(model_args%rout .gt. rlp(counter, m)) nout(m) = counter
            end if
        end do
        npts(m) = counter
    end do

    !Calculate cosdout
    do m = 1, model_args%nlp
        if(nout(m) .eq. npts(m)) then
        !Extrapolate assuming Newtonian profile
            cosdout(m) = model_args%h(m)                                       &
                / sqrt(model_args%h(m)**2 + model_args%rout**2)                &
                - model_args%h(m)                                              &
                / sqrt(model_args%h(m)**2 + rlp(npts(m), m)**2)                &
                + cosd(npts(m), m)
        else
        !Interpolate
            cosdout(m) = (cosd(nout(m) + 1, m) - cosd(nout(m), m))            &
                * (model_args%rout - rlp(nout(m), m))                          &
                / (rlp(nout(m) + 1, m) - rlp(nout(m), m))
            cosdout(m) = cosdout(m) + cosd(nout(m), m)
        end if
    end do
    !Calculate d\delta/dr on the r-grid
    do m = 1, model_args%nlp
        npts(m) = npts(m) - 1
        do k = 1, npts(m)
            dcosdr(k, m) = abs((cosd(k + 1, m) - cosd(k, m))                  &
                / (rlp(k + 1, m) - rlp(k, m)))
        end do
    end do
    !Discard the outer points as unreliable
    npts = npts - 7

    return
end subroutine getdcos
!*****************************************************************************************************

