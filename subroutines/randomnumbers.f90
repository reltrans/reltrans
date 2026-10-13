module random_numbers
    !> Random number generation for the simulation models.
    !>
    !> These are the Numerical Recipes `ran1` and `gasdev` routines. These are
    !> used rather than the intrinsic `random_number` routine because the same
    !> seed gives the same sequence on every compiler and platform.
    !>
    !> Both routines keep their state between calls, so they are not thread
    !> safe.
    use rtconstants, only: wp
    implicit none
    private
    public :: ran1, gasdev

contains

    real(wp) function gasdev(idum)
        !> Returns a normally distributed random number with zero mean and
        !> unit variance, using `ran1` as the source of uniform deviates
        !> (polar Box-Muller method). Deviates are generated in pairs; the
        !> second is saved and returned on the next call.
        integer, intent(inout) :: idum
        logical, save :: have_spare = .false.
        real(wp), save :: spare
        real(wp) :: v1, v2, rsq, fac

        if (have_spare) then
            have_spare = .false.
            gasdev = spare
            return
        end if

        ! Pick a point uniformly in the unit circle (excluding the centre) by
        ! rejecting points from the enclosing square that fall outside it.
        do
            v1 = 2.0_wp * ran1(idum) - 1.0_wp
            v2 = 2.0_wp * ran1(idum) - 1.0_wp
            rsq = v1**2 + v2**2
            if (rsq < 1.0_wp .and. rsq > 0.0_wp) exit
        end do

        fac = sqrt(-2.0_wp * log(rsq) / rsq)
        spare = v1 * fac
        have_spare = .true.
        gasdev = v2 * fac
    end function gasdev

    real(wp) function ran1(idum)
        !> Returns a uniform random number in (0, 1).
        !>
        !> This is the Park-Miller "minimal standard" generator
        !> (idum = 16807 * idum mod (2**31 - 1)) with a 32-entry Bays-Durham
        !> shuffle table to remove low-order serial correlations.
        !>
        !> Call with a negative `idum` to (re)initialise the sequence, and do
        !> not alter `idum` between calls.
        integer, intent(inout) :: idum
        integer, parameter :: NTAB = 32
        ! Number of generator values per shuffle-table slot,
        ! 1 + (IM - 1) / NTAB, written out to avoid an integer-division warning.
        integer, parameter :: NDIV = 67108864
        real(wp), parameter :: AM = 1.0_wp / 2147483647.0_wp
        ! Largest value returned, so that the result is never exactly 1.
        real(wp), parameter :: RNMX = 1.0_wp - epsilon(1.0_wp)
        integer, save :: iv(NTAB) = 0, iy = 0
        integer :: j

        if (idum <= 0 .or. iy == 0) then
            ! Initialise: discard the first 8 values, then fill the shuffle
            ! table.
            idum = max(-idum, 1)
            do j = 1, 8
                call park_miller_step(idum)
            end do
            do j = NTAB, 1, -1
                call park_miller_step(idum)
                iv(j) = idum
            end do
            iy = iv(1)
        end if

        call park_miller_step(idum)
        ! Return the table entry chosen by the previous output, and replace
        ! it with the new generator value.
        j = 1 + iy / NDIV
        iy = iv(j)
        iv(j) = idum
        ran1 = min(AM * iy, RNMX)
    end function ran1

    subroutine park_miller_step(idum)
        !> Advances the Park-Miller generator: idum = IA * idum mod IM.
        !> Uses Schrage's method (IM = IA * IQ + IR) so that IA * idum never
        !> overflows a 32-bit integer.
        integer, intent(inout) :: idum
        integer, parameter :: IA = 16807, IM = 2147483647, IQ = 127773,       &
            IR = 2836
        integer :: k

        k = idum / IQ
        idum = IA * (idum - k * IQ) - IR * k
        if (idum < 0) idum = idum + IM
    end subroutine park_miller_step

end module random_numbers
