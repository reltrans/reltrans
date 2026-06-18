#ifndef _RELTRANS_IS_DEBUG
#define _RELTRANS_IS_DEBUG 0
#endif

module rtconstants
!> This module defines constants and enumerations used throughout the reltrans
!> codebase.
!> Physical constants should not be defined here, as a constants module already
!> exists in `amodules.f90`.
!> TODO: move `constants` in `amodules.f90` here.
    implicit none
    public

    !> The `re_im` parameter enumeration. Possible values are:
    !> These are the possible output modes that reltrans can compute.
    !> The folding refers to whether the reference (REF) band is folded through
    !> the instrument response, or the subject (SUBJ) band, or both (BOTH).
    integer, parameter ::                                                      &
        ! The time-averaged spectrum.
        MODE_TIME_AVG_SPECTRUM              =  0,                              &
        MODE_CROSS_SPEC_REAL                = -1,                              &
        MODE_CROSS_SPEC_REAL_REF_FOLDED     =  1,                              &
        MODE_CROSS_SPEC_IMAG                = -2,                              &
        MODE_CROSS_SPEC_IMAG_REF_FOLDED     =  2,                              &
        MODE_CROSS_SPEC_MODULUS             = -3,                              &
        MODE_CROSS_SPEC_MODULUS_REF_FOLDED  =  3,                              &
        MODE_CROSS_SPEC_LAG                 = -4,                              &
        MODE_CROSS_SPEC_LAG_REF_FOLDED      =  4,                              &
        MODE_CROSS_SPEC_MODULUS_BOTH_FOLDED =  5,                              &
        MODE_CROSS_SPEC_LAG_BOTH_FOLDED     =  6,                              &
        MODE_LAG_FREQ                       =  7,                              &
        MODE_CROSS_SPEC_LAG_TWO_RESPONSES   =  8

contains

    logical function is_mode(reim, mode) result(ret)
        !> Checks whether `reim` is of a particular mode regardless of whether
        !> it is folded or not.
        integer, intent(in) :: reim, mode
        integer local_reim
        local_reim = reim
        ! Folding should not matter in this function, hence convert this to what
        ! it physically represents.
        select case(abs(local_reim))
        case (MODE_CROSS_SPEC_MODULUS_BOTH_FOLDED)
            local_reim = MODE_CROSS_SPEC_MODULUS
        case (MODE_CROSS_SPEC_LAG_TWO_RESPONSES,                               &
              MODE_CROSS_SPEC_LAG_BOTH_FOLDED)
            local_reim = MODE_CROSS_SPEC_LAG
        end select
        ret = abs(local_reim) == abs(mode)
    end function is_mode

    integer function parse_reim(reim) result(ret)
        !> Convert the reim values to their `mode` representation.
        integer, intent(in) :: reim
        select case(reim)
        case (-4:8)
            ret = reim
        case default
            print *, "FATAL: Invalid reim parameter: ", reim
            stop 1
        end select
    end function parse_reim

    logical function is_folded(reim) result(bool)
        !> Returns true if the `reim` parameter is a folded result, in either
        !> the subject or reference band.
        integer, intent(in) :: reim
        if (reim > 0) then
            bool = .true.
        else
            bool = .false.
        end if
    end function is_folded

    logical function is_subj_folded(reim) result(bool)
        !> Returns true if in this mode the subject band is folded.
        integer, intent(in) :: reim
        select case(reim)
        case (MODE_CROSS_SPEC_LAG_BOTH_FOLDED,                                 &
              MODE_CROSS_SPEC_LAG_TWO_RESPONSES,                               &
              MODE_CROSS_SPEC_MODULUS_BOTH_FOLDED)
            bool = .true.
        case default
            bool = .false.
        end select
    end function is_subj_folded

    logical function is_both_folded(reim) result(bool)
        !> Returns true if in this mode both the subject and reference bands are
        !> folded.
        integer, intent(in) :: reim
        select case(abs(reim))
        case (MODE_CROSS_SPEC_LAG_BOTH_FOLDED,                                 &
              MODE_CROSS_SPEC_MODULUS_BOTH_FOLDED,                             &
              MODE_CROSS_SPEC_LAG_TWO_RESPONSES)
            bool = .true.
        case default
            bool = .false.
        end select
    end function is_both_folded

    logical function is_energy_dependent(reim) result(bool)
        !> Returns true if the `reim` parameter is a folded result.
        integer, intent(in) :: reim
        if (reim .ne. MODE_LAG_FREQ) then
            bool = .false.
        else
            bool = .true.
        end if
    end function is_energy_dependent

end module rtconstants
