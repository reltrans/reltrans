! vim: cc=80 wrap tw=80

! A module for the handling various emissivity functions and datas.
module emissivities
    implicit none
    abstract interface
        function external_emissivity(re, spin) bind(C) result(em)
            double precision, intent(in) :: re, spin
            double precision :: em
        end function external_emissivity

        subroutine external_emissivity_t_extrema(re, r_tmin, r_tmax) bind(C)
            double precision, intent(in) :: re
            double precision, intent(out) :: r_tmin, r_tmax
        end subroutine external_emissivity_t_extrema

        subroutine external_get_emissivity_slice(                              &
            re, emissivity_slice, n_emissivity) bind(C)
            double precision, intent(in) :: re
            double precision, pointer, intent(out) :: emissivity_slice(:)
            integer, intent(out) :: n_emissivity
        end subroutine external_get_emissivity_slice

        double precision function external_ring_emissivity(emissivity_slice,   &
                n_emissivity, p) bind(C) result(em)
            double precision, intent (in) :: emissivity_slice(:)
            integer, intent(in) :: n_emissivity
            real, intent(in) :: p
        end function external_ring_emissivity
    end interface

    procedure(external_emissivity), pointer ::                                 &
        ext_emissivity => null()
    procedure(external_emissivity_t_extrema), pointer ::                       &
        ext_emissivity_t_extrema => null()
    procedure(external_get_emissivity_slice), pointer ::                       &
        ext_get_emissivity_slice => null()
    procedure(external_ring_emissivity), pointer ::                            &
        ext_ring_emissivity => null()

contains

    ! Used to assign the external emissivity pointer to register the callback
    ! function
    subroutine set_ext_emissivity_t_extrema(fp) bind(C)
        procedure(external_emissivity_t_extrema), pointer, intent(in) :: fp
        ext_emissivity_t_extrema => fp
    end subroutine set_ext_emissivity_t_extrema

    subroutine set_ext_get_emissivity_slice(fp) bind(C)
        procedure(external_get_emissivity_slice), pointer, intent(in) :: fp
        ext_get_emissivity_slice => fp
    end subroutine set_ext_get_emissivity_slice

    subroutine set_ext_ring_emissivity(fp) bind(C)
        procedure(external_ring_emissivity), pointer, intent(in) :: fp
        ext_ring_emissivity => fp
    end subroutine set_ext_ring_emissivity

    subroutine set_ext_emissivity(fp) bind(C)
        procedure(external_emissivity), pointer, intent(in) :: fp
        ext_emissivity => fp
    end subroutine set_ext_emissivity

    ! Used to clear the global emissivity pointer and set it to null
    subroutine free_ext_emissivity_t_extrema() bind(C)
        ext_emissivity_t_extrema => null()
    end subroutine free_ext_emissivity_t_extrema

    subroutine free_ext_get_emissivity_slice() bind(C)
        ext_get_emissivity_slice => null()
    end subroutine free_ext_get_emissivity_slice

    subroutine free_ext_ring_emissivity() bind(C)
        ext_ring_emissivity => null()
    end subroutine free_ext_ring_emissivity

    subroutine free_ext_emissivity() bind(C)
        ext_emissivity => null()
    end subroutine free_ext_emissivity

    ! Calculate the emissivity at a particular radius on the disc `re` given the
    ! black hole `spin`, the power-law index `Gamma`, the cosine factor `cosfac
    ! = |dcos\delta/dr|`, the angular emissivity `ptf`, and the redshift along a
    ! geodesic from the source to the disc `gsd`
    double precision function determine_emissivity(re, spin, Gamma, cosfac,    &
        ptf, gsd) result(em)
        use constants
        double precision, intent(in) :: ptf, spin, cosfac, re, Gamma
        real, intent(in) :: gsd

        ! functions used
        double precision :: dareafac

        if (associated(ext_emissivity)) then
            em = ext_emissivity(re, spin)
        else
            em = gsd**Gamma * 2.d0 * pi * ptf
            em = em * cosfac / dareafac(re, spin)
        end if
    end function determine_emissivity

    ! Get the extrema source-to-disc time at a particular emissision radius on
    ! the disc `re`
    subroutine get_emissivity_t_extrema(re, r_tmin, r_tmax)
        double precision, intent(in) :: re
        double precision, intent(out) :: r_tmin, r_tmax
        if (associated(ext_emissivity_t_extrema)) then
            call ext_emissivity_t_extrema(re, r_tmin, r_tmax)
        else
            print *, "WARN: ext_emissivity_t_extrema not associated"
        end if
    end subroutine get_emissivity_t_extrema

    ! Get all emissivity values at a particular radius on the disc. Points the
    ! `emissivity_slice` pointer to the start of those elements, and sets
    ! `n_emissivity` to the number of items.
    ! TODO: the backing emissivity structure **must be** column major, so that
    ! each column is a particular re, else the pointer trick won't work
    subroutine get_emissivity_slice(re, emissivity_slice, n_emissivity)
        double precision, intent(in) :: re
        double precision, pointer, intent(out) :: emissivity_slice(:)
        integer, intent(out) :: n_emissivity
        if (associated(ext_get_emissivity_slice)) then
            call ext_get_emissivity_slice(re, emissivity_slice, n_emissivity)
        else
            print *, "WARN: ext_get_emissivity_slice not associated"
        end if
    end subroutine get_emissivity_slice

    ! Given an emissivity slice for a particular radius on the disc, interpolate
    ! the emissivity at `p`, where `p ∊ [0, 1]`
    double precision function ring_emissivity(emissivity_slice,                &
            n_emissivity, p) result(em)
        double precision, intent (in) :: emissivity_slice(:)
        integer, intent(in) :: n_emissivity
        real, intent(in) :: p
        if (associated(ext_ring_emissivity)) then
            em = ext_ring_emissivity(emissivity_slice, n_emissivity, p)
        else
            print *, "WARN: ext_ring_emissivity not associated"
        end if
    end function ring_emissivity

end module emissivities

