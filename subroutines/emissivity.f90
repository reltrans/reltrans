! vim: cc=80 wrap tw=80

! A module for the handling various emissivity functions and datas.
module emissivities
    implicit none
    abstract interface
        function external_emissivity(re, spin) bind(C) result(em)
            double precision, intent(in) :: re, spin
            double precision :: em
        end function external_emissivity

        subroutine external_emissivity_time(re, phi, em, tau) bind(C)
            double precision, intent(in) :: re, phi
            double precision, intent(out) :: em, tau
        end subroutine external_emissivity_time

    end interface

    procedure(external_emissivity), pointer ::                                 &
        ext_emissivity => null()
    procedure(external_emissivity_time), pointer ::                            &
        ext_emissivity_time => null()

contains

    ! Used to assign the external emissivity pointer to register the callback
    ! function
    subroutine set_ext_emissivity_time(fp) bind(C)
        procedure(external_emissivity_time), pointer, intent(in) :: fp
        ext_emissivity_time => fp
    end subroutine set_ext_emissivity_time

    subroutine set_ext_emissivity(fp) bind(C)
        procedure(external_emissivity), pointer, intent(in) :: fp
        ext_emissivity => fp
    end subroutine set_ext_emissivity

    ! Used to clear the global emissivity pointer and set it to null
    subroutine free_ext_emissivity_time() bind(C)
        ext_emissivity_time => null()
    end subroutine free_ext_emissivity_time

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
    subroutine get_emissivity_time(re, phi, em, tau)
        double precision, intent(in) :: re, phi
        double precision, intent(out) :: em, tau
        if (associated(ext_emissivity_time)) then
            call ext_emissivity_time(re, phi, em, tau)
        else
            print *, "WARN: ext_emissivity_time not associated"
            tau = 0.d0
            em = 0.d0
        end if
    end subroutine get_emissivity_time

end module emissivities

