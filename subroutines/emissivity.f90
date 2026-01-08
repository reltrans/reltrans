! vim: cc=80 wrap tw=80

! A module for the handling various emissivity functions and datas.
module emissivities
    implicit none
    abstract interface
        function external_emissivity(re, spin) bind(C) result(em)
            double precision, intent(in) :: re, spin
            double precision :: em
        end function external_emissivity
    end interface

    procedure(external_emissivity), pointer :: external_emissivity_ptr => null()

contains

    ! Used to assign the external emissivity pointer to register the callback
    ! function
    subroutine emissivity_ptr_set(fp) bind(C)
        procedure(external_emissivity), pointer, intent(in) :: fp
        external_emissivity_ptr => fp
    end subroutine emissivity_ptr_set

    ! Used to clear the global emissivity pointer and set it to null
    subroutine emissivity_ptr_free() bind(C)
        external_emissivity_ptr => null()
    end subroutine emissivity_ptr_free

    ! Calculate the emissivity at a particular radius on the disc `re` given the
    ! black hole `spin`, the power-law index `Gamma`, the cosine factor `cosfac
    ! = |dcos\delta/dr|`, the angular emissivity `ptf`, and the redshift along a
    ! geodesic from the source to the disc `gsd`
    function determine_emissivity(re, spin, Gamma, cosfac, ptf, gsd            &
        ) result(em)
        use constants
        double precision, intent(in) :: ptf, spin, cosfac, re, Gamma
        real, intent(in) :: gsd
        double precision :: em

        ! functions used
        double precision :: dareafac

        if (associated(external_emissivity_ptr)) then
            em = external_emissivity_ptr(re, spin)
        else
            em = gsd**Gamma * 2.d0 * pi * ptf
            em = em * cosfac / dareafac(re, spin)
        end if
    end function determine_emissivity
end module emissivities

