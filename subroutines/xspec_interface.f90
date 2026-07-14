module xspec_interface
    implicit none
    interface
        subroutine xsatbl(ear,ne,params,filename,ifl,photar,photer)            &
            bind(C, name="xsatbl")
            !> An interface to the xsatbl function in libXSFunctions.
            use iso_c_binding, only: c_float, c_int, c_char
            real(c_float), dimension(*), intent(in) :: ear, params
            real(c_float), dimension(*), intent(out) :: photar, photer
            character(kind=c_char), intent(in) :: filename(*)
            integer(c_int), value, intent(in) :: ne, ifl
        end subroutine xsatbl

        subroutine c_tbabs(earx, nex, params, Ifl, absorbx, photerx, str)      &
            bind(C, name = "C_tbabs")
            !> The interface around the XSPEC C_tbabs function.
            !> It will call the C_tbabs symbol in the libXSFunctions shared
            !> library.
            use iso_c_binding, only: c_double, c_int, c_char
            real(c_double), intent(in) :: earx(nex+1)
            real(c_double), intent(in) :: params(1)
            real(c_double), intent(inout) :: absorbx(nex), photerx(nex)
            integer(c_int), value, intent(in) :: nex, Ifl
            character(kind = c_char), intent(in) :: str(*)
        end subroutine c_tbabs

        subroutine donthcomp(earx, nex, params, Ifl, absorbx, photerx)
            !> The interface around the XSPEC donthcomp. Note that this is a
            !> Fortran function and so does not need to be bound to a C symbol.
            real, intent(in) :: earx(nex+1)
            real, intent(in) :: params(5)
            real, intent(inout) :: absorbx(nex), photerx(nex)
            integer, intent(in) :: nex, Ifl
        end subroutine donthcomp
    end interface
contains

    subroutine tbabs(earx, nex, nh, Ifl, absorbx, photerx)
        !> Call the tbabs function from the XSPEC model library.
        !>
        !> This wrapper temporarily performs a runtime cast on all of the arrays
        !> to double precision, as the C_tbabs function that we
        !> will eventually call expects double precision, whilst much of reltrans
        !> still uses `real`, which maps to single precision.
        !>
        !> Once the precision has been modified, this function can be simplified.
        real, intent(in) :: earx(0:nex), nh
        real, intent(inout) :: absorbx(nex), photerx(nex)
        integer, intent(in) :: nex, Ifl

        double precision :: d_earx(0:nex), d_absorbx(nex), d_photerx(nex),     &
            d_params(1)
        integer i

        do i = 0, nex
            d_earx(i) = earx(i)
        end do

        d_params(1) = nh
        call c_tbabs(d_earx, nex, d_params, Ifl, d_absorbx, d_photerx, "")

        do i = 1, nex
            absorbx(i) = d_absorbx(i)
            photerx(i) = d_photerx(i)
        end do
    end subroutine tbabs
end module xspec_interface
