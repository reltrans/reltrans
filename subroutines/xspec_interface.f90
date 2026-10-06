module xspec_interface
    use rtconstants, only: wp
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
            integer(c_int), value, intent(in) :: nex, Ifl
            real(c_double), intent(in) :: earx(nex+1)
            real(c_double), intent(in) :: params(1)
            real(c_double), intent(inout) :: absorbx(nex), photerx(nex)
            character(kind = c_char), intent(in) :: str(*)
        end subroutine c_tbabs

        subroutine donthcomp(earx, nex, params, Ifl, absorbx, photerx)
            !> The interface around the XSPEC donthcomp. Note that this is a
            !> single-precision Fortran function and so does not need to be
            !> bound to a C symbol. Call it through `nthcomp` below.
            use rtconstants, only: sp
            integer, intent(in) :: nex, Ifl
            real(sp), intent(in) :: earx(nex+1)
            real(sp), intent(in) :: params(5)
            real(sp), intent(inout) :: absorbx(nex), photerx(nex)
        end subroutine donthcomp
    end interface
contains

    subroutine tbabs(earx, nex, nh, Ifl, absorbx, photerx)
        !> Call the tbabs function from the XSPEC model library.
        real(wp), intent(in) :: earx(0:nex), nh
        real(wp), intent(inout) :: absorbx(nex), photerx(nex)
        integer, intent(in) :: nex, Ifl

        call c_tbabs(earx, nex, [nh], Ifl, absorbx, photerx, "")
    end subroutine tbabs

    subroutine nthcomp(earx, nex, params, Ifl, photarx, photerx)
        !> Call the single-precision XSPEC donthcomp model, converting the
        !> arguments to and from the working precision.
        use rtconstants, only: sp
        integer, intent(in) :: nex, Ifl
        real(wp), intent(in) :: earx(0:nex), params(5)
        real(wp), intent(out) :: photarx(nex), photerx(nex)
        real(sp) :: s_earx(0:nex), s_photarx(nex), s_photerx(nex)

        s_earx = real(earx, sp)
        call donthcomp(s_earx, nex, real(params, sp), Ifl, s_photarx,          &
            s_photerx)
        photarx = real(s_photarx, wp)
        photerx = real(s_photerx, wp)
    end subroutine nthcomp

    subroutine table_model(ear, ne, params, filename, ifl, photar)
        !> Interpolate a spectrum from an XSPEC table model file using the
        !> single-precision xsatbl, converting the arguments to and from the
        !> working precision. `filename` must be null terminated.
        use iso_c_binding, only: c_float, c_char
        integer, intent(in) :: ne, ifl
        real(wp), intent(in) :: ear(0:ne), params(:)
        character(kind=c_char), intent(in) :: filename(*)
        real(wp), intent(out) :: photar(ne)
        real(c_float) :: s_photar(ne), s_photer(ne)

        call xsatbl(real(ear, c_float), ne, real(params, c_float), filename,   &
            ifl, s_photar, s_photer)
        photar = real(s_photar, wp)
    end subroutine table_model
end module xspec_interface
