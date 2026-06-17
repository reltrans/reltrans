module raytracing
    implicit none

contains

    subroutine raytrace(model_args)
    !> This subroutine is a shim function that allows for the substitution of the
    !> raytracing functionality without adjusting the rest of the code
    !> Inputs:
    !>   model_args: a real array containing the model parameters
    !> Outputs:
    !>   none
        use grtrace, only: GRtrace
        implicit none
        real, intent(in) :: model_args(:)
        ! Call the GRtrace subroutine
        call GRtrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
        return
    end subroutine raytrace


