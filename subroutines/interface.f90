module interface
!> This module defines the library interface for reltrans.
    contains

        subroutine wrap_disk_observer_trace(nro,nphi,rn,mueff,mu0,spin,rmin,   &
                                            rout,mudisk,d)&
            bind(C, name = "grtrace")
            use raytracing, only: disk_observer_trace
            integer, intent(in) :: nro, nphi
            double precision, intent(in) :: rn(nro), mueff, mu0, spin, rmin, rout
            double precision, intent(in) :: mudisk, d
            call disk_observer_trace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,     &
                                    mudisk,d)
            return
        end subroutine wrap_disk_observer_trace
end module interface