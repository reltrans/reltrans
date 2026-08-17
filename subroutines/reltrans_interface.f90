module reltrans_interface
!> This module defines the library interface for reltrans.
    contains
        subroutine wrap_trace_disk_observer(nro,nphi,rn,mueff,mu0,spin,rmin,   &
                   rout,mudisk,d) bind(C, name = "trace_disk_observer")
            use raytracing, only: trace_disk_observer
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            integer, intent(in) :: nro, nphi
            double precision, intent(in) :: rn(nro), mueff, mu0, spin, rmin, rout
            double precision, intent(in) :: mudisk, d
            kerr_metric = krz_KerrMetric_init(1.0d0, spin)
            call trace_disk_observer(nro,nphi,rn,mueff,mu0,spin,rmin,rout,     &
                                    mudisk,d)
            return
        end subroutine wrap_trace_disk_observer

        subroutine wrap_getlens(a_spin, h, muobs, lens, del_t, cosdelta)       &
            bind(C, name = "wrap_getlens")
            use raytracing, only: getlens
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            double precision, intent(in) :: a_spin, h, muobs
            double precision, intent(inout) :: lens, del_t, cosdelta
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            call getlens(a_spin, h, muobs, lens, del_t, cosdelta)
        end subroutine wrap_getlens

end module reltrans_interface
