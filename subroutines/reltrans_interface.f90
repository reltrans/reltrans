module reltrans_interface
!> This module defines the library interface for reltrans.
    contains
        subroutine wrap_get_raytrace_coords(p,f1234,lambda,q,sinobs,muobs,     &
             a_spin,robs, scal,radi,mu,phi,time,sigma)                         &
             bind(C, name = "get_raytrace_coords")
            use raytracing, only: get_raytrace_coords
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            double precision, intent(in)  :: p
            double precision, intent(in)  :: f1234(4), lambda, q, sinobs, muobs
            double precision, intent(in)  :: a_spin, robs, scal
            double precision, intent(out) :: radi, mu, phi, time, sigma
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            call get_raytrace_coords(p,f1234,lambda,q,sinobs,muobs,a_spin,     &
                 robs, scal,radi,mu,phi,time,sigma)
            return
        end subroutine wrap_get_raytrace_coords

        subroutine wrap_p_disk_crossing(f1234,lambda,q,sins,mus,a_spin,h,      &
          scal,mudisk,r_max,r_min, output) bind(C, name = "p_disk_crossing")
            use raytracing, only: p_disk_crossing
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            double precision, intent(in) :: f1234(4), lambda, q, sins, mus
            double precision, intent(in) :: a_spin, h, scal, mudisk, r_max, r_min
            double precision, intent(out):: output
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            output = p_disk_crossing(f1234,lambda,q,sins,mus,a_spin,h,         &
               scal,mudisk,r_max,r_min)
            return
        end subroutine wrap_p_disk_crossing
        
        subroutine wrap_p_coord_at_infinity(f1234,lambda,q,sins,mus,a_spin,    &
          h,scal, output) bind(C, name = "p_coord_at_infinity")
            use raytracing, only: p_coord_at_infinity
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            double precision, intent(in) :: f1234(4), lambda, q, sins, mus
            double precision, intent(in) :: a_spin, h, scal
            double precision, intent(out):: output(4)
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            output = p_coord_at_infinity(f1234,lambda,q,sins,mus,a_spin,h,scal)
            return
        end subroutine wrap_p_coord_at_infinity

        subroutine wrap_initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,    &
            lambda,q,f1234) bind(C, name = "initial_photon")
            use raytracing, only: initial_photon
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            double precision, intent(in)  :: pr, pt, pp, sins, mus, a_spin, h
            double precision, intent(in)  :: velocity(3)
            double precision, intent(out) :: lambda, q
            double precision, intent(out) :: f1234(4)
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            call initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,           &
            lambda,q,f1234)
            return
        end subroutine wrap_initial_photon

        subroutine wrap_constants_of_motion(alpha,beta,robs,sinobs,muobs,      &
            a_spin,scal, velocity,f1234,lambda,q)                              &
            bind(C, name = "constants_of_motion")
            use raytracing, only: constants_of_motion
            use kerrz, only: kerr_metric, krz_KerrMetric_init
            implicit none
            double precision, intent(in) :: alpha, beta, robs, sinobs, muobs
            double precision, intent(in) :: a_spin, scal
            double precision, intent(in) :: velocity(3)
            double precision, intent(out) :: f1234(4), lambda, q
            kerr_metric = krz_KerrMetric_init(1.0d0, a_spin)
            call constants_of_motion(alpha,beta,robs,sinobs,muobs,             &
            a_spin,scal, velocity,f1234,lambda,q)
            return
        end subroutine wrap_constants_of_motion
        
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
            double precision, intent(in) :: a_spin, h, muobs
            double precision, intent(inout) :: lens, del_t, cosdelta
            call getlens(a_spin, h, muobs, lens, del_t, cosdelta)
        end subroutine wrap_getlens

          
end module reltrans_interface
