! Wall-clock section timers for profiling (B1 investigation, not for upstream).
module rt_timing
    use iso_c_binding
    implicit none
    integer, parameter :: RT_NTIMERS = 16
    integer, parameter :: T_TOTAL = 1, T_RTRANS = 2, T_TRACE = 3, T_GETDCOS = 4,&
        T_SUM_GR = 5, T_SUM_FLAT = 6, T_GETLENS = 7, T_CONV_TOTAL = 8,         &
        T_RESTFRAME = 9, T_FFTCONV = 10, T_INITCONT = 11, T_POST = 12,          &
        T_GRID = 13
    real(c_double), save :: rt_elapsed(RT_NTIMERS) = 0.d0
    integer(c_int64_t), save :: rt_start(RT_NTIMERS) = 0
    integer(c_int), save :: rt_calls(RT_NTIMERS) = 0
contains
    subroutine tic(k)
        integer, intent(in) :: k
        integer(8) :: c
        call system_clock(c)
        rt_start(k) = c
    end subroutine tic
    subroutine toc(k)
        integer, intent(in) :: k
        integer(8) :: c, rate
        call system_clock(c, rate)
        rt_elapsed(k) = rt_elapsed(k) + dble(c - rt_start(k)) / dble(rate)
        rt_calls(k) = rt_calls(k) + 1
    end subroutine toc
    subroutine rt_timing_reset() bind(C, name="rt_timing_reset")
        rt_elapsed = 0.d0
        rt_calls = 0
    end subroutine rt_timing_reset
    subroutine rt_timing_get(out, calls) bind(C, name="rt_timing_get")
        real(c_double), intent(out) :: out(RT_NTIMERS)
        integer(c_int), intent(out) :: calls(RT_NTIMERS)
        out = rt_elapsed
        calls = rt_calls
    end subroutine rt_timing_get
end module rt_timing
