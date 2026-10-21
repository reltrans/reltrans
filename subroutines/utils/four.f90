module fourier_transform
    !> Numerical Recipes radix-2 complex FFT `four1` (Press et al., Numerical
    !> Recipes in Fortran 77, 2nd ed., section 12.2). Used by `E_FT` and
    !> `E_invFT` in `padcnv.f90`.
    use rtconstants, only: wp
    implicit none
    private
    public :: four1

contains

    subroutine four1(data, nn, isign)
        !> In-place FFT of `nn` complex values stored as interleaved
        !> (real, imaginary) pairs. `nn` must be a power of 2 (not checked).
        !>
        !> Computes X(k) = sum_j x(j) exp(isign * 2 pi i j k / nn), with no
        !> normalisation. Output is in standard FFT order.
        integer, intent(in) :: nn, isign
        real(wp), intent(inout) :: data(2*nn)

        call bit_reverse_reorder(data, nn)
        call danielson_lanczos(data, nn, isign)
    end subroutine four1

    subroutine bit_reverse_reorder(data, nn)
        !> Moves each complex element to the position given by reversing the
        !> bits of its index.
        integer, intent(in) :: nn
        real(wp), intent(inout) :: data(2*nn)
        real(wp) :: swap(2)
        integer :: i, j, m

        ! i and j index the real parts, so both step in 2s.
        j = 1
        do i = 1, 2*nn, 2
            if (j > i) then
                swap = data(j:j+1)
                data(j:j+1) = data(i:i+1)
                data(i:i+1) = swap
            end if
            ! Increment j in bit-reversed order.
            m = nn
            do while (m >= 2 .and. j > m)
                j = j - m
                m = m / 2
            end do
            j = j + m
        end do
    end subroutine bit_reverse_reorder

    subroutine danielson_lanczos(data, nn, isign)
        !> Combines pairs of length-L transforms into length-2L transforms,
        !> for L = 1, 2, 4, ..., nn/2.
        integer, intent(in) :: nn, isign
        real(wp), intent(inout) :: data(2*nn)
        ! 2 pi as written in the Numerical Recipes code.
        real(wp), parameter :: TWO_PI = 6.28318530717959_wp
        real(wp) :: theta, wr, wi, wpr, wpi, wtemp, tempr, tempi
        integer :: n, mmax, istep, m, i, j

        n = 2*nn
        ! Length of each sub-transform in array elements.
        mmax = 2
        do while (n > mmax)
            istep = 2*mmax
            theta = TWO_PI / (isign*mmax)
            ! Twiddle factor w = wr + i wi, advanced by the recurrence
            ! w <- w exp(i theta). wpr = cos(theta) - 1, written to avoid
            ! cancellation for small theta.
            wpr = -2.0_wp * sin(0.5_wp*theta)**2
            wpi = sin(theta)
            wr = 1.0_wp
            wi = 0.0_wp
            do m = 1, mmax, 2
                do i = m, n, istep
                    j = i + mmax
                    tempr = wr*data(j) - wi*data(j+1)
                    tempi = wr*data(j+1) + wi*data(j)
                    data(j) = data(i) - tempr
                    data(j+1) = data(i+1) - tempi
                    data(i) = data(i) + tempr
                    data(i+1) = data(i+1) + tempi
                end do
                wtemp = wr
                wr = wr*wpr - wi*wpi + wr
                wi = wi*wpr + wtemp*wpi + wi
            end do
            mmax = istep
        end do
    end subroutine danielson_lanczos

end module fourier_transform