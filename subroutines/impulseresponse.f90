! vim: cc=80 wrap tw=80
module impulseresponse
    use iso_c_binding, only: c_ptr, c_loc, c_int
    implicit none

    double precision, allocatable, target :: response(:, :)
    double precision, allocatable, target :: time_axis(:)
    double precision, allocatable, target :: energy_axis(:)

    double precision :: time_ax_min = 1.0
    double precision :: time_ax_max = 2.0e3

    integer :: n_ebins, n_tbins
contains

    ! Calculate the size of the axes deltas
    subroutine response_calculate_deltas(dlogt, dg)
        double precision, intent(out) :: dlogt, dg
        dlogt = log10( time_ax_max/time_ax_min ) / float(n_tbins)
        dg = 2.0 / float(n_ebins)
    end subroutine response_calculate_deltas

    ! Populate the axes of the impulse response matrix
    subroutine response_setup_axes(dlogt, dg)
        ! the time axis is logarithmic
        double precision, intent(in) :: dlogt, dg
        double precision :: g
        integer :: i

        do i = 0,n_tbins
            time_axis(i) = time_ax_min * 10.0**( i * dlogt )
        end do

        ! the energy axis is linear
        g = 0.0
        do i = 0,n_ebins
            energy_axis(i) = g
            g = g + dg
        end do
    end subroutine response_setup_axes

    ! Allocate the response matrix and write the delta sizes into the dlogt and
    ! dg arguments
    subroutine response_allocate(ne, nt, dlogt, dg)
        integer, intent(in) :: ne, nt
        double precision, intent(out) :: dlogt, dg

        if (.not. allocated(response)) then
            allocate(response(ne, nt))
            allocate(time_axis(0:nt))
            allocate(energy_axis(0:ne))
            n_ebins = ne
            n_tbins = nt

            call response_calculate_deltas(dlogt, dg)
            call response_setup_axes(dlogt, dg)
        else
            ! Assign in either case, as the values may be needed
            call response_calculate_deltas(dlogt, dg)
        end if

        ! Zero the matrix
        response = 0.0
    end subroutine response_allocate

    ! Zero the edges of the impulse response
    subroutine response_zero_edges()
        integer :: i
        ! Loop over t bins
        do i = 1,n_tbins
            response(1,i)  = 0.0
            response(n_ebins,i) = 0.0
        end do
        ! Loop over g bins
        do i = 1,n_ebins
            response(i,1)  = 0.0
            response(i,n_tbins) = 0.0
        end do
    end subroutine response_zero_edges

    ! External accessor for the response matrix
    subroutine response_get(ptr, eaxis_ptr, taxis_ptr, ne, nt) bind(C)
        type(c_ptr), intent(out) :: ptr, eaxis_ptr, taxis_ptr
        integer(c_int), intent(out) :: ne, nt
        ptr = c_loc(response)
        eaxis_ptr = c_loc(energy_axis)
        taxis_ptr = c_loc(time_axis)
        ne = n_ebins
        nt = n_tbins
    end subroutine response_get

end module impulseresponse
