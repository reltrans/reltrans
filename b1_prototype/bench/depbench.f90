program depbench
    implicit none
    integer, parameter :: nex = 4096, xe = 20
    integer :: nf, n, i, j, k, m, klo, khi, rb, it
    complex, allocatable :: w0(:,:,:,:,:), w1(:,:,:,:,:), w2(:,:,:,:,:), w3(:,:,:,:,:)
    double precision, allocatable :: x(:,:), tau(:), wt(:), fi(:)
    integer, allocatable :: gb(:), rbs(:)
    double precision :: t0, t1, x1, x2, x3, d31, d21, d32, c0, c1, frac, y, pi, lg
    complex :: cexp, cw
    real :: r
    character(len=16) :: arg
    pi = acos(-1.d0)
    call get_command_argument(1, arg); read(arg, *) nf
    n = 400000
    allocate(w0(1,nex,nf,1,xe), w1(1,nex,nf,1,xe), w2(1,nex,nf,1,xe), w3(1,nex,nf,1,xe))
    allocate(x(3,n), tau(n), wt(n), gb(n), rbs(n), fi(nf))
    w0 = 0; w1 = 0; w2 = 0; w3 = 0
    do i = 1, nf; fi(i) = 1d-3 * i; end do
    do i = 1, n
        call random_number(r); x(1,i) = 1000 + 2000*r
        call random_number(r); x(2,i) = x(1,i) + 3*r
        call random_number(r); x(3,i) = x(1,i) + 3*r
        call random_number(r); tau(i) = 100*r; wt(i) = 1.d0
        gb(i) = int(x(1,i)); call random_number(r); rbs(i) = 1 + int(19.99*r)
    end do
    ! (a) nearest deposit, as sum_multiple_lampposts (per sample: nf cexp, 4 adds)
    call cpu_time(t0)
    do i = 1, n
        lg = 0.1d0
        do j = 1, nf
            cexp = cmplx(cos(real(2.d0*pi*tau(i)*fi(j))), sin(real(2.d0*pi*tau(i)*fi(j))))
            w0(1,gb(i),j,1,rbs(i)) = w0(1,gb(i),j,1,rbs(i)) + real(wt(i))*cexp
            w1(1,gb(i),j,1,rbs(i)) = w1(1,gb(i),j,1,rbs(i)) + real(lg)*real(wt(i))*cexp
            w2(1,gb(i),j,1,rbs(i)) = w2(1,gb(i),j,1,rbs(i)) + real(wt(i))*cexp
            w3(1,gb(i),j,1,rbs(i)) = w3(1,gb(i),j,1,rbs(i)) + real(wt(i))*cexp
        end do
    end do
    call cpu_time(t1)
    print '(a,i4,a,f10.2,a)', 'nf=', nf, ' nearest: ', (t1-t0)/n*1d9, ' ns/sample'
    ! (b) tent deposit
    call cpu_time(t0)
    m = 0
    do i = 1, n
        x1 = minval(x(:,i)); x3 = maxval(x(:,i)); x2 = sum(x(:,i)) - x1 - x3
        d31 = max(x3-x1,1d-9); d21 = max(x2-x1,5d-10); d32 = max(x3-x2,5d-10)
        klo = floor(x1) + 1; khi = max(ceiling(x3), klo)
        rb = rbs(i)
        do j = 1, nf
            cexp = cmplx(cos(real(2.d0*pi*tau(i)*fi(j))), sin(real(2.d0*pi*tau(i)*fi(j))))
            cw = real(wt(i)) * cexp
            c0 = 0.d0
            do k = klo, khi
                y = min(max(dble(k), x1), x3)
                if (y .le. x2) then; c1 = (y-x1)**2/(d31*d21); else; c1 = 1.d0 - (x3-y)**2/(d31*d32); end if
                frac = c1 - c0; c0 = c1
                w0(1,k,j,1,rb) = w0(1,k,j,1,rb) + real(frac)*cw
                w1(1,k,j,1,rb) = w1(1,k,j,1,rb) + real(frac*0.1d0)*cw
                w2(1,k,j,1,rb) = w2(1,k,j,1,rb) + real(frac)*cw
                w3(1,k,j,1,rb) = w3(1,k,j,1,rb) + real(frac)*cw
                if (j == 1) m = m + 1
            end do
        end do
    end do
    call cpu_time(t1)
    print '(a,i4,a,f10.2,a,f6.2)', 'nf=', nf, ' tent: ', (t1-t0)/n*1d9, ' ns/triangle, bins/tri=', dble(m)/n
    print *, real(w0(1,1500,1,1,3))
end program depbench
