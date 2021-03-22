!-----------------------------------------------------------------------
subroutine propercross(nex, nf, earx, ReSraw, ImSraw, ReGraw, ImGraw)
  use telematrix
  implicit none
  integer, intent(in)  :: nex, nf
  real,    intent(in)  :: earx(0:nex), ReSraw(nex,nf), ImSraw(nex,nf)
  real,    intent(out) :: ReGraw(nex,nf), ImGraw(nex,nf)
  real,    allocatable :: ReStel(:), ImStel(:)
  real                 :: reref, imref, dum
  integer              :: i, j

!Read from response file
     if( needresp )then
        call initmatrix
     endif

!Allocate arrays
     if( .not. allocated(ReStel) ) allocate(ReStel(numchn))
     if( .not. allocated(ImStel) ) allocate(ImStel(numchn))
     
!Get energy bounds of the reference band
     if( needchans )then
        write(*,*)"Enter lower energy in reference band"
        read(*,*)Elo
        write(*,*)"Enter upper energy in reference band"
        read(*,*)Ehi
        if( Elo .gt. Ehi )then
           dum = Elo
           Elo = Ehi
           Ehi = dum
           write(*,*)"Elo>Ehi! Switched!"
        end if
        Ilo = 1
        Ihi = numchn
        do i = 0, numchn
           if( echn(i) .lt. Elo ) Ilo = i
           if( echn(i) .le. Ehi ) Ihi = i
        end do
        Ilo = Ilo + 1
        if( Ilo .gt. Ihi ) Ihi = Ilo
        needchans = .false.
     end if

     !Calculate `raw' cross-spectrum
     do j = 1, nf
        !Fold around the response matrix
        call cfold(nex, earx, ReSraw(:,j), ImSraw(:,j), ReStel, ImStel)
        !Calcluate reference band
        reref = 0.0
        imref = 0.0
        do i = Ilo, Ihi
           reref = reref + ReStel(i)
           imref = imref + ImStel(i)
        end do
        
        !Cross subject band with reference band
        do i = 1, nex
           ReGraw(i,j) = ReSraw(i,j) * reref + ImSraw(i,j) * imref
           ImGraw(i,j) = ImSraw(i,j) * reref - ReSraw(i,j) * imref
        end do
     end do
     
  return
end subroutine propercross
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine propercross_NOmatrix(nex, nf, earx, ReSraw, ImSraw, ReGraw, ImGraw)
  use telematrix
  implicit none
  integer, intent(in)  :: nex, nf
  real,    intent(in)  :: earx(0:nex), ReSraw(nex,nf), ImSraw(nex,nf)
  real,    intent(out) :: ReGraw(nex,nf), ImGraw(nex,nf)
  real,    allocatable :: ReStel(:), ImStel(:)
  real                 :: reref, imref, dum, dE
  integer              :: i, j


!Get energy bounds of the reference band
     if( needchans )then
        write(*,*)"Enter lower energy in reference band"
        read(*,*)Elo
        write(*,*)"Enter upper energy in reference band"
        read(*,*)Ehi
        if( Elo .gt. Ehi )then
           dum = Elo
           Elo = Ehi
           Ehi = dum
           write(*,*)"Elo>Ehi! Switched!"
        end if
        Ilo = 1
        Ihi = nex
        do i = 0, nex
           if( earx(i) .lt. Elo ) Ilo = i
           if( earx(i) .le. Ehi ) Ihi = i
        end do
        Ilo = Ilo + 1
        if( Ilo .gt. Ihi ) Ihi = Ilo
        needchans = .false.
     end if

     !Calculate `raw' cross-spectrum
     do j = 1, nf
        !Fold around the response matrix
!        call cfold(nex, earx, ReSraw(:,j), ImSraw(:,j), ReStel, ImStel)
        !Calcluate reference band
        reref = 0.0
        imref = 0.0
        do i = Ilo, Ihi
           dE = earx(i) - earx(i-1)
           reref = reref + ReSraw(i,j) * dE
           imref = imref + ImSraw(i,j) * dE
        end do

        !Cross subject band with reference band
        do i = 1, nex
           ReGraw(i,j) = ReSraw(i,j) * reref + ImSraw(i,j) * imref
           ImGraw(i,j) = ImSraw(i,j) * reref - ReSraw(i,j) * imref
        end do
     end do
  return
end subroutine propercross_NOmatrix
!-----------------------------------------------------------------------
