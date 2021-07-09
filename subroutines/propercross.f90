!-----------------------------------------------------------------------
subroutine propercross(nex, nf, earx, ReSraw, ImSraw, ReGraw, ImGraw, resp_matr)
  use telematrix
  use telematrix2
  implicit none
  integer, intent(in)  :: nex, nf, resp_matr
  real,    intent(in)  :: earx(0:nex), ReSraw(nex,nf), ImSraw(nex,nf)
  real,    intent(out) :: ReGraw(nex,nf), ImGraw(nex,nf)
  real,    allocatable :: ReStel(:), ImStel(:)
  real,    allocatable :: ReStel2(:), ImStel2(:)
  real                 :: reref, imref
  integer              :: i, j


  call response_and_energy_bounds(resp_matr)

  
  if (resp_matr .eq. 1) then 
     !Allocate arrays
     if( .not. allocated(ReStel) ) allocate(ReStel(numchn))
     if( .not. allocated(ImStel) ) allocate(ImStel(numchn))

     !Calculate `raw' cross-spectrum
     do j = 1, nf
        !Fold around the response matrix

        ! write(*,*) 'Folding the first matrix'
        call cfold(nex, earx, ReSraw(:,j), ImSraw(:,j), ReStel, ImStel)
        ! write(*,*) 'finished'
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
  endif
  

  if (resp_matr .eq. 2) then 
     !Allocate arrays
     if( .not. allocated(ReStel2) ) allocate(ReStel2(numchn2))
     if( .not. allocated(ImStel2) ) allocate(ImStel2(numchn2))
     do j = 1, nf
 !Fold around the response matrix
        ! write(*,*) 'Folding the second matrix'
        call cfold2(nex, earx, ReSraw(:,j), ImSraw(:,j), ReStel2, ImStel2)
        ! write(*,*) 'finished'
        !Calcluate reference band
        reref = 0.0
        imref = 0.0
        ! write(*,*) 'Calculating the refenrece band with the second matrix'
        do i = Ilo2, Ihi2
           reref = reref + ReStel2(i)
           imref = imref + ImStel2(i)
        end do
        ! write(*,*) 'finished'

        !Cross subject band with reference band
        ! write(*,*) 'Cross spectrum with the refence band (second matrix)'
        do i = 1, nex
           ReGraw(i,j) = ReSraw(i,j) * reref + ImSraw(i,j) * imref
           ImGraw(i,j) = ImSraw(i,j) * reref - ReSraw(i,j) * imref
        end do
        ! write(*,*) 'finished'
     end do
  endif
  
  return
end subroutine propercross
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine response_and_energy_bounds(resp_matr)
  use telematrix
  use telematrix2
  implicit none
  integer, INTENT(IN) :: resp_matr
  
  real     :: dum
  integer  :: i
  
!Read from response file
  if(resp_matr .eq. 1) then 
     if( needresp )then
        ! write(*,*) 'calling intmatrix'
        call initmatrix
     endif
     
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

  else if (resp_matr .eq. 2) then
     if( needresp2 )then
        ! write(*,*) 'calling intmatrix2'
        call initmatrix2
     endif
!second response matrix     
     if( needchans2 )then
        write(*,*)"Enter lower energy in reference band of the second response"
        read(*,*)Elo2
        write(*,*)"Enter upper energy in reference band of the second response"
        read(*,*)Ehi2
        if( Elo2 .gt. Ehi2 )then
           dum  = Elo2
           Elo2 = Ehi2
           Ehi2 = dum
           write(*,*)"Elo>Ehi! Switched!"
        end if
        Ilo2 = 1
        Ihi2 = numchn2
        do i = 0, numchn2
           if( echn2(i) .lt. Elo2 ) Ilo2 = i
           if( echn2(i) .le. Ehi2 ) Ihi2 = i
        end do
        Ilo2 = Ilo2 + 1
        if( Ilo2 .gt. Ihi2 ) Ihi2 = Ilo2
        needchans2 = .false.
     end if

  else
     write(*,*) 'MORE THAN 2 RESPONSES NOT YET IMPLEMENTED... TAKES THE FIRST ONE'
     if( needresp )then
        call initmatrix
     endif
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
  endif


  
  
   end subroutine response_and_energy_bounds


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
