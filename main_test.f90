program main_test
!  make -f revmakefile run
!  make -f revmakefile_laptop run
  implicit none
  integer ne, i, j, ifl
  parameter (ne=10000)
  real Emax, Emin, ear(0:ne), param(21), photar(ne), E, dE
  double precision :: dear(0:ne), dphotar(ne), dphoter(ne)
  double precision :: relxill_par(12), relxill_ion(15)
  real             :: freqmin(5), freqmax(5), single_p(7)
  character (len = 100) :: char 

!managing the time test calculation
  integer :: count,count1,count0,count3,count2,count_max,count_rate

  
  Emax  = 500.0
  Emin  = 0.1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do


!Call directly relxilllp to compare with our model 
!***************************************************************** !

  param(1)  = 6.0     !h     
  param(2)  = 0.998   !a   
  param(3)  = 30.0    !inc   
  param(4)  = -1.0    !rin   
  param(5)  = 1000.0  !rout  
  param(6)  = 0.0     !zcos  
  param(7)  = 2.0     !Gamma 
  param(8)  = 3.0     !logxi 
  param(9)  = 1.0     !Afe   
  param(10) = 15      !log(n_e)
  param(11) = 60      !kTe
  param(12) = 0.1      !kT_bb
  param(13) = 0.0     !Nh
  param(14) = 1.0     !1onB  
  param(15) = 10.0    !M     
  param(16) = 5.0     !flo   
  param(17) = 10.0    !fhi   
  param(18) = -4      !ReIm 
  param(19) = 0.0     !DelA
  param(20) = 0.0     !DelAB
  param(21) = 0.0     !gamma

  single_p(1) = 15 
  single_p(2) = 16 
  single_p(3) = 17
  single_p(4) = 18  
  single_p(5) = 19
  single_p(6) = 20
  single_p(7) = 21
  
  ! open(99, file='plots/DC_reltransD_ion1_D19_incl2.dat')
  open(99, file='../plots/reltransx/reltransx_Nion_xi_dens')
  write(99, *) 'skip on'
  ! open(98, file='plots/DC_relxilllpion.dat')
  ! write(98, *) 'skip on'

  ! param(8)  = 4.0
  ! param(10) = 19
  do j = 1, 8
  !    param(3) = single_p(j)  
     write(*,*) 'LogNe', param(10) 
     ! write(*,*) 'Inclination', param(3)
     call tdreltransx(ear,ne,param,ifl,photar)
     do i = 1, ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE = ear(i) - ear(i-1)
        if( abs(param(18)) .gt. 3.5 .and. abs(param(18)) .lt. 4.5 .or. abs(param(16)) .gt. 5.5 )then
           write(99,*)E, photar(i)/dE
        else
           write(99,*)E, E**2*photar(i)/dE
        end if
     end do
     write(99, *) 'no no'

  !    dear = dble(ear)
  !    relxill_par(10) = dble(param(10)) 
     !    call lmodrelxilllpdensf(dear, ne, relxill_par, ifl, dphotar, dphoter, char)

  !    do i = 1, ne 
  !       E  = 0.5 * ( ear(i) + ear(i-1) )
  !       dE = ear(i) - ear(i-1)
  !       write(98, *) E, E**2 * dphotar(i)/ dE
  !    enddo
  !    write(98, *) 'no no'

  enddo

  ! write(99, *) 'scr white'
  ! write(99, *) 'log x y on'
  ! write(99, *) 'lw 5'
  ! write(99, *) 'la y Flux []'
  ! write(99, *) 'la x Energy [keV]'
  ! close(99)
!***************************************************************** !
  
end program main_test
