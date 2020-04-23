program main_test
!  make -f revmakefile run
!  make -f revmakefile_laptop run
  implicit none
  integer ne, i, j, ifl
  parameter (ne=10000)
  real Emax, Emin, ear(0:ne), param(19), photar(ne), E, dE
  double precision :: dear(0:ne), dphotar(ne), dphoter(ne)
  double precision :: relxill_par(12), relxill_ion(15)
  real             :: freqmin(5), freqmax(5), single_p(5)
  character (len = 100) :: char 

!managing the time test calculation
  integer :: count,count1,count0,count3,count2,count_max,count_rate

  !----Parameters-------------------
  param(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
  param(2)  = 0.998     !a     !BH spin
  param(3)  = 30.0    !inc   !Inclination angle in degrees
  param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
  param(5)  = 1000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
  param(6)  = 0.0     !zcos  !Cosmological redshift
  param(7)  = 2.0     !Gamma !Photon index
  param(8)  = 3.0    !logxi !log10xi - ionisation parameter
  param(9)  = 1.0     !Afe   !Iron abundance      
  ! param(10) = 300.0   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
  param(10) = 15      !log(n_e) density for xillverD 
  param(11) = 0.0     !Nh
  param(12) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
  param(13) = 10.0    !M     !BH mass in solar masses
  param(14) = 5.0    !flo   !Lowest frequency in band (Hz)
  param(15) = 10.0   !fhi   !Highest frequency in band (Hz)
  param(16) = -4       !ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
  param(17) = 0.0     !DelA
  param(18) = 0.0     !DelAB
  param(19) = 0.1     !gamma
  !---------------------------------
  
  Emax  = 500.0
  Emin  = 0.1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

! ***************************************************************** !
  ! freqmin(1) = 5.0
  ! freqmax(1) = 10.0
  ! freqmin(2) = 15.0
  ! freqmax(2) = 16.0
  ! freqmin(3) = 20.0
  ! freqmax(3) = 21.0
  ! freqmin(4) = 30.0
  ! freqmax(4) = 31.0
  ! freqmin(5) = 50.0
  ! freqmax(5) = 51.0

  ! single_p(1) = 1.4
  ! single_p(2) = 1.6
  ! single_p(3) = 1.8
  ! single_p(4) = 2.0 
  ! single_p(5) = 2.2
  
  ! open(99, file='plots/lag_Nion_fixD19_gamma.dat')
  ! write(99, *) 'skip on'
  
  ! ! param(1)  = 2.0
  ! ! param(4)  = 2.0
  ! param(10) = 19
  ! do j = 1, 5
  !    ! param(10) = param(10) + 1
  !    ! param(14) = freqmin(j)
  !    ! param(15) = freqmax(j)
  !    param(7) = single_p(j)
  !    write(*,*) 'h: '    , param(1) 
  !    write(*,*) 'rin: '  , param(4) 
  !    write(*,*) 'LogNe: ', param(10) 
  !    write(*,*) 'freq: ' , param(14), param(15)
  !    call tdreltransD(ear, ne, param, ifl, photar)
  !    do i = 1, ne
  !       E  = 0.5 * ( ear(i) + ear(i-1) )
  !       dE = ear(i) - ear(i-1)
  !       if( abs(param(16)) .gt. 3.5 .and. abs(param(16)) .lt. 4.5 .or. abs(param(16)) .gt. 5.5 )then
  !          write(99,*)E, photar(i)/dE
  !       else
  !          write(99,*)E, E**2*photar(i)/dE
  !       end if
  !    end do
  !    ! param(1)  = param(1) + 10.0
  !    ! param(4)  = param(4) + 10.
  !     ! param(19) = param(19) + 0.05
  !    write(99, *) 'no no'

  ! enddo

  ! write(99, *) 'scr white'
  ! write(99, *) 'log x on'
  ! write(99, *) 'lw 5'
  ! write(99, *) 'la y Lag (s)'
  ! write(99, *) 'la x Energy [keV]'
  ! close(99)
!***************************************************************** !


!***************************************************************** !


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
  param(10) = 19      !log(n_e)
  param(11) = 0.0     !Nh
  param(12) = 1.0     !1onB  
  param(13) = 10.0    !M     
  param(14) = 0.0     !flo   
  param(15) = 0.0    !fhi   
  param(16) = 1      !ReIm 
  param(17) = 0.0     !DelA
  param(18) = 0.0     !DelAB
  param(19) = 0.0     !gamma
  ! relxill_par(1)  = dble(param(1) )
  ! relxill_par(2)  = dble(param(2) )
  ! relxill_par(3)  = dble(param(3) )
  ! relxill_par(4)  = dble(param(4) )
  ! relxill_par(5)  = dble(param(5) )
  ! relxill_par(6)  = dble(param(6) )
  ! relxill_par(7)  = dble(param(7) )
  ! relxill_par(8)  = dble(param(8) )
  ! relxill_par(9)  = dble(param(9) )
  ! relxill_par(10) = dble(param(10)) 
  ! relxill_par(11) = dble(param(12)) 
  ! relxill_par(12) = 1.d0 

  single_p(1) = 15 
  single_p(2) = 30 
  single_p(3) = 45
  single_p(4) = 60  
  single_p(5) = 75
  
  ! open(99, file='plots/DC_reltransD_ion1_D19_incl2.dat')
  open(99, file='../check_plots/DC_reltransD_zones.dat')
  write(99, *) 'skip on'
  ! open(98, file='plots/DC_relxilllpion.dat')
  ! write(98, *) 'skip on'

  param(8)  = 4.0
  ! param(10) = 19
  ! do j = 1, 5
  !    param(3) = single_p(j)  
     write(*,*) 'LogNe', param(10) 
     ! write(*,*) 'Inclination', param(3)
     call tdreltransD(ear,ne,param,ifl,photar)
     do i = 1, ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE = ear(i) - ear(i-1)
        if( abs(param(16)) .gt. 3.5 .and. abs(param(16)) .lt. 4.5 .or. abs(param(16)) .gt. 5.5 )then
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

  ! enddo

  ! write(99, *) 'scr white'
  ! write(99, *) 'log x y on'
  ! write(99, *) 'lw 5'
  ! write(99, *) 'la y Flux []'
  ! write(99, *) 'la x Energy [keV]'
  ! close(99)
!***************************************************************** !




!***************************************************************** !
  ! param(1)  = 6.0     !h     
  ! param(2)  = 0.998   !a   
  ! param(3)  = 30.0    !inc   
  ! param(4)  = -1.0    !rin   
  ! param(5)  = 1000.0  !rout  
  ! param(6)  = 0.0     !zcos  
  ! param(7)  = 2.0     !Gamma 
  ! param(8)  = 4.0     !logxi 
  ! param(9)  = 1.0     !Afe   
  ! param(10) = 15.0      !log(n_e)
  ! param(11) = 0.0     !Nh
  ! param(12) = 1.0     !1onB  
  ! param(13) = 10.0    !M     
  ! param(14) = 0.0     !flo   
  ! param(15) = 0.0    !fhi   
  ! param(16) = 1      !ReIm 
  ! param(17) = 0.0     !DelA
  ! param(18) = 0.0     !DelAB
  ! param(19) = 0.0     !gamma
  ! relxill_ion(1)  = dble(param(1) )
  ! relxill_ion(2)  = dble(param(2) )
  ! relxill_ion(3)  = dble(param(3) )
  ! relxill_ion(4)  = dble(param(4) )
  ! relxill_ion(5)  = dble(param(5) )
  ! relxill_ion(6)  = dble(param(6) )
  ! relxill_ion(7)  = dble(param(7) )
  ! relxill_ion(8)  = dble(param(8) )
  ! relxill_ion(9)  = dble(param(9) )
  ! relxill_ion(10) = 300.0 ! Ecut 
  ! relxill_ion(11) = 0.d0  ! velocity of the source  
  ! relxill_ion(12) = 2.d0  ! ion_grad_type
  ! relxill_ion(13) = 1.d0  ! xi_index
  ! relxill_ion(14) = 1.d0  ! refl_frac
  ! relxill_ion(15) = 1.d0  ! fixReflFrac
  
  ! open(99, file='plots/DC_reltransD_fixD_ion0_100zones.dat')
  ! write(99, *) 'skip on'
  ! open(98, file='plots/DC_relxilllpion.dat')
  ! write(98, *) 'skip on'

  ! param(8) = 2.0
  ! do j = 1, 5
  !    param(8) = param(8) + 0.5 
  !    write(*,*) 'LogXi', param(8) 
  !    call tdreltransD(ear,ne,param,ifl,photar)
  !    do i = 1, ne
  !       E  = 0.5 * ( ear(i) + ear(i-1) )
  !       dE = ear(i) - ear(i-1)
  !       if( abs(param(16)) .gt. 3.5 .and. abs(param(16)) .lt. 4.5 .or. abs(param(16)) .gt. 5.5 )then
  !          write(99,*)E, photar(i)/dE
  !       else
  !          write(99,*)E, E**2*photar(i)/dE
  !       end if
  !    end do
  !    write(99, *) 'no no'

     ! dear = dble(ear)
     ! relxill_ion(8) = dble(param(8)) 
     ! call lmodrelxilllpionf(dear, ne, relxill_ion, ifl, dphotar, dphoter, char)

     ! do i = 1, ne 
     !    E  = 0.5 * ( ear(i) + ear(i-1) )
     !    dE = ear(i) - ear(i-1)
     !    write(98, *) E, E**2 * dphotar(i)/ dE
     ! enddo
     ! write(98, *) 'no no'

  ! enddo
!***************************************************************** !


  
!Call directly relxilllp to compare with our model 
!***************************************************************** !

  ! relxill_par(1)  = dble(param(1) )
  ! relxill_par(2)  = dble(param(2) )
  ! relxill_par(3)  = dble(param(3) )
  ! relxill_par(4)  = dble(param(4) )
  ! relxill_par(5)  = dble(param(5) )
  ! relxill_par(6)  = dble(param(6) )
  ! relxill_par(7)  = dble(param(7) )
  ! relxill_par(8)  = dble(param(8) )
  ! relxill_par(9)  = dble(param(9) )
  ! relxill_par(10) = 300.d0 
  ! relxill_par(11) = dble(param(12)) 
  ! relxill_par(12) = 1.d0 
  

  ! dear = dble(ear)
  ! call lmodrelxilllpf(dear, ne, relxill_par, ifl, dphotar, dphoter, char)

  ! do i = 1, ne 
  !    E  = 0.5 * ( ear(i) + ear(i-1) )
  !    dE = ear(i) - ear(i-1)
  !    write(88, *) E, E**2 * dphotar(i)/ dE
  ! enddo

! ***************************************************************** !



  
end program main_test
