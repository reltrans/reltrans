program main_simple_call
 implicit none
 integer ne,i,ifl
 parameter (ne=6000)
 real Emax,Emin,ear(0:ne),param(19),photar(ne),E,dE
 real dist_par(24)
 
!----Parameters of distance model-------------------
  dist_par(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
  dist_par(2)  = 0.9     !a     !BH spin
  dist_par(3)  = 30.0    !inc   !Inclination angle in degrees
  dist_par(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
  dist_par(5)  = 20000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
  dist_par(6)  = 0.0247  !zcos  !Cosmological redshift
  dist_par(7)  = 2.0     !Gamma !Photon index
  dist_par(8)  = 1e5 !3.0     !Dkpc  !Distance in kpc
  dist_par(9)  = 1.0     !Afe   !Iron abundance
  dist_par(10) = 15.0  !20.0  !logne !Log10 of minimum disc electron number density
  dist_par(11) = 300.0   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
  dist_par(12) = 0.0     !Nh
  dist_par(13) = 1.0     !1onB  !(1/\mathcal{B}): "proper" boosting factor
  dist_par(14) = 1e7  !10.0  !M     !BH mass in solar masses
  dist_par(15) = 0.01    !h/r   !Disc scaleheight
  dist_par(16) = 0.0 !0.429   !b1    !Linear term in angular emissivity function
  dist_par(17) = 0.0 !-2.0    !b2    !Quadratic term in angular emissivity function
  dist_par(18) = 3e-5 !30.0    !flo   !Lowest frequency in band (Hz)
  dist_par(19) = 1e-4 !100.0   !fhi   !Highest frequency in band (Hz)
  dist_par(20) = 4       !ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
  dist_par(21) = 0.0     !DelA
  dist_par(22) = 0.0     !DelAB
  dist_par(23) = 0.0     !gamma
  dist_par(24) = 6.615e-05 !0.0735  !Anorm !Normalisation -- takes the place of the XSPEC norm (fix that to 1)
!---------------------------------

! Set energy grid
  Emax  = 500.0
  Emin  = 0.1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

! Call the distance model
  call tdrtdist(ear, ne, dist_par, ifl, photar)
 
! Write out
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     if( dist_par(20) .gt. 3.5 .and. dist_par(20) .lt. 4.5 .or. dist_par(20) .gt. 5.5 )then
        write(99,*)E,photar(i)/dE
     else
        write(99,*)E,E**2*photar(i)/dE
     end if
  end do
  write(99,*)"no no" 
 
! Write out DC component
  dist_par(18) = 0.0
  dist_par(19) = 0.0
  dist_par(20) = 1.0
  call tdrtdist(ear, ne, dist_par, ifl, photar)
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     if( dist_par(20) .gt. 3.5 .and. dist_par(20) .lt. 4.5 .or. dist_par(20) .gt. 5.5 )then
        write(99,*)E,photar(i)/dE
     else
        write(99,*)E,E**2*photar(i)/dE
     end if
  end do
  write(99,*)"no no" 
 
!  !----Parameters-------------------
!   param(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
!   param(2)  = 0.9     !a     !BH spin
!   param(3)  = 30.0    !inc   !Inclination angle in degrees
!   param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
!   param(5)  = 20000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
!   param(6)  = 0.0     !zcos  !Cosmological redshift
!   param(7)  = 2.0     !Gamma !Photon index
!   param(8)  = 3.75    !logxi !log10xi - ionisation parameter
!   param(9)  = 1.0     !Afe   !Iron abundance      
!   param(10) = 300.0   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
!   param(11) = 0.0     !Nh
!   param(12) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
!   param(13) = 10.0    !M     !BH mass in solar masses
!   param(14) = 30.0    !flo   !Lowest frequency in band (Hz)
!   param(15) = 100.0   !fhi   !Highest frequency in band (Hz)
!   param(16) = 4       !ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
!   param(17) = 0.0     !DelA
!   param(18) = 0.0     !DelAB
!   param(19) = 0.0     !gamma
!   !---------------------------------
!   
! 
! 
! ! Call one of the existing models
!   call tdreltransCp(ear,ne,param,ifl,photar)
  ! do i = 1,ne
  !    E  = 0.5 * ( ear(i) + ear(i-1) )
  !    dE = ear(i) - ear(i-1)
  !    if( param(16) .gt. 3.5 .and. param(16) .lt. 4.5 .or. param(16) .gt. 5.5 )then
  !       write(99,*)E,photar(i)/dE
  !    else
  !       write(99,*)E,E**2*photar(i)/dE
  !    end if
  ! end do
  ! write(99,*)"no no"

  
  write(*,*)"---------------------------------------------------------"
  write(*,*)"Outputs:"
  write(*,*)"fort.188: r,logxi,logxieff"
  write(*,*)"fort.99: spectrum / lag spectrum etc"
  write(*,*)"---------------------------------------------------------"

  
end program main_simple_call
