program main_simple_call
 implicit none
 integer ne,i,ifl,neb
 parameter (ne=6000,neb=25)
 real Emax,Emin,ear(0:ne),param(19),photar(ne),E,dE
 real dist_par(24),spar(26),earb(0:neb)
 
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
 
! Set broad energy grid for simulation
  Emax  = 10.0
  Emin  = 0.3
  do i = 0,neb
    earb(i) = Emin * (Emax/Emin)**(real(i)/real(neb))
  end do 
 
!----Parameters of simulation model-------------------
  spar(1)  = 25.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
  spar(2)  = 0.998     !a     !BH spin
  spar(3)  = 37.0    !inc   !Inclination angle in degrees
  spar(4)  = 1.87  !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
  spar(5)  = 20000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
  spar(6)  = 0.0!0.0247  !zcos  !Cosmological redshift
  spar(7)  = 1.969     !Gamma !Photon index
  spar(8)  = 8.6     !Dkpc  !Distance in kpc
  spar(9)  = 0.644     !Afe   !Iron abundance
  spar(10) = 18.927  !20.0  !logne !Log10 of minimum disc electron number density
  spar(11) = 34.73   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
  spar(12) = 8.0     !Nh
  spar(13) = 3.55     !1onB  !(1/\mathcal{B}): "proper" boosting factor
  spar(14) = 12.4  !M     !BH mass in solar masses
  spar(15) = 0.054    !h/r   !Disc scaleheight
  spar(16) = 0.0 !0.429   !b1    !Linear term in angular emissivity function
  spar(17) = 0.0 !-2.0    !b2    !Quadratic term in angular emissivity function
  spar(18) = 10.0!3e-5 !30.0    !flo   !Lowest frequency in band (Hz)
  spar(19) = 30.0!1e-4 !100.0   !fhi   !Highest frequency in band (Hz)
  spar(20) = 0.8      !gamma_c^2 -- squared cohenrence
  spar(21) = 0.0     !DelA
  spar(22) = 0.0     !DelAB
  spar(23) = 0.0     !gamma
  spar(24) = 2.40398E-02  !Anorm !Normalisation -- takes the place of the XSPEC norm (fix that to 1)
  spar(25) = 130.0e3   !Exposure time (s)
  spar(26) = 0.0025    !Power in [rms/mean]^2/Hz units ( |A(nu)|^2/A_0^2 )ss
!--------------------------------- 
  call simrtdist(earb, neb, spar, ifl, photar) 
 
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
