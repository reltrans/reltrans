program main_test
!  make -f revmakefile run
  implicit none
  integer ne,i,ifl
  parameter (ne=6000)
  real Emax,Emin,ear(0:ne),param(20),photar(ne),E,dE

  !----Parameters-------------------
  param(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
  param(2)  = 0.9     !a     !BH spin
  param(3)  = 30.0    !inc   !Inclination angle in degrees
  param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
  param(5)  = 20000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
  param(6)  = 0.0     !zcos  !Cosmological redshift
  param(7)  = 2.0     !Gamma !Photon index
  param(8)  = 3.75     !logxi !log10xi - ionisation parameter
  param(9)  = 1.0     !Afe   !Iron abundance      
  param(10) = 0.2     !kTbb in observer's restframe
  param(11) = 300.0   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
  param(12) = 0.0     !Nh
  param(13) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
  param(14) = 4.6e7   !M     !BH mass in solar masses
  param(15) = 0.0     !flo   !Lowest frequency in band (Hz)
  param(16) = 0.0     !fhi   !Highest frequency in band (Hz)
  param(17) = 1       !ReIm  !1=Re, 2=Im, 3=Modulus, 4=phase lag (cycles), 5=time lag (s)
  param(18) = 0.0     !phiA
  param(19) = 0.0     !phiB
  param(20) = 0.0     !gamma
  !---------------------------------
  
  Emax  = 500.0
  Emin  = 0.1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

  call tdreltransCpT(ear,ne,param,ifl,photar)

  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     write(99,*)E,E**2*photar(i)/dE
  end do
  
end program main_test
