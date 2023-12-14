program main_simple_call
 implicit none
 integer ne,i,ifl,neb
 parameter (ne=5000,neb=25)
 real Emax,Emin,ear(0:ne),param(21),photar(ne),E,dE
 real dist_par(24),spar(26),earb(0:neb)
 
! Set energy grid
  Emax  = 500.0
  Emin  = 0.1
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

! !----Parameters-------------------
  param(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
  param(2)  = 0.998     !a     !BH spin
  param(3)  = 30.0    !inc   !Inclination angle in degrees
  param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
  param(5)  = 1e3 !rout  !Disk outer radius in Rg - will probably hardwire this
  param(6)  = 0.0     !zcos  !Cosmological redshift
  param(7)  = 2.0     !Gamma !Photon index
  param(8)  = 1.0    !logxi !log10xi - ionisation parameter
  param(9)  = 1.0     !Afe   !Iron abundance      
  param(10) = 15.0   !density
  param(11) = 60.0   !kTe   !Electron temperature ***IN OBSERVER'S RESTFRAME***
  param(12) = 0.0     !Nh
  param(13) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
  param(14) = 10.0    !M     !BH mass in solar masses
  param(15) = 0.0    !flo   !Lowest frequency in band (Hz)
  param(16) = 0.0   !fhi   !Highest frequency in band (Hz)
  param(17) = 1       !ReIm  !1=Re, 2=Im, 3=modulus, 4=time lag (s), 5=folded modulus, 6=folded time lag (s)
  param(18) = 0.0     !DelA
  param(19) = 0.0     !DelAB
  param(20) = 0.0     !gamma
  param(21) = 1       !telescope response= 
! !---------------------------------

! !----Parameters-------------------
!   param(1)  = 29.7014   
!   param(2)  = 0.5       
!   param(3)  = 44.2792   
!   param(4)  = -1        
!   param(5)  = 1000      
!   param(6)  = 0         
!   param(7)  = 1.79279   
!   param(8)  = 1.44735   
!   param(9)  = 1         
!   param(10) = 19.9874   
!   param(11) = 60        
!   param(12) = 0.0       
!   param(13) = 1         
!   param(14) = 10        
!   param(15) = 0.0         
!   param(16) = 0.0         
!   param(17) = 1         
!   param(18) = 0.0
!   param(19) = 0.0
!   param(20) = 0.0   
!   param(21) = 1         
! !---------------------------------

  ! Call one of the existing models

  if (param(13) .gt. 0.0) then
     write(99,*)"skip on"
     call tdreltransDCp(ear,ne,param,ifl,photar)
     do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE = ear(i) - ear(i-1)
        if( param(16) .gt. 3.5 .and. param(16) .lt. 4.5 .or. param(16) .gt. 5.5 )then
           write(99,*)E,photar(i)/dE
        else
           write(99,*)E,E**2*photar(i)/dE
        end if
     end do
     write(99,*)"no no"

  else
     write(99,*)"skip on"
     call tdreltransDCp(ear,ne,param,ifl,photar)
     do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE = ear(i) - ear(i-1)
        if( param(16) .gt. 3.5 .and. param(16) .lt. 4.5 .or. param(16) .gt. 5.5 )then
           write(88,*)E,photar(i)/dE
        else
           write(88,*)E,E**2*photar(i)/dE
        end if
     end do
     write(99,*)"no no"


     ! call xsatbl(ear, ne, param_xil, trim(pathname_xillver), ifl, photar, photer)


     
  endif
     
!   write(*,*)"---------------------------------------------------------"
!   write(*,*)"Outputs:"
!   write(*,*)"fort.188: r,logxi,logxieff"
!   write(*,*)"fort.99: spectrum / lag spectrum etc"
!   write(*,*)"---------------------------------------------------------"

  
end program main_simple_call
