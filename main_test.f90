PROGRAM  MAIN
! Settings:
! nro,nphi  The observer's camera has nro*nphi pixels: nro radial and nphi azimuthal
! nr        The transfer function is calculated for nr radial binsfrom rmin to emax.
!           Transfer functions for different rin can then be calculated by interpolation
! nex       The number of (logarithmic) energy bins used for the transfer function
! nt        The number of (logarithmic) time bins used for the transfer function
! ndelta    The number of geodisics used for the emissivity calculation
! nup       The number of points along each geodesic used for ray tracing from observer to disk
! nup_lp    The number of points along each geodesic used for the emissivity calculation
! rinmax    The maximum value of rin considered. This is typically a lot smaller than rout.
      IMPLICIT NONE
      real param(19),numin,numax
      integer k,ne,i,ifl,kmax
      parameter (ne=1000)
      real ear(0:ne),emin,emax,t0,t1,photar(ne),E,dE
      real :: freq1(5),freq2(5),vector(5)
      character (len=200) name
      
      !----Parameters-------------------
      param(1)  = 6.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
      param(2)  = 0.9     !a     !BH spin
      param(3)  = 30.0    !inc   !Inclination angle in degrees
      param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
      param(5)  = 20000.0 !rout  !Disk outer radius in Rg - will probably hardwire this
      param(6)  = 0.0     !zcos  !Cosmological redshift
      param(7)  = 2.0     !Gamma !Photon index
      param(8)  = 3.0     !logxi !log10xi - ionisation parameter
      param(9)  = 1.0     !Afe   !Iron abundance      
      param(10) = 300.0   !Ecut  !High energy exponential cut off ***IN OBSERVER'S RESTFRAME***
      param(11) = 0.0     !Nh    !Hydrogen absorption column (using tbabs)
      param(12) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
      param(13) = 4.6e7   !M     !BH mass in solar masses
      param(14) = 1e-5    !flo   !Lowest frequency in band (Hz)
      param(15) = 2e-5    !fhi   !Highest frequency in band (Hz)
      param(16) = 1       !ReIm  !1=Re, 2=Im, 3=Modulus, 4=time lag (s), 5=Modulus (with response), 6=time lag (with response)
      param(17) = 0.0     !phiA  !Frequency-dependent phase normalisation (radians) - calculate self-consistently in full version of the model
      param(18) = 0.0     !phiB  !Frequency-dependent phase of the power-law index oscillation 
      param(19) = 0.0     !g     !Ratio of the two phases amplitude 
      !---------------------------------
      
      Emax  = 500.0 !300.0
      Emin  = 0.1   !1.0
      do i = 0,ne
        ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
      end do


      numin = 1e-8
      numax = 1e-4
      kmax  = 100


        ! call CPU_TIME(t0)
        ! call tdreltrans(ear,ne,param,ifl,photar)
        ! call CPU_TIME(t1)
        ! write(*,*)"Total CPU time=",t1-t0


!*********************TEST**************************
      
      param(1)  = 3.0
      param(3)  = 30.0
      param(13) = 14.7   !M     !BH mass in solar masses
      param(14) = 1    !flo   !Lowest frequency in band (Hz)
      param(15) = 2   !fhi   !Highest frequency in band (Hz)
      param(17) = 0.0     !phiA  
      param(18) = 0.2     !phiB  
      param(19) = 0.1     !g     

      kmax = 1000
      
      do k=1,kmax

         param(1) = param(1) + float(k)/(float(kmax)*10)
      write(*,*) "h = ",param(1)
      
        call CPU_TIME(t0)
        call tdreltrans(ear,ne,param,ifl,photar)
        call CPU_TIME(t1)
        write(*,*)"Total CPU time=",t1-t0

     enddo
!*******************************************
      

!Print the model in test_prova.dat
!*******************************************
 !       name = '../sim_data/200_Re_expar.dat'
 ! !     name = '../sim_data/050_lag.dat'
 !      open(99,file=name)
      
 !      write(99,*) 'skip on' 

 !      param(16) = 1

 !      param(1)  = 3.0
 !      param(3)  = 80.0
 !      param(13) = 10   !M     !BH mass in solar masses
 !      param(14) = 1    !flo   !Lowest frequency in band (Hz)
 !      param(15) = 2   !fhi   !Highest frequency in band (Hz)
 !      param(17) = 0.0     !phiA  
 !      param(18) = 0.2     !phiB  
 !      param(19) = 0.1     !g     
      
 !        call CPU_TIME(t0)
 !        call tdreltrans(ear,ne,param,ifl,photar)
 !        call CPU_TIME(t1)
 !        write(*,*)"Total CPU time=",t1-t0
        
 !        if( param(16) .lt. 4 )then
 !          do i = 1,ne
 !            E  = 0.5 * ( ear(i) + ear(i-1) )
 !            dE =         ear(i) - ear(i-1)
 !            write(99,*)E,E**2*photar(i)/dE
 !          end do
 !        else
 !          do i = 1,ne
 !            E  = 0.5 * ( ear(i) + ear(i-1) )
 !            dE =         ear(i) - ear(i-1)
 !            write(99,*)E,photar(i)/dE
 !          end do
 !        end if
          
        ! write(99,*)"log"
!*******************************************


!Print the model 5 times with different frequency values 
!*******************************************
!       name = '../sim_test/test_freq.dat'
!       open(99,file=name)
!       write(99,*) 'skip on'
      
! !Set the frequency ranges       
!       freq1(1) = 0.05
!       freq2(1) = 0.06 
!       freq1(2) = 0.5
!       freq2(2) = 0.6
!       freq1(3) = 1.0
!       freq2(3) = 2.0
!       freq1(4) = 10.0
!       freq2(4) = 11.0
!       freq1(5) = 30.0
!       freq2(5) = 31.0

! !changing some of the fiducial parametes to compare with Mastroserio et al 2018       
!       param(1)  = 10.0     !h     
!       param(2)  = 0.98     !a     
!       param(3)  = 45.0    !inc   
!       param(4)  = 10.0    !rin   
!       param(8)  = 3.1     !logxi 
!       param(13) = 10.0  !Set the mass for a black hole binary
!       param(16) = 4 !set the lag as output
!       param(17) = 0.0     !phiA  
!       param(18) = 0.2     !phiB  
!       param(19) = 0.1     !g     

!       do k = 1,5

!         param(14) = freq1(k)
!         param(15) = freq2(k)
         
!         call CPU_TIME(t0)
!         call tdreltrans(ear,ne,param,ifl,photar)
!         call CPU_TIME(t1)
!         write(*,*)"Total CPU time=",t1-t0

        
!         if( param(16) .lt. 4 )then
!           do i = 1,ne
!             E  = 0.5 * ( ear(i) + ear(i-1) )
!             dE =         ear(i) - ear(i-1)
!             write(99,*)E,E**2*photar(i)/dE
!           end do
!         else
!           do i = 1,ne
!             E  = 0.5 * ( ear(i) + ear(i-1) )
!             dE =         ear(i) - ear(i-1)
!             write(99,*)E,photar(i)/dE
!           end do
!         end if
          
!         write(99,*)"no no"
!       end do
!*******************************************


!Testing the time 
!*******************************************
      ! name = '../sim_test/test_prova1.dat'
      ! open(99,file=name)
      ! vector(1) = -1.0 
      ! vector(2) = -2.0
      ! vector(3) = -2.5
      ! vector(4) = -3.0
      ! vector(5) = -5.0

      ! param(13) = 10.0  !Set the mass for a black hole binary
      ! param(14) = 0.05    !flo   !Lowest frequency in band (Hz)
      ! param(15) = 0.06   !fhi   !Highest frequency in band (Hz)

      ! do k = 1,5

      ! param(4)  = vector(k)     !a 
      ! call CPU_TIME(t0)
      ! call tdreltrans(ear,ne,param,ifl,photar)
      ! call CPU_TIME(t1)
      ! write(*,*)"Total CPU time=",t1-t0
        
      !   if( param(16) .lt. 4 )then
      !     do i = 1,ne
      !       E  = 0.5 * ( ear(i) + ear(i-1) )
      !       dE =         ear(i) - ear(i-1)
      !       write(99,*)E,E**2*photar(i)/dE
      !     end do
      !   else
      !     do i = 1,ne
      !       E  = 0.5 * ( ear(i) + ear(i-1) )
      !       dE =         ear(i) - ear(i-1)
      !       write(99,*)E,photar(i)/dE
      !     end do
      !   end if
          
      !   write(99,*)"no no"
      ! end do
!*******************************************


      
      
      end program main
