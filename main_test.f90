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
      real ear(0:ne),emin,emax,t0,t1,photar(ne),E,dE,cont(ne)
      real :: freq1(5),freq2(5),vector(5),av
      character (len=200) name
      character (len=100) char
      double precision dear(0:ne),drelxillpar(12),dphotar(ne),dphoter(ne)
      double precision dcont(ne)
      
      !----Parameters-------------------
      param(1)  = 10.0     !h     !Source height **-ve means in units of BH horizon, +ve means in Rg***
      param(2)  = 0.9     !a     !BH spin
      param(3)  = 30.0    !inc   !Inclination angle in degrees
      param(4)  = -1.0    !rin   !Disk inner radius **-ve means in units of ISCO, +ve means in Rg***
      param(5)  = 400.0   !rout  !Disk outer radius in Rg - will probably hardwire this
      param(6)  = 0.0     !zcos  !Cosmological redshift
      param(7)  = 2.0     !Gamma !Photon index
      param(8)  = 3.0     !logxi !log10xi - ionisation parameter
      param(9)  = 1.0     !Afe   !Iron abundance      
      param(10) = 300.0   !Ecut  !High energy exponential cut off ***IN OBSERVER'S RESTFRAME***
      param(11) = 0.0     !Nh    !Hydrogen absorption column (using tbabs)
      param(12) = 1.0     !1onB  !(1/\mathcal{B}): boosting fudge factor that lowers normalisation of reflection spectrum
      param(13) = 4.6e7   !M     !BH mass in solar masses
      param(14) = 0.0     !flo   !Lowest frequency in band (Hz)
      param(15) = 0.0     !fhi   !Highest frequency in band (Hz)
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
      kmax  = 30

      do k = 1,1 != 1,kmax

!        param(1) = 6.0 + real(k-1) * 94.0 / real(kmax-1)
!         param(3) = 10.0 + real(k-1) * 75.0 / real(kmax-1)
!         write(*,*)"param(3)=",param(3)

!         param(2) = real(k-1) * 0.998 / real(kmax-1)
!         param(8) = 2.0 + real(k-1) * 1.5 / real(kmax-1)
         
      !Write out full model of reltrans
      param(12) = 1.0
        
      call CPU_TIME(t0)
      call tdreltrans(ear,ne,param,ifl,photar)
     write(99,*) 'skip on'
     write(99,*) 'read serr 1 2'
      
      if( param(16) .lt. 4 )then
        do i = 1,ne
          E  = 0.5 * ( ear(i) + ear(i-1) )
          dE =         ear(i) - ear(i-1)
          write(99,*)E,0.5*dE,E**2*photar(i)/dE
        end do
      else
        do i = 1,ne
          E  = 0.5 * ( ear(i) + ear(i-1) )
          dE =         ear(i) - ear(i-1)
          write(99,*)E,0.5*dE,photar(i)/dE
        end do
      end if
      write(99,*) 'no no'
      call CPU_TIME(t1)
      write(*,*)"Total CPU time=",t1-t0

      !Now just continuum

      param(12) = 0.0
      call tdreltrans(ear,ne,param,ifl,cont)
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,E**2*cont(i)/dE
        photar(i) = photar(i) - cont(i)
      end do
      write(99,*)"no no"
      
      !Now write out just reflection
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,E**2*photar(i)/dE
      end do
      write(99,*)"no no"
      
!Call relxilllp with the same model parameters

      !Energy grid
      dear = dble( ear )

      !Same parameters as reltrans
      do i = 1,10
        drelxillpar(i) = dble( param(i) )
      end do
      drelxillpar(11) = 0.0
      drelxillpar(12) = 1
      
      call lmodrelxilllpf(dear,ne,drelxillpar,ifl,dphotar,dphoter,char)

      !Write out total
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,E**2*real(dphotar(i))/dE
      end do
      write(99,*)"no no"
      
      !Now just the continuum
      drelxillpar(11) = 0.0
      drelxillpar(12) = 0
      call lmodrelxilllpf(dear,ne,drelxillpar,ifl,dcont,dphoter,char)
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,E**2*real(dcont(i))/dE
        dphotar(i) = dphotar(i) - dcont(i)
      end do
      write(99,*)"no no"
      
      !Finally just reflection
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,E**2*real(dphotar(i))/dE
      end do
      write(99,*)"no no"

! Write out the ratio between the two
      av = 0.0
      do i = 1,ne
        E  = 0.5 * ( ear(i) + ear(i-1) )
        dE =         ear(i) - ear(i-1)
        write(99,*)E,0.5*dE,real(dphotar(i))/photar(i)
        av = av + real(dphotar(i))/photar(i)
      end do
      write(99,*)"no no"
      
      write(99,*)"skip on"
      write(99,*)"log"
      write(99,*)"li s on"
      write(99,*)"lw 5"
      write(99,*)"co 1 on 1,2,3"
      write(99,*)"co 2 on 4,5,6"
      write(99,*)"co off 1,4"
      write(99,*)"co 15 on 7"
      
      av = av / real(ne)
      write(*,*)"Average ratio=",av
      write(64,*)param(2),av

      end do
      
      end program main
