!-----------------------------------------------------------------------
subroutine initialiser(firstcall,Emin,Emax,nex,dloge,earx,rnmax,d,needtrans,check&
     ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim)
! Initialises the model and writes the header
      implicit none
      integer nex,i,spin_dim,mu_dim,nphi,nro,nphi_grid,nro_grid,lrec,irec,check
      logical firstcall,needtrans
      real Emin,Emax,dloge,earx(0:nex)
      double precision :: honr_grid,rout_grid,d_grid
      double precision :: d,rnmax,spin_start,spin_end,mu_start,mu_end
      character (len=500) gridname
      needtrans = .false.
      call get_environment_variable("GRID",gridname,check)

      if( firstcall )then
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS: a transfer function model for"
        write(*,*)"X-ray reverberation mapping written by Adam Ingram."
        write(*,*)"Please cite Ingram et al (2018)."
        write(*,*)"----------------------------------------------------"
        !Create *logarithmic* working energy grid
        !Will need to evaluate xillver on this grid to use the FT convolution code 
        Emax  = 1e3
        Emin  = 1e-1
        dloge = log10( Emax / Emin ) / float(nex)
        do i = 0,nex
          earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
        end do


        if (check .ne. 0) then 
           write(*,*) 'This code uses a grid to compute the kernel trasfer funcion' 
      
!open the grid 
           lrec = 8*nphi*nro         
           open(98,file=gridname,access='direct',form='unformatted',status='old',recl=lrec)
           irec = 1 
           read(98,rec=1) rnmax,nphi_grid,nro_grid,honr_grid,rout_grid,d_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim
!remember that rout_grid is the rout used to make the grid. It is not the same as param 5
!         write(*,*) 'rout of the grid', rout_grid


           
!check if the grid has the correct values 
           if (nphi_grid .ne. nphi .or. nro_grid .ne. nro) then
              write(*,*) 'Not compatible grid dimentions'
              stop
           endif

           
! Set sensible distance for observer from the BH now that we took rnmax from the grid 
           d = max( 1.0d4 , 2.0d2 * rnmax**2 )
         
!check if the grid distance has been calculated in the same way
           if (d_grid .ne. d ) then
              write(*,*) 'The distance has been computed differently in the grid'
              stop
           endif
        else
           rnmax= 300.d0
! Set sensible distance for observer from the BH
           d = max( 1.0d4 , 2.0d2 * rnmax**2 )
           
        endif

        firstcall = .false.
     end if
     return
    end subroutine initialiser
!-----------------------------------------------------------------------
