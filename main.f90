program relwrap
    implicit none

    real    :: params_file(32)      !fulll list of all model parameters, indipendent of model flavour, from the input file
    real    :: params_reltrans(27)  !stuff to be used from the parameters in reltransDbl
    real    :: emin, emax           !minimum, maximum energy and increment to set model grid
    integer :: ne                   !energy grid resolution
    integer :: ifl                  !weird integer thing that equals 1 for some reason I don't understand. Just pass it to wrapper and
                                    !don't ask questions
    integer :: i , j                !loop variables 
    real time_start,time_end        !runtime stuff
                                    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array
    
    logical single_lp               !switch to change between single and double LP models
    
    single_lp = .true. 

    !set up energy and countrate arrays  
    ne = 1000                             
    allocate(ear(0:ne))  
    allocate(photar(ne))       
    emin = 0.1
    emax = 200.
    do i=0,ne
        ear(i) = emin * (emax/emin)**(real(i)/real(ne))
    end do                                             
    
    print *,"Reltrans wrapper code."
    print *,"--------------------------------------------------------------------------------------------------------------"
    
    open(1,file="Input/ip.dat",status='old')
    read(1,*) params_file
    close(1)
    
    if( single_lp .eqv. .true. ) then
        !set up parameters for single LP model
        params_reltrans(1) = params_file(1)
        params_reltrans(2) = params_file(3)
        params_reltrans(3) = params_file(4)
        params_reltrans(4) = params_file(5)
        params_reltrans(5) = params_file(6)
        params_reltrans(6) = params_file(7)
        params_reltrans(7) = params_file(8)
        params_reltrans(8) = params_file(9)
        params_reltrans(9) = params_file(10)
        params_reltrans(10) = params_file(11)
        params_reltrans(11) = params_file(12)
        params_reltrans(12) = params_file(16)
        params_reltrans(13) = params_file(17)
        params_reltrans(14) = params_file(18)
        params_reltrans(15) = params_file(19)
        params_reltrans(16) = params_file(23)
        params_reltrans(17) = params_file(24)
        params_reltrans(18) = params_file(25)
        params_reltrans(19) = params_file(26)
        params_reltrans(20) = params_file(27)
        params_reltrans(21) = params_file(28)
        params_reltrans(22) = params_file(32)
        
        print *, "Single LP model: "
        print *, "Height: ", params_reltrans(1)
        print *, "Spin: ", params_reltrans(2)
        print *, "Inclination: ", params_reltrans(3)
        print *, "Rin: ", params_reltrans(4)
        print *, "Rout: ", params_reltrans(5)
        print *, "redshift: ", params_reltrans(6)
        print *, "Gamma: ", params_reltrans(7)
        print *, "Log Xi: ", params_reltrans(8)
        print *, "Afe: ", params_reltrans(9)
        print *, "Log ne: ", params_reltrans(10)
        print *, "Ecut/Te: ", params_reltrans(11)
        print *, "Nh: ", params_reltrans(12)
        print *, "boost: ", params_reltrans(13)
        print *, "diffusion time: ",params_reltrans(14)
        print *, "Mass: ", params_reltrans(15)
        print *, "flo: ", params_reltrans(16)
        print *, "fhi: ", params_reltrans(17)
        print *, "ReIm: ", params_reltrans(18)
        print *, "DelA: ", params_reltrans(19)
        print *, "DelAB: ", params_reltrans(20)
        print *, "g: ", params_reltrans(21)   
        print *, "# of rsp:", params_reltrans(22)  
    else
        !set up parameters for double LP model
        params_reltrans(1) = params_file(1)
        params_reltrans(2) = params_file(2)
        params_reltrans(3) = params_file(3)
        params_reltrans(4) = params_file(4)
        params_reltrans(5) = params_file(5)
        params_reltrans(6) = params_file(6)
        params_reltrans(7) = params_file(7)
        params_reltrans(8) = params_file(8)
        params_reltrans(9) = params_file(9)
        params_reltrans(10) = params_file(10)
        params_reltrans(11) = params_file(11)
        params_reltrans(12) = params_file(12)
        params_reltrans(13) = params_file(13)
        params_reltrans(14) = params_file(14)
        params_reltrans(15) = params_file(15)
        params_reltrans(16) = params_file(16)
        params_reltrans(17) = params_file(17)
        params_reltrans(18) = params_file(19)
        params_reltrans(19) = params_file(23)
        params_reltrans(20) = params_file(24)
        params_reltrans(21) = params_file(25)
        params_reltrans(22) = params_file(26)
        params_reltrans(23) = params_file(27)
        params_reltrans(24) = params_file(28)
        params_reltrans(25) = params_file(29)
        params_reltrans(26) = params_file(30)
        params_reltrans(27) = params_file(32)
        
        print *, "Double LP model: "      
        print *, "Height 1: ", params_reltrans(1)
        print *, "Height 2: ", params_reltrans(2)
        print *, "Spin: ", params_reltrans(3)
        print *, "Inclination: ", params_reltrans(4)
        print *, "Rin: ", params_reltrans(5)
        print *, "Rout: ", params_reltrans(6)
        print *, "redshift: ", params_reltrans(7)
        print *, "Gamma: ", params_reltrans(8)
        print *, "Log Xi: ", params_reltrans(9)
        print *, "Afe: ", params_reltrans(10)
        print *, "Log ne: ", params_reltrans(11)
        print *, "Ecut/Te: ", params_reltrans(12)
        print *, "eta_0:", params_reltrans(13)
        print *, "eta:", params_reltrans(14)
        print*,  "beta_p:", params_reltrans(15)
        print *, "Nh: ", params_reltrans(16)
        print *, "boost: ", params_reltrans(17)
        print *, "Mass: ", params_reltrans(18)
        print *, "flo: ", params_reltrans(19)
        print *, "fhi: ", params_reltrans(20)
        print *, "ReIm: ", params_reltrans(21)
        print *, "DelA: ", params_reltrans(22)
        print *, "DelAB1: ", params_reltrans(23)
        print *, "g1: ", params_reltrans(24) 
        print *, "DelAB2: ", params_reltrans(25)
        print *, "g2: ", params_reltrans(26)   
        print *, "# of rsp:", params_reltrans(27)         
    end if
        
    !set up frequency array appropriately if using lag/frequency spectra, depending on BH mass:
    if( params_file(25) .eq. 7 ) then
        if( params_file(19) .gt. 1.e4 ) then
            emin = 1.e-4
            emax = 3.e-2
            do i=0,ne
                ear(i) = emin * (emax/emin)**(real(i)/real(ne))
            end do    
        else
            emin = 0.1
            emax = 700.
            do i=0,ne
                ear(i) = emin * (emax/emin)**(real(i)/real(ne))
            end do   
        end if
    end if
      
    call CPU_TIME (time_start)
    if( single_lp .eqv. .true. ) then
        call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    else
        call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)    
    end if    

    call CPU_TIME (time_end)
    print *, 'Call runtime: ', time_end - time_start, ' seconds'
    print *,"--------------------------------------------------------------------------------------------------------------"
            
    !uncomment this to plot the output and all the components
    call system("python3 Plot.py")     
                
    deallocate(ear)
    deallocate(photar)
    
end program
