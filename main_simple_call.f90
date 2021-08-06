program relwrap
    implicit none

    real    :: params_file(25)      !fulll list of all model parameters, indipendent of model flavour
    real    :: params_DCp(21)       !stuff to be used from the parameters in reltransDCp
    real    :: emin, emax           !minimum, maximum energy and increment to set model grid
    integer :: ne                   !energy grid resolution
    integer :: ifl                  !weird integer thing that equals 1 for some reason I don't understand. Just pass it to wrapper and
                                    !don't ask questions
    integer :: i                    !loop variable 
    real time_start,time_end        !runtime stuff
                                    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array
    
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
    
    params_DCp(1) = params_file(1)
    params_DCp(2) = params_file(2)
    params_DCp(3) = params_file(3)
    params_DCp(4) = params_file(4)
    params_DCp(5) = params_file(5)
    params_DCp(6) = params_file(6)
    params_DCp(7) = params_file(7)
    params_DCp(8) = params_file(8)
    params_DCp(9) = params_file(9)
    params_DCp(10) = params_file(10)
    params_DCp(11) = params_file(11)
    params_DCp(12) = params_file(12)
    params_DCp(13) = params_file(13)
    params_DCp(14) = params_file(14)
    params_DCp(15) = params_file(18)
    params_DCp(16) = params_file(19)
    params_DCp(17) = params_file(20)
    params_DCp(18) = params_file(21)
    params_DCp(19) = params_file(22)
    params_DCp(20) = params_file(23)
    params_DCp(21) = params_file(25)
    
    print *, "Height: ", params_DCp(1)
    print *, "Spin: ", params_DCp(2)
    print *, "Inclination: ", params_DCp(3)
    print *, "Rin: ", params_DCp(4)
    print *, "Rout: ", params_DCp(5)
    print *, "redshift: ", params_DCp(6)
    print *, "Gamma: ", params_DCp(7)
    print *, "Log Xi: ", params_DCp(8)
    print *, "Afe: ", params_DCp(9)
    print *, "Log ne: ", params_DCp(10)
    print *, "Ecut/Te: ", params_DCp(11)
    print *, "Nh: ", params_DCp(12)
    print *, "boost: ", params_DCp(13)
    print *, "Mass: ", params_DCp(14)
    print *, "flo: ", params_DCp(15)
    print *, "fhi: ", params_DCp(16)
    print *, "ReIm: ", params_DCp(17)
    print *, "DelA: ", params_DCp(18)
    print *, "DelAB: ", params_DCp(19)
    print *, "g: ", params_DCp(20)   
    print *, "# of rsp:", params_DCp(21) 
    
    call CPU_TIME (time_start)
    call tdreltransDCp(ear,ne,params_DCp,ifl,photar)
    call CPU_TIME (time_end)
    print *, 'Call runtime: ', time_end - time_start, ' seconds'
    print *,"--------------------------------------------------------------------------------------------------------------"
    
    !uncomment this to plot the output and all the components
    call system("python3 Plot.py") 
             
    deallocate(ear)
    deallocate(photar)
    
end program
