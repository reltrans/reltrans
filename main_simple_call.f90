program relwrap
    implicit none

    real    :: params_file(21)      !thing to read in from files
    real    :: emin, emax           !minimum, maximum energy and increment to set model grid
    integer :: ne                   !energy grid resolution
    integer :: ifl                  !weird integer thing that equals 1 for some reason I don't understand. Just pass it to wrapper and
                                    !don't ask questions
    integer :: i                    !loop variable 
    real time_start,time_end        !runtime stuff
    real test                       !check to see if Fortran is silly
                                    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array? maybe? 
    
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
    
    print *, "Height: ", params_file(1)
    print *, "Spin: ", params_file(2)
    print *, "Inclination: ", params_file(3)
    print *, "Rin: ", params_file(4)
    print *, "Rout: ", params_file(5)
    print *, "redshift: ", params_file(6)
    print *, "Gamma: ", params_file(7)
    print *, "Log Xi: ", params_file(8)
    print *, "Afe: ", params_file(9)
    print *, "Log ne: ", params_file(10)
    print *, "Ecut/Te: ", params_file(11)
    print *, "Nh: ", params_file(12)
    print *, "boost: ", params_file(13)
    print *, "Mass: ", params_file(14)
    print *, "flo: ", params_file(15)
    print *, "fhi: ", params_file(16)
    print *, "ReIm: ", params_file(17)
    print *, "DelA: ", params_file(18)
    print *, "DelAB: ", params_file(19)
    print *, "g: ", params_file(20)   
    print *, "# of rsp:", params_file(21) 
    
    call CPU_TIME (time_start)
    call tdreltransDCp(ear,ne,params_file,ifl,photar)
    call CPU_TIME (time_end)
    print *, 'Call runtime: ', time_end - time_start, ' seconds'
    print *,"--------------------------------------------------------------------------------------------------------------"
    
    !uncomment this to plot the output and all the components
    call system("python3 Plot.py") 
             
    deallocate(ear)
    deallocate(photar)
    
end program
