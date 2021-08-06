program test_reltrans

    !this program runs reltrans using the same input parameters as Mastroserio et al. 2021, and checks whether the total model,
    !as well as the contribution of individual transfer functions, is identical to that paper. Note that in order to run this,
    !you NEED to set the environment variable VERBOSE to 2 and the call to FININT uncommentated

    real    :: params_file(25)      !fulll list of all model parameters, indipendent of model flavour
    real    :: params_reltrans(21)  !stuff to be used from the parameters in reltransDCp; change this is if you want to try 
                                    !a different model flavour
    real    :: emin, emax           !minimum, maximum energy and increment to set model grid
    integer :: ne                   !energy grid resolution
    integer :: ifl                  !weird integer thing that equals 1 for some reason I don't understand. Just pass it to wrapper and
                                    !don't ask questions
    integer :: i                    !loop variable 
    real time_start,time_end        !runtime stuff
    real test_model                 !temporary value to test the model output
    real test_precision             !accuracy to which we want the tests to check the model; by default test_precision, ie 0.1% change will
                                    !be flagged as a failed test
    logical test_bool               !boolean to check whether tests are passed or not
    logical full_print              !option to print extra information beyond the results of the individual tests
                                    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array
    
    real    :: total_benchmark(1000,2)      !total model from the benchmark file
    real    :: pivoting_benchmark(1000,2)   !continuum lag from the benchmark file
    real    :: refl_benchmark(1000,2)       !pivoting reflection from the benchmark file
    real    :: reverb_benchmark(1000,2)     !reverberation from the benchmark file
    real    :: total_model(1000,2)      !total model from the model file
    real    :: pivoting_model(1000,2)   !continuum lag from the model file
    real    :: refl_model(1000,2)       !pivoting reflection from the model file
    real    :: reverb_model(1000,2)     !reverberation from the model file
    
    full_print = .false.
    test_precision = 1.001
    
    !set up energy and countrate arrays. Note that the grid HAS to be 1000 in lenght and be defined between 0.1 and 200 keV for
    !the comparison with the benchmark lag-energy spectra to make sense
    ne = 1000                             
    allocate(ear(0:ne))  
    allocate(photar(ne))       
    emin = 0.1
    emax = 200.
    do i=0,ne
        ear(i) = emin * (emax/emin)**(real(i)/real(ne))
    end do     
    
    print *,"----------------------------------------------------------------"
    print *,"Reltrans test; 0.12-0.25 Hz"  
    !get parameter file fron the benchmark folder for the first frequency interval
    open(1,file="Benchmarks/ip_0,12_0,25.dat",status='old')
    read(1,*) params_file
    close(1)    
    
    params_reltrans(1)  = params_file(1)
    params_reltrans(2)  = params_file(2)
    params_reltrans(3)  = params_file(3)
    params_reltrans(4)  = params_file(4)
    params_reltrans(5)  = params_file(5)
    params_reltrans(6)  = params_file(6)
    params_reltrans(7)  = params_file(7)
    params_reltrans(8)  = params_file(8)
    params_reltrans(9)  = params_file(9)
    params_reltrans(10) = params_file(10)
    params_reltrans(11) = params_file(11)
    params_reltrans(12) = params_file(12)
    params_reltrans(13) = params_file(13)
    params_reltrans(14) = params_file(14)
    params_reltrans(15) = params_file(18)
    params_reltrans(16) = params_file(19)
    params_reltrans(17) = params_file(20)
    params_reltrans(18) = params_file(21)
    params_reltrans(19) = params_file(22)
    params_reltrans(20) = params_file(23)
    params_reltrans(21) = params_file(25)   
    
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !read benchmark files and output files 
    open(20,file="Benchmarks/Total_0,12_0,25.dat",status='old')
    open(21,file="Benchmarks/PivotingPL_0,12_0,25.dat",status='old')
    open(22,file="Benchmarks/PivotingReflection_0,12_0,25.dat",status='old')
    open(23,file="Benchmarks/LightTravelTime_0,12_0,25.dat",status='old')
    open(30,file="Output/Total.dat",status='old')
    open(31,file="Output/PivotingPL.dat",status='old')
    open(32,file="Output/PivotingReflection.dat",status='old')
    open(33,file="Output/LightTravelTime.dat",status='old')
    do i=1,ne
        read(20,*) (total_benchmark(i,j), j=1, 2)
        read(21,*) (pivoting_benchmark(i,j), j=1, 2)
        read(22,*) (refl_benchmark(i,j), j=1, 2)
        read(23,*) (reverb_benchmark(i,j), j=1, 2)        
        read(30,*) (total_model(i,j), j=1, 2)
        read(31,*) (pivoting_model(i,j), j=1, 2)
        read(32,*) (refl_model(i,j), j=1, 2)
        read(33,*) (reverb_model(i,j), j=1, 2)                           
    end do
    close(20) 
    close(21) 
    close(22) 
    close(23)     
    close(30) 
    close(31) 
    close(32) 
    close(33) 
   
    !compare total and individual components
    print *,"----------------------------------------------------------------"
    print *, "Comparing total model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = total_benchmark(i,2)/total_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Total model output different at energy bin", total_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", total_benchmark(i,2), total_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Total model output test passed"
    else                                    
        print*, "Total model output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting continuum output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = pivoting_benchmark(i,2)/pivoting_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting continuum output different at energy bin", pivoting_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", pivoting_benchmark(i,2), pivoting_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting continuum output test passed"
    else                                    
        print*, "Pivoting continuum output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting reflection output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = refl_benchmark(i,2)/refl_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting reflection output different at energy bin", refl_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", refl_benchmark(i,2), refl_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting reflection output test passed"
    else                                    
        print*, "Pivoting reflection output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing reverberation output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = reverb_benchmark(i,2)/reverb_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Reverberation output different at energy bin", reverb_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", reverb_benchmark(i,2), reverb_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Reverberation output test passed"
    else                                    
        print*, "Reverberation output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *,"----------------------------------------------------------------"
    print *,"Reltrans test; 0.31-0.73 Hz"  
    !get parameter file fron the benchmark folder for the first frequency interval
    open(1,file="Benchmarks/ip_0,31_0,73.dat",status='old')
    read(1,*) params_file
    close(1)    
    
    params_reltrans(1)  = params_file(1)
    params_reltrans(2)  = params_file(2)
    params_reltrans(3)  = params_file(3)
    params_reltrans(4)  = params_file(4)
    params_reltrans(5)  = params_file(5)
    params_reltrans(6)  = params_file(6)
    params_reltrans(7)  = params_file(7)
    params_reltrans(8)  = params_file(8)
    params_reltrans(9)  = params_file(9)
    params_reltrans(10) = params_file(10)
    params_reltrans(11) = params_file(11)
    params_reltrans(12) = params_file(12)
    params_reltrans(13) = params_file(13)
    params_reltrans(14) = params_file(14)
    params_reltrans(15) = params_file(18)
    params_reltrans(16) = params_file(19)
    params_reltrans(17) = params_file(20)
    params_reltrans(18) = params_file(21)
    params_reltrans(19) = params_file(22)
    params_reltrans(20) = params_file(23)
    params_reltrans(21) = params_file(25)   
    
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !read benchmark files and output files 
    open(20,file="Benchmarks/Total_0,31_0,73.dat",status='old')
    open(21,file="Benchmarks/PivotingPL_0,31_0,73.dat",status='old')
    open(22,file="Benchmarks/PivotingReflection_0,31_0,73.dat",status='old')
    open(23,file="Benchmarks/LightTravelTime_0,31_0,73.dat",status='old')
    open(30,file="Output/Total.dat",status='old')
    open(31,file="Output/PivotingPL.dat",status='old')
    open(32,file="Output/PivotingReflection.dat",status='old')
    open(33,file="Output/LightTravelTime.dat",status='old')
    do i=1,ne
        read(20,*) (total_benchmark(i,j), j=1, 2)
        read(21,*) (pivoting_benchmark(i,j), j=1, 2)
        read(22,*) (refl_benchmark(i,j), j=1, 2)
        read(23,*) (reverb_benchmark(i,j), j=1, 2)        
        read(30,*) (total_model(i,j), j=1, 2)
        read(31,*) (pivoting_model(i,j), j=1, 2)
        read(32,*) (refl_model(i,j), j=1, 2)
        read(33,*) (reverb_model(i,j), j=1, 2)                           
    end do
    close(20) 
    close(21) 
    close(22) 
    close(23)     
    close(30) 
    close(31) 
    close(32) 
    close(33) 
   
    !compare total and individual components
    print *,"----------------------------------------------------------------"
    print *, "Comparing total model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = total_benchmark(i,2)/total_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Total model output different at energy bin", total_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", total_benchmark(i,2), total_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Total model output test passed"
    else                                    
        print*, "Total model output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting continuum output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = pivoting_benchmark(i,2)/pivoting_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting continuum output different at energy bin", pivoting_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", pivoting_benchmark(i,2), pivoting_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting continuum output test passed"
    else                                    
        print*, "Pivoting continuum output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting reflection output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = refl_benchmark(i,2)/refl_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting reflection output different at energy bin", refl_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", refl_benchmark(i,2), refl_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting reflection output test passed"
    else                                    
        print*, "Pivoting reflection output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing reverberation output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = reverb_benchmark(i,2)/reverb_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Reverberation output different at energy bin", reverb_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", reverb_benchmark(i,2), reverb_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Reverberation output test passed"
    else                                    
        print*, "Reverberation output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *,"----------------------------------------------------------------"
    print *,"Reltrans test; 0.80-2.1 Hz"  
    !get parameter file fron the benchmark folder for the first frequency interval
    open(1,file="Benchmarks/ip_0,80_2,1.dat",status='old')
    read(1,*) params_file
    close(1)    
    
    params_reltrans(1)  = params_file(1)
    params_reltrans(2)  = params_file(2)
    params_reltrans(3)  = params_file(3)
    params_reltrans(4)  = params_file(4)
    params_reltrans(5)  = params_file(5)
    params_reltrans(6)  = params_file(6)
    params_reltrans(7)  = params_file(7)
    params_reltrans(8)  = params_file(8)
    params_reltrans(9)  = params_file(9)
    params_reltrans(10) = params_file(10)
    params_reltrans(11) = params_file(11)
    params_reltrans(12) = params_file(12)
    params_reltrans(13) = params_file(13)
    params_reltrans(14) = params_file(14)
    params_reltrans(15) = params_file(18)
    params_reltrans(16) = params_file(19)
    params_reltrans(17) = params_file(20)
    params_reltrans(18) = params_file(21)
    params_reltrans(19) = params_file(22)
    params_reltrans(20) = params_file(23)
    params_reltrans(21) = params_file(25)   
    
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !read benchmark files and output files 
    open(20,file="Benchmarks/Total_0,80_2,1.dat",status='old')
    open(21,file="Benchmarks/PivotingPL_0,80_2,1.dat",status='old')
    open(22,file="Benchmarks/PivotingReflection_0,80_2,1.dat",status='old')
    open(23,file="Benchmarks/LightTravelTime_0,80_2,1.dat",status='old')
    open(30,file="Output/Total.dat",status='old')
    open(31,file="Output/PivotingPL.dat",status='old')
    open(32,file="Output/PivotingReflection.dat",status='old')
    open(33,file="Output/LightTravelTime.dat",status='old')
    do i=1,ne
        read(20,*) (total_benchmark(i,j), j=1, 2)
        read(21,*) (pivoting_benchmark(i,j), j=1, 2)
        read(22,*) (refl_benchmark(i,j), j=1, 2)
        read(23,*) (reverb_benchmark(i,j), j=1, 2)        
        read(30,*) (total_model(i,j), j=1, 2)
        read(31,*) (pivoting_model(i,j), j=1, 2)
        read(32,*) (refl_model(i,j), j=1, 2)
        read(33,*) (reverb_model(i,j), j=1, 2)                           
    end do
    close(20) 
    close(21) 
    close(22) 
    close(23)     
    close(30) 
    close(31) 
    close(32) 
    close(33) 
   
    !compare total and individual components
    print *,"----------------------------------------------------------------"
    print *, "Comparing total model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = total_benchmark(i,2)/total_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Total model output different at energy bin", total_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", total_benchmark(i,2), total_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Total model output test passed"
    else                                    
        print*, "Total model output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting continuum output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = pivoting_benchmark(i,2)/pivoting_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting continuum output different at energy bin", pivoting_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", pivoting_benchmark(i,2), pivoting_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting continuum output test passed"
    else                                    
        print*, "Pivoting continuum output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting reflection output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = refl_benchmark(i,2)/refl_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting reflection output different at energy bin", refl_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", refl_benchmark(i,2), refl_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting reflection output test passed"
    else                                    
        print*, "Pivoting reflection output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing reverberation output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = reverb_benchmark(i,2)/reverb_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Reverberation output different at energy bin", reverb_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", reverb_benchmark(i,2), reverb_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Reverberation output test passed"
    else                                    
        print*, "Reverberation output test failed, recommend re-running the test with full_print set to true!"
    endif    
     
     
    print *,"----------------------------------------------------------------"
    print *,"Reltrans test; 2.1-5.8 Hz"  
    !get parameter file fron the benchmark folder for the first frequency interval
    open(1,file="Benchmarks/ip_2,1_5,8.dat",status='old')
    read(1,*) params_file
    close(1)    
    
    params_reltrans(1)  = params_file(1)
    params_reltrans(2)  = params_file(2)
    params_reltrans(3)  = params_file(3)
    params_reltrans(4)  = params_file(4)
    params_reltrans(5)  = params_file(5)
    params_reltrans(6)  = params_file(6)
    params_reltrans(7)  = params_file(7)
    params_reltrans(8)  = params_file(8)
    params_reltrans(9)  = params_file(9)
    params_reltrans(10) = params_file(10)
    params_reltrans(11) = params_file(11)
    params_reltrans(12) = params_file(12)
    params_reltrans(13) = params_file(13)
    params_reltrans(14) = params_file(14)
    params_reltrans(15) = params_file(18)
    params_reltrans(16) = params_file(19)
    params_reltrans(17) = params_file(20)
    params_reltrans(18) = params_file(21)
    params_reltrans(19) = params_file(22)
    params_reltrans(20) = params_file(23)
    params_reltrans(21) = params_file(25)   
    
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !read benchmark files and output files 
    open(20,file="Benchmarks/Total_2,1_5,8.dat",status='old')
    open(21,file="Benchmarks/PivotingPL_2,1_5,8.dat",status='old')
    open(22,file="Benchmarks/PivotingReflection_2,1_5,8.dat",status='old')
    open(23,file="Benchmarks/LightTravelTime_2,1_5,8.dat",status='old')
    open(30,file="Output/Total.dat",status='old')
    open(31,file="Output/PivotingPL.dat",status='old')
    open(32,file="Output/PivotingReflection.dat",status='old')
    open(33,file="Output/LightTravelTime.dat",status='old')
    do i=1,ne
        read(20,*) (total_benchmark(i,j), j=1, 2)
        read(21,*) (pivoting_benchmark(i,j), j=1, 2)
        read(22,*) (refl_benchmark(i,j), j=1, 2)
        read(23,*) (reverb_benchmark(i,j), j=1, 2)        
        read(30,*) (total_model(i,j), j=1, 2)
        read(31,*) (pivoting_model(i,j), j=1, 2)
        read(32,*) (refl_model(i,j), j=1, 2)
        read(33,*) (reverb_model(i,j), j=1, 2)                           
    end do
    close(20) 
    close(21) 
    close(22) 
    close(23)     
    close(30) 
    close(31) 
    close(32) 
    close(33) 
   
    !compare total and individual components
    print *,"----------------------------------------------------------------"
    print *, "Comparing total model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = total_benchmark(i,2)/total_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Total model output different at energy bin", total_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", total_benchmark(i,2), total_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Total model output test passed"
    else                                    
        print*, "Total model output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting continuum output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = pivoting_benchmark(i,2)/pivoting_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting continuum output different at energy bin", pivoting_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", pivoting_benchmark(i,2), pivoting_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting continuum output test passed"
    else                                    
        print*, "Pivoting continuum output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting reflection output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = refl_benchmark(i,2)/refl_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting reflection output different at energy bin", refl_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", refl_benchmark(i,2), refl_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting reflection output test passed"
    else                                    
        print*, "Pivoting reflection output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing reverberation output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = reverb_benchmark(i,2)/reverb_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Reverberation output different at energy bin", reverb_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", reverb_benchmark(i,2), reverb_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Reverberation output test passed"
    else                                    
        print*, "Reverberation output test failed, recommend re-running the test with full_print set to true!"
    endif  
    
    print *,"----------------------------------------------------------------"
    print *,"Reltrans test; 5.8-16 Hz"  
    !get parameter file fron the benchmark folder for the first frequency interval
    open(1,file="Benchmarks/ip_5,8_16.dat",status='old')
    read(1,*) params_file
    close(1)    
    
    params_reltrans(1)  = params_file(1)
    params_reltrans(2)  = params_file(2)
    params_reltrans(3)  = params_file(3)
    params_reltrans(4)  = params_file(4)
    params_reltrans(5)  = params_file(5)
    params_reltrans(6)  = params_file(6)
    params_reltrans(7)  = params_file(7)
    params_reltrans(8)  = params_file(8)
    params_reltrans(9)  = params_file(9)
    params_reltrans(10) = params_file(10)
    params_reltrans(11) = params_file(11)
    params_reltrans(12) = params_file(12)
    params_reltrans(13) = params_file(13)
    params_reltrans(14) = params_file(14)
    params_reltrans(15) = params_file(18)
    params_reltrans(16) = params_file(19)
    params_reltrans(17) = params_file(20)
    params_reltrans(18) = params_file(21)
    params_reltrans(19) = params_file(22)
    params_reltrans(20) = params_file(23)
    params_reltrans(21) = params_file(25)   
    
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !read benchmark files and output files 
    open(20,file="Benchmarks/Total_5,8_16.dat",status='old')
    open(21,file="Benchmarks/PivotingPL_5,8_16.dat",status='old')
    open(22,file="Benchmarks/PivotingReflection_5,8_16.dat",status='old')
    open(23,file="Benchmarks/LightTravelTime_5,8_16.dat",status='old')
    open(30,file="Output/Total.dat",status='old')
    open(31,file="Output/PivotingPL.dat",status='old')
    open(32,file="Output/PivotingReflection.dat",status='old')
    open(33,file="Output/LightTravelTime.dat",status='old')
    do i=1,ne
        read(20,*) (total_benchmark(i,j), j=1, 2)
        read(21,*) (pivoting_benchmark(i,j), j=1, 2)
        read(22,*) (refl_benchmark(i,j), j=1, 2)
        read(23,*) (reverb_benchmark(i,j), j=1, 2)        
        read(30,*) (total_model(i,j), j=1, 2)
        read(31,*) (pivoting_model(i,j), j=1, 2)
        read(32,*) (refl_model(i,j), j=1, 2)
        read(33,*) (reverb_model(i,j), j=1, 2)                           
    end do
    close(20) 
    close(21) 
    close(22) 
    close(23)     
    close(30) 
    close(31) 
    close(32) 
    close(33) 
   
    !compare total and individual components
    print *,"----------------------------------------------------------------"
    print *, "Comparing total model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = total_benchmark(i,2)/total_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Total model output different at energy bin", total_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", total_benchmark(i,2), total_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Total model output test passed"
    else                                    
        print*, "Total model output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting continuum output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = pivoting_benchmark(i,2)/pivoting_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting continuum output different at energy bin", pivoting_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", pivoting_benchmark(i,2), pivoting_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting continuum output test passed"
    else                                    
        print*, "Pivoting continuum output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing pivoting reflection output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = refl_benchmark(i,2)/refl_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Pivoting reflection output different at energy bin", refl_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", refl_benchmark(i,2), refl_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Pivoting reflection output test passed"
    else                                    
        print*, "Pivoting reflection output test failed, recommend re-running the test with full_print set to true!"
    endif
    
    print *, "Comparing reverberation output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = reverb_benchmark(i,2)/reverb_model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                print*, "Reverberation output different at energy bin", reverb_benchmark(i,1), "keV;", &
                        "Benchmark and model output:", reverb_benchmark(i,2), reverb_model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        print*, "Reverberation output test passed"
    else                                    
        print*, "Reverberation output test failed, recommend re-running the test with full_print set to true!"
    endif       
    print *,"----------------------------------------------------------------"
                              
end program
