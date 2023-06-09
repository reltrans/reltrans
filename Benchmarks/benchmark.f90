program relbench
    implicit none    
    !This interfaces the C routine that sets the environment variables correctly
    !Necessary because the fortran "system" call somehow can't do this
    interface
        subroutine c_setgetenv (emin,emax) bind(c)
            use, intrinsic :: iso_c_binding
            implicit none
            CHARACTER(c_char) :: emin, emax
            !import ! use declarations from host (implicit none, iso_c_binding)
        end subroutine c_setgetenv
    end interface
    
    real    :: emin, emax           !minimum, maximum energy and increment to set model grid
    integer :: i , j                !loop variables 
    real time_start,time_end        !runtime stuff                                    
                           
    character(len=:), allocatable   :: mode, frange  !strings to call the model appriopriately
    logical spec_flag               !flag to call time-averaged model too
                                         
    call CPU_TIME (time_start)
    open(60,file='Benchmarks/Benchmark_result.txt',status='replace', action = 'write')
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Setting up" 
    write (60,*) "MAKE SURE TO RUN THE TEST WITH FNINT COMMENTED OUT FROM SREVTRANS.F90!!!!!!!!!!!!!!!!!"
    write (60,*) "Note: this test is very strict, and some benchmarks may report as being failed due to numerical"
    write (60,*) "precision issues in fftw. If the test reports discrepancy in only a handful of bins, consider the "
    write (60,*) "test passed."    
      
    call c_setgetenv("0.3","10.")
                            
    !test one: xrb, maxi j1820, single lamp post 
    write (60,*)  "--------------------------------------------------------------------------------------------------------"
    write (60,*)  "XRB test: "
    spec_flag = .true.  
    mode = "xrb"
    frange = "0,12_0,25" 
    call model_singleLP(mode,frange,spec_flag)    
    frange = "0,31_0,73" 
    call model_singleLP(mode,frange,spec_flag)          
    frange = "0,80_2,10" 
    call model_singleLP(mode,frange,spec_flag)    
    frange = "2,10_5,80" 
    call model_singleLP(mode,frange,spec_flag)    
    frange = "5,80_16,0" 
    call model_singleLP(mode,frange,spec_flag)
    
    !test two: agn, ark 564, single lamp post, xmm first,nustar second 
    write (60,*)  "--------------------------------------------------------------------------------------------------------"   
    write (60,*)  "AGN test: "
    spec_flag = .true.  
    mode = "agn"
    frange = "xmm_01_06"
    call model_singleLP(mode,frange,spec_flag)    
    frange = "xmm_06_20"
    call model_singleLP(mode,frange,spec_flag)        
    frange = "xmm_20_90" 
    call model_singleLP(mode,frange,spec_flag)
    
    frange = "nus_01_06"
    call model_singleLP(mode,frange,spec_flag)    
    frange = "nus_06_20"
    call model_singleLP(mode,frange,spec_flag)        
    frange = "nus_20_90" 
    call model_singleLP(mode,frange,spec_flag) 
    
    !test three: dual lamp post 
    write (60,*)  "--------------------------------------------------------------------------------------------------------"  
    write (60,*)  "Dual LP test: "
    spec_flag = .true.
    mode = "dbl"
    frange = "0,10_0,40"
    call model_doubleLP(mode,frange,spec_flag)
    frange = "0,50_0,60"
    call model_doubleLP(mode,frange,spec_flag)
    frange = "1,10_1,40"
    call model_doubleLP(mode,frange,spec_flag)
    frange = "3,00_4,20"
    call model_doubleLP(mode,frange,spec_flag)
    frange = "4,30_15,6"
    call model_doubleLP(mode,frange,spec_flag) 
                   
    call CPU_TIME (time_end)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Benchmarks runtime:", time_end - time_start, "seconds"
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    close(60) 

end program

subroutine model_singleLP(mode,frange,spec_flag)
    implicit none
    
    character(len=:), allocatable   :: first_path, second_path, par_path, mtype
    character(len=3)                :: mode
    character(len=9)                :: frange
    
    integer     :: ne, i                   
    integer     :: ifl
    
    real        :: emin, emax
    real        :: params_reltrans(22)
    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array
    
    logical spec_flag    
    
    !from the path, find out automatically the path to the parameter files, load into reltrans, run time-averaged spectrum+
    !lags+cross spectral amplitude and compare to benchmark (also in the path)
    !then compare the results to the path in the benchmarks folder
    first_path = "Benchmarks/" // trim(mode) // "/ip_"
    second_path = '.dat'
    par_path = trim(first_path) // trim(frange) // trim(second_path)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) par_path
    
    open(50,file=par_path,status='old')
    read(50,*) params_reltrans
    close(50)
    ne = 1000                             
    allocate(ear(0:ne))  
    allocate(photar(ne))       
    emin = 0.1
    emax = 200.
    do i=0,ne
        ear(i) = emin * (emax/emin)**(real(i)/real(ne))
    end do  
    
    mtype = "/Lags/"
    params_reltrans(18) = 4
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running lag test: "
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)    
    
    !then compare modulus
    params_reltrans(18) = 3
    mtype = "/Mods/"
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running modulus test: "
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)
    
    !then imaginary part:
    mtype = "/Imag/"
    params_reltrans(18) = 2
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running imaginary test: "
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)    
    
    !then real:
    mtype = "/Real/"
    params_reltrans(18) = 1
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running real test: "
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)  
    
    !finally compare spectra and kernel 
    if( spec_flag .eqv. .true. ) then
        params_reltrans(16) = 0
        params_reltrans(17) = 0
        params_reltrans(18) = 1
        mtype = "/Spec/"
        write (60,*) "--------------------------------------------------------------------------------------------------------"
        write (60,*) "Running time-averaged spectrum test: " 
        call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
        write (60,*) "--------------------------------------------------------------------------------------------------------"
        call compare_spectrum(ne,mode,mtype)
        call compare_kernel(mode,mtype) 
        spec_flag = .false.
    end if

    deallocate(ear)
    deallocate(photar)

    return 
end subroutine

subroutine model_doubleLP(mode,frange,spec_flag)
    implicit none
    
    character(len=:), allocatable   :: first_path, second_path, par_path, mtype
    character(len=3)                :: mode
    character(len=9)                :: frange
    
    integer     :: ne, i                   
    integer     :: ifl
    
    real        :: emin, emax
    real        :: params_reltrans(27)
    
    real, dimension(:), allocatable      :: ear        !photon energy grid 
    real, dimension(:), allocatable      :: photar     !spectrum array
    
    logical spec_flag    
    
    !from the path, find out automatically the path to the parameter files, load into reltrans, run time-averaged spectrum+
    !lags+cross spectral amplitude and compare to benchmark (also in the path)
    !then compare the results to the path in the benchmarks folder
    first_path = "Benchmarks/" // trim(mode) // "/ip_"
    second_path = '.dat'
    par_path = trim(first_path) // trim(frange) // trim(second_path)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) par_path
    
    open(50,file=par_path,status='old')
    read(50,*) params_reltrans
    close(50)
    ne = 1000                             
    allocate(ear(0:ne))  
    allocate(photar(ne))       
    emin = 0.1
    emax = 200.
    do i=0,ne
        ear(i) = emin * (emax/emin)**(real(i)/real(ne))
    end do  
    
    mtype = "/Lags/"
    params_reltrans(21) = 4
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running lag test: "
    call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)    
    
    !then compare modulus
    params_reltrans(21) = 3
    mtype = "/Mods/"
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running modulus test: "
    call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)
    
    !then imaginary part:
    mtype = "/Imag/"
    params_reltrans(21) = 2
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running imaginary test: "
    call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)    
    
    !then real:
    mtype = "/Real/"
    params_reltrans(21) = 1
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    write (60,*) "Running real test: "
    call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)
    write (60,*) "--------------------------------------------------------------------------------------------------------"
    call compare_timing(ne,mode,mtype,frange)  
    
    !finally compare spectra
    if( spec_flag .eqv. .true. ) then
        params_reltrans(19) = 0
        params_reltrans(20) = 0
        params_reltrans(21) = 1
        mtype = "/Spec/"
        write (60,*) "--------------------------------------------------------------------------------------------------------"
        write (60,*) "Running time-averaged spectrum test: " 
        call tdreltransDbl(ear,ne,params_reltrans,ifl,photar)
        write (60,*) "--------------------------------------------------------------------------------------------------------"
        call compare_spectrum(ne,mode,mtype)
        call compare_kernel(mode,mtype) 
        spec_flag = .false.
    end if

    deallocate(ear)
    deallocate(photar)

    return 
end subroutine

subroutine compare_timing(ne,mode,mtype,frange)
    implicit none

    integer ne, i  
    real test_precision, test_model     
    logical test_bool, full_print  

    character(len=3) mode
    character(len=9) frange
    character(len=5) mtype
    character(len=50) full_path
    character(len=50) output_path

    real    :: benchmark(1000,2)      !total model from the benchmark file
    real    :: model(1000,2)          !total model from the model file
    
    output_path = "Output/Total.dat"
    open(50,file=output_path,status='old')
    do i=1,ne 
        read(50,*) model(i,1),model(i,2)
    end do
    close(50)
    
    full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/Total_" // trim(frange) // ".dat"
    open(50,file=full_path,status='old')
    do i=1,ne 
        read(50,*) benchmark(i,1),benchmark(i,2)
    end do
    close(50)
    test_precision = 1.001 
    full_print = .true.  
    
    write (60,*) "Comparing timing model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = benchmark(i,2)/model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                write (60,*) "Total model output different at energy bin", benchmark(i,1), "keV; ", &
                             "Benchmark and model output:", benchmark(i,2), model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        write (60,*) "Timing model output test passed"
    else                                    
        write (60,*) "Timing model output test failed, comparing model components"
    endif  
    if( test_bool .eqv. .false. ) then  
        output_path = "Output/IonVariations.dat"
        open(50,file=output_path,status='old')
        do i=1,ne 
            read(50,*) model(i,1),model(i,2)
        end do
        close(50)
        
        full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/IonVar_" // trim(frange) // ".dat"
        open(50,file=full_path,status='old')
        do i=1,ne 
            read(50,*) benchmark(i,1),benchmark(i,2)
        end do
        close(50)  
        write (60,*) "Comparing ionization model output: "
        test_bool = .true.
        do i=1,ne  
            test_model = benchmark(i,2)/model(i,2)
            if (abs(test_model) .ge. test_precision) then
                if (full_print .eqv. .true.) then
                    write (60,*) "Ionization changes output different at energy bin", benchmark(i,1), "keV; ", &
                                 "Benchmark and model output:", benchmark(i,2), model(i,2)
                endif
                test_bool = .false.
            endif
        end do    
        if (test_bool .eqv. .true.) then
            write (60,*) "Ionization output test passed"
        else                                    
            write (60,*) "Ionization output test failed"
        endif  

        output_path = "Output/PivotingPL.dat"
        open(50,file=output_path,status='old')
        do i=1,ne 
            read(50,*) model(i,1),model(i,2)
        end do
        close(50)
        
        full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/PivPL_" // trim(frange) // ".dat"
        open(50,file=full_path,status='old')
        do i=1,ne 
            read(50,*) benchmark(i,1),benchmark(i,2)
        end do
        close(50)    
        write (60,*) "Comparing pivoting PL output: "
        test_bool = .true.
        do i=1,ne  
            test_model = benchmark(i,2)/model(i,2)
            if (abs(test_model) .ge. test_precision) then
                if (full_print .eqv. .true.) then
                    write (60,*) "Pivoting changes output different at energy bin", benchmark(i,1), "keV; ", &
                                 "Benchmark and model output:", benchmark(i,2), model(i,2)
                endif
                test_bool = .false.
            endif
        end do    
        if (test_bool .eqv. .true.) then
            write (60,*) "Pivoting PL output test passed"
        else                                    
            write (60,*) "Pivoting PL output test failed"
        endif  
         !!!       
        output_path = "Output/PivotingReflection.dat"
        open(50,file=output_path,status='old')
        do i=1,ne 
            read(50,*) model(i,1),model(i,2)
        end do
        close(50)
        
        full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/PivRef_" // trim(frange) // ".dat"
        open(50,file=full_path,status='old')
        do i=1,ne 
            read(50,*) benchmark(i,1),benchmark(i,2)
        end do
        close(50) 
        write (60,*) "Comparing pivoting reflection output: "
        test_bool = .true.
        do i=1,ne  
            test_model = benchmark(i,2)/model(i,2)
            if (abs(test_model) .ge. test_precision) then
                if (full_print .eqv. .true.) then
                    write (60,*) "Pivoting reflection changes output different at energy bin", benchmark(i,1), "keV; ", &
                                 "Benchmark and model output:", benchmark(i,2), model(i,2)
                endif
                test_bool = .false.
            endif
        end do    
        if (test_bool .eqv. .true.) then
            write (60,*) "Pivoting reflection test passed"
        else                                    
            write (60,*) "Pivoting reflection output test failed"
        endif  
        
        output_path = "Output/LightTravelTime.dat"
        open(50,file=output_path,status='old')
        do i=1,ne 
            read(50,*) model(i,1),model(i,2)
        end do
        close(50)
        
        full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/Reverb_" // trim(frange) // ".dat"
        open(50,file=full_path,status='old')
        do i=1,ne 
            read(50,*) benchmark(i,1),benchmark(i,2)
        end do
        close(50)
        write(60,*) "Comparing reverberation output: "
        test_bool = .true.
        do i=1,ne  
            test_model = benchmark(i,2)/model(i,2)
            if (abs(test_model) .ge. test_precision) then
                if (full_print .eqv. .true.) then
                    write (60,*) "Reverberation changes output different at energy bin", benchmark(i,1), "keV; ", &
                                 "Benchmark and model output:", benchmark(i,2), model(i,2)
                endif
                test_bool = .false.
            endif
        end do    
        if (test_bool .eqv. .true.) then
            write (60,*) "Reverberation output test passed"
        else                                    
            write (60,*) "Reverberation output test failed"
        endif               
    end if  

    return 
end subroutine

subroutine compare_spectrum(ne,mode,mtype)
    implicit none

    integer ne, i   
    real test_precision, test_model     
    logical test_bool, full_print      

    character(len=3) mode
    character(len=9) frange
    character(len=5) mtype
    character(len=50) full_path
    character(len=50) output_path

    real    :: benchmark(1000,2)      !total model from the benchmark file
    real    :: model(1000,2)          !total model from the model file
    
    output_path = "Output/Total.dat"
    open(50,file=output_path,status='old')
    do i=1,ne 
        read(50,*) model(i,1),model(i,2)
    end do
    close(50)
    
    full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/Total.dat" 
    open(50,file=full_path,status='old')
    do i=1,ne 
        read(50,*) benchmark(i,1),benchmark(i,2)
    end do
    close(50)
    
    test_precision = 1.001 
    full_print = .true. 
    
    write (60,*) "Comparing spectrum model output: "
    test_bool = .true.     
    do i=1,ne  
        test_model = benchmark(i,2)/model(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                write (60,*) "Total model output different at energy bin", benchmark(i,1), model(i,1), "keV; ", &
                             "Benchmark and model output:", benchmark(i,2), model(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        write (60,*) "Model spectrum output test passed"
    else                                    
        write (60,*) "Model spectrum output test failed"
    endif     

    return 
end subroutine

subroutine compare_kernel(mode,mtype)
    implicit none

    !tbd: hardcode size of transfer functions
    integer i   
    integer, parameter :: ne = 2**12
    integer, parameter :: nt = 2**9
    real test_precision, test_model     
    logical test_bool, full_print      

    character(len=3) mode
    character(len=9) frange
    character(len=5) mtype
    character(len=60) full_path
    character(len=60) output_path

    real    :: benchmark_energy(ne,2)      
    real    :: model_energy(ne,2) 
    real    :: benchmark_time(ne,2)      
    real    :: model_time(ne,2)              
    
    output_path = "Output/Impulse_1dImpulseVsEnergy.dat"  
    open(50,file=output_path,status='old')
    do i=1,ne 
        read(50,*) model_energy(i,1),model_energy(i,2)
    end do
    close(50)
    
    full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/Impulse_1dImpulseVsEnergy.dat" 
    open(50,file=full_path,status='old')
    do i=1,ne 
        read(50,*) benchmark_energy(i,1),benchmark_energy(i,2)
    end do
    close(50)
    
    test_precision = 1.001 
    full_print = .true. 
    
    write (60,*) "Comparing kernel energy dependence: "
    test_bool = .true.     
    do i=1,ne  
        test_model = benchmark_energy(i,2)/model_energy(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                write (60,*) "Total model output different at g-factor", benchmark_energy(i,1), model_energy(i,1), "; ", &
                             "Benchmark and model output:", benchmark_energy(i,2), model_energy(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        write (60,*) "Kernel energy dependence test passed"
    else                                    
        write (60,*) "Kernel energy dependence test failed"
    endif     
    
    output_path = "Output/Impulse_1dImpulseVsTime.dat"  
    open(50,file=output_path,status='old')
    do i=1,nt 
        read(50,*) model_time(i,1),model_time(i,2)
    end do
    close(50)
    
    full_path = "Benchmarks/" // trim(mode) // trim(mtype) // "/Impulse_1dImpulseVsTime.dat" 
    open(50,file=full_path,status='old')
    do i=1,nt 
        read(50,*) benchmark_time(i,1),benchmark_time(i,2)
    end do
    close(50)
    
    write (60,*) "Comparing kernel time dependence: "
    test_bool = .true.     
    do i=1,nt  
        test_model = benchmark_time(i,2)/model_time(i,2)
        if (abs(test_model) .ge. test_precision) then
            if (full_print .eqv. .true.) then
                write (60,*) "Total model output different at g-factor", benchmark_time(i,1), model_time(i,1), "; ", &
                             "Benchmark and model output:", benchmark_time(i,2), model_time(i,2)
            endif
            test_bool = .false.
        endif
    end do    
    if (test_bool .eqv. .true.) then
        write (60,*) "Kernel time dependence test passed"
    else                                    
        write (60,*) "Kernel time dependence test failed"
    endif     

end subroutine 
