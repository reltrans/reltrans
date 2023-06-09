program relbenchgen
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
    !tbd: check here to make sure users actually want to run the benchmark!   
      
    call c_setgetenv("0.3","10.")
                            
    !test one: xrb, maxi j1820, single lamp post 
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
    !spec_flag = .true.  
    !mode = "agn"
    !frange = "xmm_01_06"
    !call model_singleLP(mode,frange,spec_flag)    
    !frange = "xmm_06_20"
    !call model_singleLP(mode,frange,spec_flag)        
    !frange = "xmm_20_90" 
    !call model_singleLP(mode,frange,spec_flag)
    
    !frange = "nus_01_06"
    !call model_singleLP(mode,frange,spec_flag)    
    !frange = "nus_06_20"
    !call model_singleLP(mode,frange,spec_flag)        
    !frange = "nus_20_90" 
    !call model_singleLP(mode,frange,spec_flag) 
    
    !test three: dual lamp post
    !spec_flag = .true.
    !mode = "dbl"
    !frange = "0,10_0,40"
    !call model_doubleLP(mode,frange,spec_flag)
    !frange = "0,50_0,60"
    !call model_doubleLP(mode,frange,spec_flag)
    !frange = "1,10_1,40"
    !call model_doubleLP(mode,frange,spec_flag)
    !frange = "3,00_4,20"
    !call model_doubleLP(mode,frange,spec_flag)
    !frange = "4,30_15,6"
    !call model_doubleLP(mode,frange,spec_flag) 
                   
    call CPU_TIME (time_end)
    print*, "Benchmark generation runtime:", time_end - time_start, "seconds"
    
end program

subroutine model_singleLP(mode,frange,spec_flag)
    implicit none
    
    character(len=:), allocatable   :: first_path, second_path, par_path, start_path, bench_path, command, mtype
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
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    !move the total over
    start_path = "Output/Total.dat "
    bench_path = "Benchmarks/" // trim(mode) // trim(mtype) // "Total_" // trim(frange) // "_temp.dat"
    command = "cp " // trim(start_path) // " " // trim(bench_path)
    call system(command)
    !then each component 
    start_path = "Output/IonVariations.dat "
    bench_path = "Benchmarks/" // trim(mode) // trim(mtype) // "IonVar_" // trim(frange) // "_temp.dat"
    command = "cp " // trim(start_path) // " " // trim(bench_path)
    call system(command)
    start_path = "Output/LightTravelTime.dat "
    bench_path = "Benchmarks/" // trim(mode) // trim(mtype) // "Reverb_" // trim(frange) // "_temp.dat"
    command = "cp " // trim(start_path) // " " // trim(bench_path)
    call system(command)  
    start_path = "Output/PivotingPL.dat "
    bench_path = "Benchmarks/" // trim(mode) // trim(mtype) // "PivPL_" // trim(frange) // "_temp.dat"
    command = "cp " // trim(start_path) // " " // trim(bench_path)
    call system(command)     
    start_path = "Output/PivotingReflection.dat "
    bench_path = "Benchmarks/" // trim(mode) // trim(mtype) // "PivRef_" // trim(frange) // "_temp.dat"
    command = "cp " // trim(start_path) // " " // trim(bench_path)
    call system(command)    
    
    !then compare modulus
    params_reltrans(18) = 3
    mtype = "/Mods/"
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !then imaginary part:
    mtype = "/Imag/"
    params_reltrans(18) = 2
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !then real:
    mtype = "/Real/"
    params_reltrans(18) = 1
    call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
    
    !finally calculate spectra and kernel
    if( spec_flag .eqv. .true. ) then
        params_reltrans(16) = 0
        params_reltrans(17) = 0
        params_reltrans(18) = 1
        mtype = "/Spec/"
        call tdreltransDCp(ear,ne,params_reltrans,ifl,photar)
        spec_flag = .false.
        !tbd: copy spectral and kernel files
    end if

    deallocate(ear)
    deallocate(photar)

    return 
end subroutine
