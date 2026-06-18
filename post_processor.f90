! gfortran post_processor.f90 -J ./build/cache -L"${HEADAS}/lib" -lcfitsio -lreltrans -L./build/lib

program ppmain
  implicit none
  real param(15)
  integer xe,adensity
  logical chainmode
  character (len=500) chainfile,newchainfile
  integer iparam(15),i

! Input parameters
  param(1)  = 11.4144          !h
  param(2)  = 0.998            !a
  param(3)  = 42.2816          !inc
  param(4)  = -2.95166         !rin
  param(5)  = 2e4              !rout
  param(6)  = 0.0              !zcos
  param(7)  = 1.70545          !Gamma
  param(8)  = 3.39811          !logxi
  param(9)  = 2.64304          !Afe
  param(10) = 19.0404          !lognep
  param(11) = 359.389          !kTe
  param(12) = 0.0              !Nh
  param(13) = 0.294556         !boost
  param(14) = 14.8223          !Mass
  param(15) = 0.0999336        !Anorm

! Settings
  chainmode = .true.  !Reading in a chain (true) or just entering one parameter set (false)
  xe        = 20       !Number of radial zones
  adensity  = 1        !1 = zone A ne; 0 = const ne
  
! Name of input chain
  chainfile = '/Users/nai47/Dropbox/Patrick_CygX1_RXTE/post_referee_fits/contour_plots/'
  chainfile = trim(chainfile) // 'adam_dcp_final.out'

! Name of output chain (with distance column added)
  newchainfile = 'adam_dcp_final_dist.out'

! Column number corresponding to each parameter in the chain
  !(not used if chainmode=false)
  !0 means that the parameter was fixed to the value in param(:)
  iparam(1)  = 5           !h
  iparam(2)  = 0           !a
  iparam(3)  = 6           !inc
  iparam(4)  = 7           !rin
  iparam(5)  = 0           !rout
  iparam(6)  = 0           !zcos
  iparam(7)  = 8           !Gamma
  iparam(8)  = 9           !logxi
  iparam(9)  = 10          !Afe
  iparam(10) = 11          !lognep
  iparam(11) = 12          !kTe
  iparam(12) = 0           !Nh
  iparam(13) = 13          !boost
  iparam(14) = 14          !Mass
  iparam(15) = 15          !Anorm

  call post_processor(param,xe,adensity,chainmode,chainfile,newchainfile,iparam)

end program ppmain
  

!-----------------------------------------------------------------------
subroutine post_processor(param,xe,adensity,chainmode,chainfile,newchainfile,  &
     iparam)
!> 
!> Takes a set of reltransDCp parameters and calculates the distance in kpc. 
!> If chainmode = false, it just writes the distance to terminal.
!> If chainmode = true, it also reads in a chain, calculates the distance for
!> every step of the chain, and creates a new fits file with a distance column.
!> Inputs:
!> param(15)     ReltransDCp parameters
!> xe            Number of radial zones
!> adensity      Radial density profile (=1 means zone A, =0 means constant)
!> chainmode     Just calculate once (false) or for an entire chain (true).
!> chainfile     Name of input chain file (must already exist)
!> newchainfile  Name of output chain file.
!> iparam(15)    Column number of that parameter in the chain. If parameter i
!>               was frozen for the chain, set iparam(i) = 0. The frozen
!>               parameters are then read from param(:), so the main parameter
!>               array is also required for chainmode=true.
!>
    implicit none
    real, intent(in)    :: param(15)
    integer, intent(in) :: xe,adensity,iparam(15)
    logical, intent(in) :: chainmode
    character (len=500), intent(in) :: chainfile,newchainfile
    real par(32)
    integer Cp,nlp
    double precision distance,Dkpc
    logical anynull
    integer unit,status,readwrite,blocksize
    character (len=500) comment
    integer steps,columns,col,k,newunit
    real chisquared,parray(15)
    integer i,colnum,ncols,nhdu,hdutype
    character (len=16) ttype(1),tform(1),tunit(1)
  
!   Calculate distance for input parameters
    call pack_reltransDCp(param,par,Cp,nlp)
    Dkpc = distance(Cp, nlp, xe, adensity, par)
    write(*,*)"reltransDCp distance (kpc) = ",Dkpc
  
!   Read in chain and append distance
    if( chainmode )then

        write(*,*)"Input chain file=",trim(chainfile)
        write(*,*)"Output chainfile=",trim(newchainfile)
     
        !Open chain file
        status = 0
        call ftgiou(unit,status)
        readwrite = 0
        call ftopen(unit,chainfile,readwrite,blocksize,status)
        if( status .ne. 0 ) stop 'cannot open chain file'
     
        !Shift to  extension "CHAIN"
        status = 0
        call ftmnhd(unit,2,'CHAIN',0,status)
        if( status .ne. 0 ) stop 'cannot shift to extension CHAIN'
     
        !Read number of rows and columns
        call ftgkyj(unit,'NAXIS2',steps,comment,status)
        if(status .ne. 0) stop 'Cannot determine No of rows'
        call ftgkyj(unit,'TFIELDS',columns,comment,status)
        if(status .ne. 0) stop 'Cannot determine No of columns'
     
        !Copy chain file to new chain file
        status = 0
        call copyhdu(unit,newchainfile,newunit)
        if( status .ne. 0 ) stop 'Could not copy chain file'
     
        !Add a new column to the new CHAIN binary table
        !First move to the CHAIN extension
        status = 0
        call ftmnhd(newunit,2,'CHAIN',0,status)
        if( status .ne. 0 ) stop 'Could not move to CHAIN ext of newchainfile'
        !Insert the name and type of the new column.
        status = 0
        colnum = columns + 1
        ncols  = 1
        ttype(1) = 'Dkpc'
        tform(1) = '1D'
        tunit(1) = 'kpc'
        call FTICLS(newunit,colnum,ncols,ttype,tform,status)
        if( status .ne. 0 ) stop 'Something wrong in FTICLS'
  
        !Go through each step and calculate distance for each one
        write(*,*)"Calculating distance and appending to chain..."
        do k = 1,steps
            status  = 0
            !Read in parameters
            call ftgcve(unit,columns,k,1,1,-1.0,chisquared,anynull,status)
            do i = 1,15
                if( iparam(i) .gt. 0 )then
                    call ftgcve(unit,iparam(i),k,1,1,-1.0,parray(i),anynull    &
                    ,status)
                else
                    parray(i) = param(i)
                end if
            end do
            !Calculate distance
            call pack_reltransDCp(parray,par,Cp,nlp)
            Dkpc = distance(Cp, nlp, xe, adensity, par)
            !Append to new chain file
            call ftpcld(newunit,columns+1,k,1,1,Dkpc,status)
        end do

    end if
  
!   Close chain files
    call ftclos(unit, status)
    call ftfiou(unit, status)
    call ftclos(newunit, status)
    call ftfiou(newunit, status)

    end subroutine post_processor
!-----------------------------------------------------------------------







! -----------------------------------------------------------------------
function distance(Cp, nlp, xe_in, adensity_in, param)
!
! Calculates the distance given a set of model parameters for reltransDCp.
! Can't do this for reltransPL because the density is hardwired to 15.
! In principle could also do for the double lamppost model, but I haven't
! yet implemented that.  
! Inputs:
! Cp        = PL/nthcomp flag. Only works for Cp=2.
! nlp       = No of lampposts. Code doesn't work for nlp=2.  
! param(32) = full genreltrans parameter array
! Outputs:
! distance  = distance in kpc
!
    use dyn_gr
    use conv_mod
    use radial_grids
    use gr_continuum
    use m_genreltrans
    use env_variables
    use saved_variables
    use telematrix2
    implicit none
    ! Constants
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, dlogf = 0.09 !This is a resolution parameter (base 10)
    ! Args:
    integer, intent(in) :: Cp, nlp, xe_in, adensity_in
    real, intent(inout) :: param(32)
    ! Variables of the subroutine
    ! initializer
    integer :: m, prev_nf, Cpsave, i, j, Cp_cont, dset, get_env_int
    double precision :: d, distance
    ! relativistic parameters and limit on rin and h
    ! lens needs to be allocatable to save it.
    double precision, allocatable :: frobs(:), frrel(:)
    double precision :: fhisave, flosave, fcons, contx_temp

    double precision, allocatable :: logxir_dens(:)

    real time_start, time_end !runtime stuff

    ! SAVE
    real, save :: paramsave(32)

    data Cpsave/2/
    data prev_nf /-1/
    ! Save the first call variables
    save d, frobs, frrel

    type(t_config), pointer :: config
    type(t_arrays), save :: arrays
    type(t_model_arguments) :: model_args
    ! make arrays static so its values are kept between function calls

    ! Setup
    config => global_config
    call unwrap_arguments(model_args, nlp, dset, param, Cp)
    call config_frequency(config, model_args)
    call arguments_check(config, model_args)

    if (config%firstcall) then
        call init_fftw_allconv()
        config%firstcall = .false.
        config%needtrans = .true.
        config%needconv  = .true.
        prev_nf = 0 !this is needed to reallocate arrays with realloc_arrays, if firstcall is set to true externally
        ! set sensible distance for observer from the BH
        d = max(1.0d4, 2.0d2 * config%rnmax**2)
        spinsav = -2.d0 !this is needed to force the run of the GRtrace routine
        config%verbose = get_env_int("REV_VERB", 0)
        ! initialise environment and allocate all arrays
        config%xe      = xe_in
        config%me      = 1
        adensity       = adensity_in        
        call setup_global_arrays(config, model_args%nlp)
    end if

    ! Allocate arrays that need to be re-allocated every time incase xe changed
    config%xe      = xe_in
    adensity       = adensity_in        
    if (allocated(dfer_arr)) deallocate(dfer_arr)
    if (allocated(logxir)) deallocate(logxir)
    if (allocated(gsdr)) deallocate(gsdr)
    if (allocated(logner)) deallocate(logner)
    allocate(dfer_arr(config%xe))
    allocate(logxir (config%xe))
    allocate(gsdr (config%xe))
    allocate(logner (config%xe))    
    call setup_arrays(config, arrays, model_args%nlp)
    
    ! reallocated frequency dependent arrays
    call realloc_arrays(config, model_args, arrays, prev_nf)

    ! Note: the two different calls are because for the double lP we set the
    ! temperature from the coronal frame(s), but for the single
    ! LP we use the temperature in the observer frame

    ! Decide if this is the DC component/time averaged spectrum or not
    if (config%flo .lt. tiny(config%flo) .or. config%fhi .lt. tiny(config%fhi))then
        config%DC = 1
        model_args%g = 0.0
        model_args%DelAB = 0.0
        model_args%DelA = 0.0
        model_args%ReIm = 1
        model_args%eta = model_args%eta_0
        ! this is an ugly hack for the double LP model to calculate the time-
        ! averaged spectrum
        model_args%beta_p = 1.
    else
        config%DC = 0
        model_args%boost = abs(model_args%boost)
    end if

    ! Determine if I need to calculate the kernel
    call need_check(model_args%Cp, Cpsave, param, paramsave, config%fhi,       &
         config%flo, fhisave, flosave, config%nf, prev_nf, config%needtrans,   &
         config%needconv)

    ! Calculate the kernel
    if (config%verbose .gt. 2) call CPU_TIME (time_start)
    if (config%needtrans)then
       ! allocate lensing/reflection fraction arrays if necessary
       if (allocated(lens)) deallocate(lens)
       allocate (lens(nlp))
       if (allocated(frobs)) deallocate(frobs)
       allocate (frobs(nlp))
       if (allocated(frrel)) deallocate(frrel)
       allocate (frrel(nlp))
       ! Calculate the Kernel for the given parameters
       call rtrans(config, model_args, arrays, dset, d, nex, frobs, frrel)
       ! print *, 'gso ', gso(1)
    end if
    if (config%verbose .gt. 2) then
       call CPU_TIME (time_end)
       print *, 'Transfer function runtime: ', time_end - time_start, ' seconds'
    end if

! Calculate ionization profile measured by reltransDCp
    dset = 0
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    call radfunctions_dens(config, model_args, arrays)
    allocate( logxir_dens(config%xe) )
    logxir_dens = logxir

! Calculate ionization profile assuming isotropic corona seen from D=1kpc
    dset            = 1
    model_args%Dkpc = 1.0
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    call radfuncs_dist(config%xe, model_args%rin, rnmax,model_args%b1,         &
         model_args%b2, model_args%qboost, fcons,                              &
         & dble(model_args%lognep), model_args%a, model_args%h(1),             &
         model_args%honr, rlp, dcosdr, cosd, ndelta, config%rmin,npts(1),      &
         & logxir, gsdr, logner, pnorm)
    
! Caculate the re-scaling distance
    Distance = 10.0**( 0.5 * ( logxir_dens(2)-logxir(2) )  )
    Distance = Distance / sqrt(model_args%boost)

! Save parameters
    fhisave = config%fhi
    flosave = config%flo
    prev_nf = config%nf
    paramsave = param
    Cpsave = model_args%Cp
    
  end function distance  
! -----------------------------------------------------------------------




  


!-----------------------------------------------------------------------
subroutine pack_reltransDCp(param,par,Cp,nlp)
  implicit none
  real, intent(in)     :: param(15)
  real, intent(out)    :: par(32)
  integer, intent(out) :: Cp,nlp
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter
  nlp  = 1   !Number of lampposts
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = 0.0              !floHz
  par(24) = 0.0              !fhiHz
  par(25) = 1.0              !ReIm
  par(26) = 0.0              !DelA
  par(27) = 0.0              !DelAB
  par(28) = 0.0              !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(15)        !Anorm
  par(32) = 1.0              !telescope response  
  return
end subroutine pack_reltransDCp
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine copyhdu(inunit,outfile,outunit)
! Inputs
! inunit:  unit number of (alreay opened) input file.
! outfile: name of output file to be created
! Outputs
! outunit: unit number of output file
  implicit none
  integer status,inunit,outunit,blocksize,morekeys,hdutype
  character (len=500) outfile
  
! The STATUS parameter must always be initialized.
  status=0

! Delete the file if it already exists, so we can then recreate it
! The deletefile subroutine is listed at the end of this file.
  call deletefile(outfile,status)

! Get  unused Logical Unit Numbers to use to open the FITS file.
  call ftgiou(outunit,status)

! Create the new empty FITS file (value of blocksize is ignored)
  blocksize=1
  call ftinit(outunit,outfile,blocksize,status)

! Skip to the 2nd extension in the input file
  call ftmahd(inunit,2,hdutype,status)
  
! FTCOPY copies the current HDU from the input FITS file to the output
! file.  The MOREKEY parameter allows one to reserve space for additional
! header keywords when the HDU is created.   FITSIO will automatically
! insert more header space if required, so programmers do not have to
! reserve space ahead of time, although it is more efficient to do so if
! it is known that more keywords will be appended to the header.
  morekeys=0
  call ftcopy(inunit,outunit,morekeys,status)

 end subroutine copyhdu
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine deletefile(filename,status)

!C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!C  Simply return if status is greater than zero
      if (status .gt. 0)return

!C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
!C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!C  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
!-----------------------------------------------------------------------
