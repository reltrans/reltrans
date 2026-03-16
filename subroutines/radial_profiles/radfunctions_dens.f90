!-----------------------------------------------------------------------
subroutine radfunctions_dens(config, model_args, arrays)
    ! In  : config      - global configuration (e.g. number of radial grid points: config%xe)
    !       model_args  - model parameters (e.g. geometry, logxi, lognep, honr, number of LPs: model_args%nlp)
    !       arrays      - input arrays bundled in t_arrays (if applicable; main radial profiles come from modules)
    !      Uses/updates module arrays from radial_grids:
    !       logxir(1:config%xe) - log10 of ionization parameter as a function of radius
    !       gsdr(1:config%xe)   - source-to-disc blueshift factor as a function of radius
    !       logner(1:config%xe) - log10 of electron density as a function of radius
    !       dfer_arr(1:config%xe) - emissivity-related radial scaling array
    use common_types
    use env_variables
    use dyn_gr, only: ndelta, rlp, dcosdr, cosd, npts
    use radial_grids, only: logxir, gsdr, logner, dfer_arr
    implicit none
    type(t_config)         , intent(in)    :: config
    type(t_model_arguments), intent(in)    :: model_args
    type(t_arrays)         , intent(in)    :: arrays

    real                           :: gso(model_args%nlp)
    double precision :: rlp_column(ndelta),dcosdr_column(ndelta),cosd_column(ndelta), dgsofac
    integer          :: i, kk, get_index, get_env_int, l, m
    double precision :: rp, logxinorm, lognenorm,  mus, interper, newtex, mui, dinang, gsd(model_args%nlp), dglpfacthick
    double precision :: xi_lp(config%xe,model_args%nlp), logxi_lp(config%xe,model_args%nlp), logxip_lp(model_args%nlp)
    double precision :: xitot, xiraw, mylogne, mudisk, gsd_temp
    double precision, allocatable :: rad(:)

    ! Set disk opening angle
    mudisk   = model_args%honr / sqrt( model_args%honr**2 + 1.d0  )
    
    allocate(rad(config%xe))
    !Now calculate logxi itself
    ! The loop calculates the raw xi and raw n_e.
    ! This means they are without normalization: only to find the maximum and the minimum. Remember that the max of the ionisation is not the same as the minumim in the density because the flux depends on r
    !The loops calculates also the correction factor mui

    !TBD: include luminosity ratio between LPs 
    do i = 1, config%xe
        rad(i) = (config%rnmax/model_args%rin)**(real(i-1) / real(config%xe))
        rad(i) = rad(i) + (config%rnmax/model_args%rin)**(real(i) / real(config%xe))
        rad(i) = rad(i) * model_args%rin * 0.5
        !Initialize total ionization tracker
        xitot = 0. 
        gsd_temp = 0.
        !Now calculate the raw density (this matters only for high dens model reltransD)
        logner(i) = adensity * mylogne(rad(i), model_args%rin)
        do m=1,model_args%nlp
            do l=1,ndelta
                rlp_column(l) = rlp(l,m)
                dcosdr_column(l) = dcosdr(l,m)
                cosd_column(l) = cosd(l,m)
            end do    
            gso(m) = real( dgsofac(model_args%a, model_args%h(m)) )
            xi_lp(i,m) = xiraw(rad(i), model_args%a, model_args%h(m),         &
                model_args%honr, rlp_column, dcosdr_column, ndelta,            &
                config%rmin, npts(m), mudisk, gsd(m))
            if (m .eq. 2) xi_lp(i,m) = model_args%eta_0 * xi_lp(i,m)
            !Calculate the incident angle for this bin
            kk = get_index(rlp_column, ndelta, rad(i), config%rmin, npts(m))
            mus = interper(rlp_column, cosd_column, ndelta, rad(i), kk)
            if( kk .eq. npts(m) ) mus = newtex(rlp_column, cosd_column,        &
                ndelta, rad(i), model_args%h(m), model_args%honr, kk)
            mui = dinang(model_args%a, rad(i), model_args%h(m), mus)
            !Correction to account for the radial dependence of incident angle, and for the g factors
            xi_lp(i,m) = xi_lp(i,m) / (sqrt(2.) * mui)                        &
                * arrays%contx_int(m) * (gso(m))**(model_args%Gamma - 2)
            xitot = xitot + xi_lp(i,m)
            gsd_temp = gsd_temp + gsd(m)*xi_lp(i,m)
        end do 
        !This and the line above calculate the gsd factor along the disk, averaging over the flux the disk sees from each LP 
        gsdr(i) = gsd_temp/xitot
        logxir(i) = log10(xitot) - logner(i)     
    end do
    !After the loop calculate the max and the min - ionization renormalized wrt to the first LP
    logxinorm = maxval(logxir)
    lognenorm = minval(logner)
    logxir = logxir - (logxinorm - dble(model_args%logxi))
    logner = logner - (lognenorm - dble(model_args%lognep))
    
    do m=1,model_args%nlp
        do i=1,config%xe
            logxi_lp(i,m) = log10(xi_lp(i,m)) - logner(i) - lognenorm         &
                - logxinorm + dble(model_args%lognep) + dble(model_args%logxi)
        end do
        logxip_lp(m) = max(maxval(logxi_lp(:,m)),0.)
    end do
    
    !Write radii, ionisation (for both and each LP), gamma factors, and log(xi(r))+log(ne(r)) (which is nearly the same as
    !epsilon(r) for identical coronal spectrra and gamma=2) to file. 
    !note 1) we need to do this before the ionisation array is set to have a minimum of 0, in order
    !to recover the correct scaling of the emissivity at large radii
    !2) in order to correctly compare the dfer_arr array with the single LP case, it has to be renormalized by (1+eta_0)
    if( config%verbose .gt. 1 ) then
        print*, "Peak ionisations from each LP: first " , logxip_lp(1), " second ", logxip_lp(2)
        open (unit = 27, file = 'Output/RadialScalings.dat', status='replace', action = 'write')
            do i = 1, config%xe
                write(27,*) rad(i), logxir(i), gsdr(i), logxir(i)+logner(i), logxi_lp(i,1), logxi_lp(i,2), dfer_arr(i) 
            end do 
        close(27)    
    end if
    
    !check max and min for ionisation 
    logxir = max( logxir , 0.d0  )
    logxir = min( logxir , 4.7d0 )
    
    deallocate(rad)

    return
end subroutine radfunctions_dens
!-----------------------------------------------------------------------
