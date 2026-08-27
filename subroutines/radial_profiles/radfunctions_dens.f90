!-----------------------------------------------------------------------
subroutine radfunctions_dens(verbose,xe,rin,rnmax,eta_0,logxip,lognep,spin,h,Gamma,honr,rlp,dcosdr&
     &,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner,dfer_arr)
    ! In  : xe,rin,rnmax,eta_0,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
    ! logxir(xe),gsdr(xe)   logxi (ionization parameter) and gsd (source to disc blueshift) as a function of radius
    ! Out : logxir(1:xe), gsdr(1:xe), logner(1:xe)
    use env_variables
    use isco 
    implicit none
    integer         , intent(IN)   :: xe, ndelta, nlp, npts(nlp)
    double precision, intent(IN)   :: rin,rmin,rnmax,eta_0,logxip,lognep,spin,h(nlp),honr,Gamma,dfer_arr(xe)
    real                           :: gso(nlp)
    double precision, intent(IN)   :: rlp(ndelta,nlp), dcosdr(ndelta,nlp), cosd(ndelta,nlp), contx_int(nlp)
    double precision :: rlp_column(ndelta),dcosdr_column(ndelta),cosd_column(ndelta), dgsofac
    double precision, intent(INOUT):: logxir(xe), gsdr(xe), logner(xe)
    integer          :: i, kk, get_index, get_env_int, l, m, verbose
    double precision :: rp, logxinorm, lognenorm,  cosd_interp, interper, newtex, mui, dinang, gsd(nlp), dglpfacthick
    double precision :: xi_lp(xe,nlp), logxi_lp(xe,nlp), logxip_lp(nlp), xitot, xiraw, mylogne_zoneB, mudisk, gsd_temp
    double precision, allocatable :: rad(:)
    double precision :: xi_isco, logxir_isco, mui_isco, xitot_isco, rnorm
    
    ! Set disk opening angle
    mudisk   = honr / sqrt( honr**2 + 1.d0  )
    
    allocate(rad(xe))
    !Now calculate logxi itself
    ! The loop calculates the raw xi and raw n_e.
    ! This means they are without normalization: only to find the maximum and the minimum. Remember that the max of the ionisation is not the same as the minumim in the density because the flux depends on r
    !The loops calculates also the correction factor mui
    
    !TBD: include luminosity ratio between LPs 
    ! write(*,*) "DEBUG: ------- in the profile calculation loop ----------"
    !     write(*,*) "DEBUG:    rad(i),   xitot,     logxir(i),    logner(i)"
    do i = 1, xe        
        rad(i) = (rnmax/rin)**(real(i-1) / real(xe))
        rad(i) = rad(i) + (rnmax/rin)**(real(i) / real(xe))
        rad(i) = rad(i) * rin * 0.5
        !Initialize total ionization tracker
        xitot = 0. 
        gsd_temp = 0.

        !Now calculate the raw density (this matters only for high dens model reltransD)

        if (rad(i) .ge. risco) then
           logner(i) = adensity * mylogne_zoneB(rad(i), rin)
        else 
           logner(i) = adensity * mylogne_zoneB(risco, rin) + log10( (risco/rad(i))**(6.0/((5.0/3.0)+1.0))*&
                ((1.d0/(10.0**(-2.0d0)*sqrt(1.5d0*risco)))*(risco/rad(i)-1.d0)**1.5d0 +1.d0)**(-2.d0/((5.0/3.0)+1.0)))

        endif

        ! if (adensity .lt. 2) then
        !    logner(i) = adensity * mylogne_zoneB(rad(i), rin)
        ! else
        !    delete = ((1.d0/(10.0**(-2.0d0)*sqrt(1.5d0*risco)))*(risco/rad(i)-1.d0)**1.5d0 +1.d0)**(-2.d0/((5.0/3.0)+1.0))
        !    write(*,*) "DEBUG: risco", risco, (risco/rad(i))**(6.0/((5.0/3.0)+1.0)), delete
                
        !    logner(i) = log10( (risco/rad(i))**(6.0/((5.0/3.0)+1.0))*&
        !         ((1.d0/(10.0**(-2.0d0)*sqrt(1.5d0*risco)))*(risco/rad(i)-1.d0)**1.5d0 +1.d0)**(-2.d0/((5.0/3.0)+1.0)))
        !    ! radial density profile of Mummery and Balbus 2023 eq. 55
        ! endif
        
        do m=1,nlp
            do l=1,ndelta
                rlp_column(l) = rlp(l,m)
                dcosdr_column(l) = dcosdr(l,m)
                cosd_column(l) = cosd(l,m)
            end do    
            gso(m) = real( dgsofac(spin,h(m)) )     
            if (m .eq. 2) xi_lp(i,m) = eta_0*xi_lp(i,m)
            !Calculate the incident angle for this bin
            kk = get_index(rlp_column, ndelta, rad(i), rmin, npts(m))
            cosd_interp = interper(rlp_column, cosd_column, ndelta, rad(i), kk)
            if( kk .eq. npts(m) ) cosd_interp = newtex(rlp_column, cosd_column, ndelta, rad(i), h(m), honr, kk)
            xi_lp(i,m) = xiraw(rad(i),spin,h(m),honr,rlp_column,dcosdr_column,ndelta,rmin,npts(m),mudisk,cosd_interp, gsd(m))            
            mui = dinang(spin, rad(i), h(m), cosd_interp)
            !Correction to account for the radial dependence of incident angle, and for the g factors
            xi_lp(i,m) = xi_lp(i,m)/(sqrt(2.)*mui)*contx_int(m)*(gso(m))**(Gamma-2)  
            xitot = xitot + xi_lp(i,m)
            gsd_temp = gsd_temp + gsd(m)*xi_lp(i,m)
        end do 
        !This and the line above calculate the gsd factor along the disk, averaging over the flux the disk sees from each LP 
        gsdr(i) = gsd_temp/xitot
        logxir(i) = log10(xitot) - logner(i)
        ! write(*,*) "DEBUG: ", rad(i), xitot,  logxir(i), logner(i)
     end do

     !ionisation value at ISCO to renormalise the ionisation profile
     ! do m=1,nlp
     !    !Calculate the incident angle for this bin
     !    kk = get_index(rlp_column, ndelta, risco, rmin, npts(m))
     !    cosd_interp = interper(rlp_column, cosd_column, ndelta, risco, kk)
     !    if( kk .eq. npts(m) ) cosd_interp = newtex(rlp_column, cosd_column, ndelta, risco, h(m), honr, kk)
     !    xi_isco = xiraw(risco,spin,h(m),honr,rlp_column,dcosdr_column,ndelta,rmin,npts(m),mudisk,cosd_interp, gsd(m))
     !    mui_isco = dinang(spin, risco, h(m), cosd_interp)
     !    !Correction to account for the radial dependence of incident angle, and for the g factors
     !    xi_isco = xi_isco/(sqrt(2.)*mui_isco)*contx_int(m)*(gso(m))**(Gamma-2)
     !    xitot_isco = xitot_isco + xi_isco
     ! enddo
     xitot = 0. 
     gsd_temp = 0.
     rnorm = risco
     do m=1,nlp
        do l=1,ndelta
           rlp_column(l) = rlp(l,m)
           dcosdr_column(l) = dcosdr(l,m)
           cosd_column(l) = cosd(l,m)
        end do
        gso(m) = real( dgsofac(spin,h(m)) )     
        ! if (m .eq. 2) xi_lp(i,m) = eta_0*xi_lp(i,m)
        !Calculate the incident angle for this bin
        kk = get_index(rlp_column, ndelta, rnorm, rmin, npts(m))
        cosd_interp = interper(rlp_column, cosd_column, ndelta, rnorm, kk)
        if( kk .eq. npts(m) ) cosd_interp = newtex(rlp_column, cosd_column, ndelta, rnorm, h(m), honr, kk)
        xi_isco = xiraw(rnorm,spin,h(m),honr,rlp_column,dcosdr_column,ndelta,rmin,npts(m),mudisk,cosd_interp, gsd(m))            
        mui = dinang(spin, rnorm, h(m), cosd_interp)
        !Correction to account for the radial dependence of incident angle, and for the g factors
        xi_isco = xi_isco/(sqrt(2.)*mui)*contx_int(m)*(gso(m))**(Gamma-2)  
        ! xitot = xitot + xi_isco
        ! gsd_temp = gsd_temp + gsd(m)*xi_isco
     end do
     !This and the line above calculate the gsd factor along the disk, averaging over the flux the disk sees from each LP 
     ! gsdr(i) = gsd_temp/xitot
     logxir_isco = log10(xi_isco) - ( adensity *  mylogne_zoneB(rnorm, rin))
    ! write(*,*) "DEBUG: xi_isco, logxir_isco, mylogne_zoneB(rnorm, rin) ", rnorm, xi_isco, logxir_isco, mylogne_zoneB(rnorm, rin)
     
     logxinorm = logxir_isco
     lognenorm = adensity * mylogne_zoneB(rnorm, rin)
    ! write(*,*) "DEBUG: logxinorm lognenorm ", risco, logxinorm, lognenorm


     ! if (adensity .lt. 2) then 
     !    logxinorm = maxval(logxir)
     !    lognenorm = minval(logner)
     ! else 
     !    logxinorm = log10(xi_isco) - mylogne_zoneB(risco, rin)
     !    lognenorm = mylogne_zoneB(risco, rin)
     ! endif
     
    write(*,*) "DEBUG: logxinorm lognenorm ", risco, logxinorm, lognenorm
    write(*,*) "DEBUG: -----------------"
    !After the loop calculate the max and the min - ionization renormalized wrt to the first LP    
    ! logxinorm = maxval(logxir)
    ! lognenorm = minval(logner)
    
    logxir = logxir - (logxinorm - logxip) 
    logner = logner - (lognenorm - lognep)    
    
    ! write(*,*) "DEBUG: -----------------"
    ! open(99, file='delete_density_profiles.qdp')
    ! write(99,*) "skip on"
    !     do i=1,xe
    !         write(99,*) rad(i), logxir(i), logner(i)
    !         write(*,*) "DEBUG: ", rad(i), logxir(i), logner(i)
    !     end do
    !  write(99,*) "log x on"
    !  close(99)
    !  write(*,*) "DEBUG: -----------------"

    do m=1,nlp 
        do i=1,xe
            logxi_lp(i,m) = log10(xi_lp(i,m)) - logner(i) - lognenorm - logxinorm + lognep + logxip        
        end do
        logxip_lp(m) = max(maxval(logxi_lp(:,m)),0.)
    end do
    
    !Write radii, ionisation (for both and each LP), gamma factors, and log(xi(r))+log(ne(r)) (which is nearly the same as
    !epsilon(r) for identical coronal spectrra and gamma=2) to file. 
    !note 1) we need to do this before the ionisation array is set to have a minimum of 0, in order
    !to recover the correct scaling of the emissivity at large radii
    !2) in order to correctly compare the dfer_arr array with the single LP case, it has to be renormalized by (1+eta_0)
    if( verbose .gt. 1 ) then
        print*, "Peak ionisations from each LP: first " , logxip_lp(1), " second ", logxip_lp(2)
        open (unit = 27, file = 'Output/RadialScalings.dat', status='replace', action = 'write')
            do i = 1, xe
                write(27,*) rad(i), logxir(i), gsdr(i), logxir(i)+logner(i), logxi_lp(i,1), logxi_lp(i,2), dfer_arr(i) 
            end do 
        close(27)    
    end if
    
    !check max and min for ionisation 
    logxir = max( logxir , 0.d0  )
    logxir = min( logxir , 4.7d0 )
    logner = max( logner , 15.d0  )
    logner = min( logner , 20.d0 )
    
    deallocate(rad)

    return
end subroutine radfunctions_dens
!-----------------------------------------------------------------------
