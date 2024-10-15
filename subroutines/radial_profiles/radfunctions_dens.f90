!-----------------------------------------------------------------------
subroutine radfunctions_dens(verbose,xe,rin,rnmax,eta_0,logxip,lognep,spin,h,Gamma,honr,rlp,dcosdr&
     &,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner,dfer_arr)
    ! In  : xe,rin,rnmax,eta_0,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
    ! logxir(xe),gsdr(xe)   logxi (ionization parameter) and gsd (source to disc blueshift) as a function of radius
    ! Out : logxir(1:xe), gsdr(1:xe), logner(1:xe)
    use env_variables
    implicit none
    integer         , intent(IN)   :: xe, ndelta, nlp, npts(nlp)
    double precision, intent(IN)   :: rin,rmin,rnmax,eta_0,logxip,lognep,spin,h(nlp),honr,Gamma,dfer_arr(xe)
    real                           :: gso(nlp)
    double precision, intent(IN)   :: rlp(ndelta,nlp), dcosdr(ndelta,nlp), cosd(ndelta,nlp), contx_int(nlp)
    double precision :: rlp_column(ndelta),dcosdr_column(ndelta),cosd_column(ndelta), dgsofac
    double precision, intent(INOUT):: logxir(xe), gsdr(xe), logner(xe)
    integer          :: i, kk, get_index, get_env_int, l, m, verbose
    double precision :: rp, logxinorm, lognenorm,  mus, interper, newtex, mui, dinang, gsd(nlp), dglpfacthick
    double precision :: xi_lp(xe,nlp), logxi_lp(xe,nlp), logxip_lp(nlp), xitot, xiraw, mylogne, mudisk, gsd_temp
    double precision, allocatable :: rad(:)

    ! Set disk opening angle
    mudisk   = honr / sqrt( honr**2 + 1.d0  )
    
    allocate(rad(xe))
    !Now calculate logxi itself
    ! The loop calculates the raw xi and raw n_e.
    ! This means they are without normalization: only to find the maximum and the minimum. Remember that the max of the ionisation is not the same as the minumim in the density because the flux depends on r
    !The loops calculates also the correction factor mui

    !TBD: include luminosity ratio between LPs 
    do i = 1, xe        
        rad(i) = (rnmax/rin)**(real(i-1) / real(xe))
        rad(i) = rad(i) + (rnmax/rin)**(real(i) / real(xe))
        rad(i) = rad(i) * rin * 0.5
        !Initialize total ionization tracker
        xitot = 0. 
        gsd_temp = 0.
        !Now calculate the raw density (this matters only for high dens model reltransD)
        logner(i) = adensity * mylogne(rad(i), rin)
        do m=1,nlp
            do l=1,ndelta
                rlp_column(l) = rlp(l,m)
                dcosdr_column(l) = dcosdr(l,m)
                cosd_column(l) = cosd(l,m)
            end do    
            gso(m) = real( dgsofac(spin,h(m)) )     
            xi_lp(i,m) = xiraw(rad(i),spin,h(m),honr,rlp_column,dcosdr_column,ndelta,rmin,npts(m),mudisk,gsd(m))            
            if (m .eq. 2) xi_lp(i,m) = eta_0*xi_lp(i,m)
            !Calculate the incident angle for this bin
            kk = get_index(rlp_column, ndelta, rad(i), rmin, npts(m))
            mus = interper(rlp_column, cosd_column, ndelta, rad(i), kk)
            if( kk .eq. npts(m) ) mus = newtex(rlp_column, cosd_column, ndelta, rad(i), h(m), honr, kk)
            mui = dinang(spin, rad(i), h(m), mus)
            !Correction to account for the radial dependence of incident angle, and for the g factors
            xi_lp(i,m) = xi_lp(i,m)/(sqrt(2.)*mui)*contx_int(m)*(gso(m))**(Gamma-2)  
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
    logxir = logxir - (logxinorm - logxip) 
    logner = logner - (lognenorm - lognep)    
    
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
    
    deallocate(rad)

    return
end subroutine radfunctions_dens
!-----------------------------------------------------------------------
