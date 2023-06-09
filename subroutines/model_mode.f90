subroutine model_mode(ReIm,floHz,flo,fhiHz,fhi,Mass,nf,fc,DC,nlp,g,DelAB,DelA,eta_0,eta,beta_p,boost)
    use dyn_en
    implicit none
    integer, intent(in)     ::      nlp
    integer, intent(inout)  ::      ReIm    
    integer, intent(out)    ::      nf,DC
    integer                 ::      fbinx
    real, intent(in)        ::      Mass   
    real, intent(out)       ::      floHz,fhiHz,g(nlp),DelAB(nlp),DelA,beta_p,boost
    double precision, intent(out) :: flo,fhi,fc,eta_0,eta
    double precision, parameter ::  dlogf = 0.09 !This is a resolution parameter (base 10)      
    
    if( ReIm .eq. 7 .and. floHz .gt. tiny(floHz) .and. fhiHz .gt. tiny(fhiHz)) then
        !set up frequency grid in Hz if using lag frequency mode, depending on whether we're looking at AGN or XRBs
        if( Mass .gt. 1000 ) then
            floHz = 1.e-5
            fhiHz = 5.e-2 
        else
            floHz = 0.07
            fhiHz = 700.
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
        !Note that the frequency grid is using a higher resolution since it's what we care about in this mode 
        nf = ceiling( log10(fhiHz/floHz) / 0.01 )
        allocate(fix(0:nf))
        do fbinx = 0, nf 
            fix(fbinx) = floHz * (fhiHz / floHz)**( real(fbinx) / real(nf) )
        end do
    else 
        !if doing lag-energy spectra, just work out how many frequencies to average over 
        fc = 0.5d0 * ( floHz + fhiHz )
        nf = ceiling( log10(fhiHz/floHz) / dlogf )
        if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
            fhiHz = 0.d0
            floHz = 0.d0
            nf    = 1
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
    end if 
    
        !Decide if this is the DC component/time averaged spectrum or not
    if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
        DC     = 1
        g      = 0.0
        DelAB  = 0.0
        DelA   = 0.0
        ReIm   = 1
        eta    = eta_0
        beta_p = 1. !this is an ugly hack for the double LP model, to calculate the time-averaged spectrum
    else
        DC     = 0
        boost  = abs(boost)
    end if

end subroutine model_mode
