!-----------------------------------------------------------------------
function mudiff(cosdelta,par)
    use blcoordinate
    !Calculates muobs (cosine of distant inclination angle) when given
    !cos(delta) (cosine of angle between initial photon trajectory and -z)
    implicit none
    double precision mudiff,cosdelta,par(3)
    double precision a_spin,h,muobs
    double precision scal,mus,sins
    double precision velocity(3),sindelta,pp,pr,pt,lambda,q,f1234(4)
    double precision ptotal,x,y,z,xprev,yprev,zprev,delx,dely,delz
    double precision ra,mua,phya,timea,sigmaa,p,cosdum
    a_spin = par(1)
    h      = par(2)
    muobs  = par(3)
    velocity = 0.0D0
    sindelta = sqrt( 1.d0 - cosdelta**2 )
    if ( sindelta .eq. 0.d0 )then
        cosdum = 1.d0
    else
        !Calculate 4-momentum in source rest frame tetrad
        pp= sindelta
        pr= cosdelta
        pt= 0.d0
        !Convert to LNRF (locally non-rotating reference frame)
        scal = 1.d0
        mus  = 1.d0
        sins = 0.d0
        call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
        ptotal = p_total(f1234(1),lambda,q,sins,mus,a_spin,h,scal)
        p = 0.9998d0 * ptotal
        call YNOGK(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
                   ra,mua,phya,timea,sigmaa)
        xprev = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*cos(phya)
        yprev = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*sin(phya)
        zprev = ra*mua
        p = 0.9999d0 * ptotal
        call YNOGK(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
                   ra,mua,phya,timea,sigmaa)
        x = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*cos(phya)
        y = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*sin(phya)
        z = ra*mua
        delx = x - xprev
        dely = y - yprev
        delz = z - zprev
        cosdum = delz / sqrt( delx**2 + dely**2 + delz**2 )
    end if
    
    mudiff = cosdum - muobs
    return
end function mudiff

!-----------------------------------------------------------------------
