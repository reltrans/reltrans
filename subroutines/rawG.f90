subroutine rawG(nex,earx,nf,flo,fhi,nlp,contx,absorbx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,h,z,Gamma,eta,&
                boost,ReIm,g,DelAB,ionvar,DC,resp_matr,ReGraw,ImGraw)
                
    use constants
    implicit none
    integer nex,nf,ionvar,DC,nlp
    complex W0,W1,W2,W3,Sraw(nlp,nex,nf),cexp_d,cexp_phi,Stemp   
    real earx(0:nex),absorbx(nex),contx(nex,nlp),tauso(nlp),ReW0(nlp,nex,nf),ImW0(nlp,nex,nf)
    real ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)
    real DelAB(nlp),g(nlp),boost,z,gso(nlp),Gamma,eta,ReSraw(nlp,nex,nf),ImSraw(nlp,nex,nf),h(nlp) 
    real tau_d,phase_d,E,fac,f,flo,fhi,ReGtemp(nlp,nex,nf),ImGtemp(nlp,nex,nf),ReGraw(nex,nf),ImGraw(nex,nf)
    integer i,j,m,ReIm,resp_matr
    
    ReSraw = 0.
    ImSraw = 0.
    ReGraw = 0.
    ImGraw = 0.
    Sraw = 0.
    
    phase_d = 0.
    tau_d = 0.

    do m=1,nlp 
        if (boost .lt. 0 .and. DC .eq. 1) then             
            do j = 1,nf
                do i = 1,nex
                    !if (m .gt. 1) ReW0(m,i,j) = eta*ReW0(m,i,j)
                    ReSraw(m,i,j) = ReSraw(m,i,j) + (-boost) * ReW0(m,i,j)
                enddo
            enddo  
        else
            do j = 1,nf
                if (DC .eq. 1) then
                    f = 0.
                else 
                    f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )
                endif
                do i = 1,nex
                    E   = 0.5 * ( earx(i) + earx(i-1) )
                    fac = log(gso(m)/((1.0+z)*E))
                    if (m .gt. 1) then
                        tau_d = tauso(m)-tauso(1)
                        phase_d = 2.*pi*tau_d*f  
                    endif
                    cexp_d = cmplx(cos(phase_d),sin(phase_d))     
                    cexp_phi = cmplx(cos(DelAB(m)),sin(DelAB(m)))
                    !set up transfer functions 
                    W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                    W1 = (1-DC) * boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                    W2 = (1-DC) * boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                    W3 = ionvar * (1-DC) * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                    !calculate complex covariance
                    !note: the reason we use complex here is to ease the calculations 
                    !when we add all the extra phases from the double lamp post 
                    Stemp = g(m)*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m)) + W0 + W3 + cexp_d*contx(i,m)
                    Sraw(m,i,j) = Sraw(m,i,j) + Stemp
                enddo 
            enddo    
        endif 
    end do

    !include absorption and separate - this is a bit awful but I think it's the only way?
    do m=1,nlp 
        do j=1,nf 
            do i=1,nex 
                ReSraw(m,i,j) = real(Sraw(m,i,j))*absorbx(i)
                ImSraw(m,i,j) = aimag(Sraw(m,i,j))*absorbx(i)
            end do
        end do
    end do

    do m=1,nlp 
        if (ReIm .gt. 0.0) then 
            call propercross(nex,nf,earx,ReSraw(m,:,:),ImSraw(m,:,:),ReGtemp(m,:,:),ImGtemp(m,:,:),resp_matr)
        else 
            call propercross_NOmatrix(nex,nf,earx,ReSraw(m,:,:),ImSraw(m,:,:),ReGtemp(m,:,:),ImGtemp(m,:,:))
        endif
        do j=1,nf 
            do i=1,nex 
                if (m .gt. 1) then
                    ReGtemp(m,i,j) = eta**2.*ReGtemp(m,i,j)
                    ImGtemp(m,i,j) = eta**2.*ImGtemp(m,i,j)
                endif
                ReGraw(i,j) = ReGraw(i,j) + ReGtemp(m,i,j)
                ImGraw(i,j) = ImGraw(i,j) + ImGtemp(m,i,j)
            end do
        end do
    end do   

    return 
end subroutine rawG
