!-----------------------------------------------------------------------
      function ave_weight2D(par1,par1_lo,par1_hi,par2,par2_lo,par2_hi&
           ,val_1lo_2lo,val_1lo_2hi,val_1hi_2lo,val_1hi_2hi)
        implicit none
        
        double precision :: par1,par1_lo,par1_hi,par2,par2_lo,par2_hi,val_1lo_2lo&
             ,val_1lo_2hi,val_1hi_2lo,val_1hi_2hi,temp1,temp2,ave_weight2D


        if (par2_hi .eq. par2_lo) then
           temp1 = val_1lo_2lo
           temp2 = val_1hi_2lo
        else 
           temp1 = (val_1lo_2hi - val_1lo_2lo) / (par2_hi - par2_lo) * &
                (par2 - par2_lo) + val_1lo_2lo

           temp2 = (val_1hi_2hi - val_1hi_2lo) / (par2_hi - par2_lo) * &
                (par2 - par2_lo) + val_1hi_2lo
        endif

        if (par1_hi .eq. par1_lo) then
           ave_weight2D = temp1
        else
           ave_weight2D = (temp2 - temp1) / (par1_hi - par1_lo) * &
                (par1 - par1_lo) + temp1
        endif
        
      end function ave_weight2D
!-----------------------------------------------------------------------
