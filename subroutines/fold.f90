!-----------------------------------------------------------------------
      subroutine fold(NENMAX,CHNMAX,NENERG,photar,RESP,NGRP,FCHAN,LCHAN,S)
!     Fold the model around the response
      implicit none
      integer NENMAX,CHNMAX,NENERG,NGRP(NENMAX),FCHAN(NENMAX,CHNMAX)
      integer LCHAN(NENMAX,CHNMAX),I,J,K
      real photar(NENMAX),RESP(CHNMAX,NENMAX)
      real S(CHNMAX)
      S = 0.0
      do J = 1,NENERG
        do K = 1,NGRP(J)
          do I = FCHAN(J,K)+1,LCHAN(J,K)
            S(I) = S(I) + photar(J) * RESP(I,J)
          end do
        end do
      end do
      RETURN
      END
!-----------------------------------------------------------------------
