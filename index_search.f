c-----------------------------------------------------------------------
c
      subroutine index_search(ki,kpiil,ko,ll,kpool,ks,ms,kl)
      integer k,ki,kpiil,ko,ll,kup,kdown
      integer kpool(0:ks,ms)

      ko=-1
      kup=-1
      kdown=-1

      do k=0,kl
        if(kpiil.eq.kpool(k,ll)) then
          kup=k
          exit
        endif
      enddo

      do k=kl,0,-1
        if(kpiil.eq.kpool(k,ll)) then
          kdown=k
          exit
        endif
      enddo

      if(kup.eq.-1.and.kdown.eq.-1) return

      if(abs(kup-ki).lt.abs(kdown-ki)) then
        ko=kup
      else
        ko=kdown
      endif

      return
      end


c-----------------------------------------------------------------------
