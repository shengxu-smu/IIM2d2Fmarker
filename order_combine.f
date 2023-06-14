c-----------------------------------------------------------------------
c
      subroutine order_combine(a,na,la,b,nb,lb,ab,fa,fb,fab,num,nsame)
      integer na,la,nb,lb,num,nsame,nc
      integer i,m,ipa,ipb,ipab,ip
      real*8 abc,diff
      real*8 a(0:la),b(0:lb),ab(0:(la+lb))
      real*8 fa(num,0:la),fb(num,0:lb),fab(num,0:(la+lb))

      ipa=0
      ipb=0
      ipab=0
      nsame=0
5     continue
      if(a(ipa).lt.b(ipb)) then
        ab(ipab)=a(ipa)
        do m=1,num
          fab(m,ipab)=fa(m,ipa)
        enddo
        ipa=ipa+1
        ipab=ipab+1
      elseif(a(ipa).gt.b(ipb)) then
        ab(ipab)=b(ipb)
        do m=1,num
          fab(m,ipab)=fb(m,ipb)
        enddo
        ipb=ipb+1
        ipab=ipab+1
      elseif(abs(a(ipa)-b(ipb)).eq.0.0d0) then
        ab(ipab)=a(ipa)
        do m=1,num
          fab(m,ipab)=fa(m,ipa)
        enddo
        ipa=ipa+1
        ipb=ipb+1
        ipab=ipab+1
        nsame=nsame+1
      endif
      if(ipa.le.na.and.ipb.le.nb) goto 5
      if(ipa.eq.na+1) then
        do i=ipb,nb
          ab(ipab)=b(i)
          do m=1,num
            fab(m,ipab)=fb(m,i)
          enddo
          ipab=ipab+1
        enddo
      endif
      if(ipb.eq.nb+1) then
        do i=ipa,na
          ab(ipab)=a(i)
          do m=1,num
            fab(m,ipab)=fa(m,i)
          enddo
          ipab=ipab+1
        enddo
      endif

      diff=1.0d-3
      ip=1
      abc=ab(0)
      nc=na+nb+1-nsame
      do i=1,nc
        if((ab(i)-abc).lt.diff) then
          nsame=nsame+1
        else
          abc=ab(i)
          ab(ip)=ab(i)
          do m=1,num
            fab(m,ip)=fab(m,i)
          enddo
          ip=ip+1
        endif
      enddo

      return
      end


c-----------------------------------------------------------------------
