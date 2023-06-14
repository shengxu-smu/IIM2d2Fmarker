c-----------------------------------------------------------------------
c
      subroutine fdf(a,b,fa,fb,ca,cb,c,f,df,num)
      integer num,n
      real*8 a,b,c,h
      real*8 fa(num),fb(num),ca(num),cb(num),f(num),df(num)
    
      h=b-a
      do n=1,num
        f(n)=(ca(n)*(b-c)**3.0d0+cb(n)*(c-a)**3.0d0)/(6.0d0*h)+
     .       (fa(n)*(b-c)+fb(n)*(c-a))/h-
     .        h*(ca(n)*(b-c)+cb(n)*(c-a))/6.0d0
        df(n)=(-ca(n)*(b-c)**2.0d0+cb(n)*(c-a)**2.0d0)/(2.0d0*h)+
     .        (fb(n)-fa(n))/h-h*(cb(n)-ca(n))/6.0d0
      enddo

      return
      end

      
c-----------------------------------------------------------------------
