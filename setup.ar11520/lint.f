      subroutine lint(xa,ya,n,x,y)
c-----linear interpolation at x
      implicit none
      real*8 xa(n),ya(n)
      real*8 x,y,h,a,b
      integer n, k, klo, khi

      klo=1
      khi=n

      if((x.lt.xa(1)).or.(x.gt.xa(n))) then
         write(6,63) x, xa(1), xa(n)
 63      format('x = ',e12.3,' x(1) = ',e12.3,' x(n) = ',e12.3)
         write(6,10)
 10      format(/,'BOGUS! x out of range in lint')
         stop
      endif

 5    if((khi-klo).gt.1) then
         k = (khi+klo)/2
         if(xa(k).gt.x) then
            khi = k
         else
            klo = k
         endif
         go to 5
      endif

      h = xa(khi) - xa(klo)
      if(h.eq.0.d0) then
         write(6,20)
 20      format(/,'BOGUS! xa not strictly increasing in splint')
         stop
      endif
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)

      return
      end
