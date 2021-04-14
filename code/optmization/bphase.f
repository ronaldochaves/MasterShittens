c *************************************************** 
c *                                                 *
c *           BOUNDING PHASE METHOD                 *
c *                                                 *
c ***************************************************
c Developed by Dr. Kalyanmoy Deb
c Indian Institute of Technology, Kanpur
c All rights reserved. This listing is for personal use.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      Change the function funct() for a new function
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      implicit real*8 (a-h,o-z)
      nfun = 0
      call bphase(a,b,nfun)
      write(*,5) a,b,nfun
  5   format(2x,’The bounds are (’,f10.3,’, ’,f10.3,’)’,
     -          /,2x,’Total function evaluations:’,i6)
      stop
      end


      subroutine bphase(a,b,nfun)
c     bounding phase algorithm
c   ****************************************************
c   a and b are lower and upper bounds (output)
c   nfun is the number of function evaluations (output)
c   ****************************************************
      implicit real*8 (a-h,o-z)
c.....step 1 of the algorithm
  1   write(*,*) ’enter x0, delta’
      read(*,*) x0,delta
      call funct(x0-delta,fn,nfun)
      call funct(x0,f0,nfun)
      call funct(x0+delta,fp,nfun)
      write(*,*) ’enter 1 for intermediate results’
      read(*,*) iprint
c.....step 2 of the algorithm
      if (fn .ge. f0) then
        if (f0 .ge. fp) then
         delta = 1 * delta
        else
          a = x0 - delta
          b = x0 + delta
        endif
      elseif ((fn .le. f0) .and. (f0 .le. fp))  then
         delta = -1 * delta
      else
        go to 1
      endif
      k=0
      xn = x0 - delta
c.....step 3 of the algorithm
  3   x1 = x0 + (2**k) * delta
      call funct(x1,f1,nfun)
      if (iprint .eq. 1) then
        write(*,4) x1, f1
  4     format(2x,’Current point ’,f10.4,
     -            ’ function value ’,1pe15.4)
      endif
c.....step 4 of the algorithm
      if (f1 .lt. f0) then
        k = k+1
        xn = x0
        fn = f0
        x0 = x1
        f0 = f1
        go to 3
      else
        a = xn
        b = x1
      endif
      if (b .lt. a) then
        temp = a
        a = b
        b = temp
      endif
      return
      end

      subroutine funct(x,f,nfun)
c   ****************************************************
c   x is the current point (input)
c   f is the function value (output)
c   nfun if the number of function evaluations (output)
c   ****************************************************
      implicit real*8 (a-h,o-z)
      nfun = nfun + 1
      f = x*x + 54.0/x
      return
      end