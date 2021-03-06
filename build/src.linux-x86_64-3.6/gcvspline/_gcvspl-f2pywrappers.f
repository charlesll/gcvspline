C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapsplc (splcf2pywrap, m, n, k, y, ny, wx, w
     &y, mode, val, p, eps, c, nc, stat, b, we, el, bwe)
      external splc
      integer m
      integer n
      integer k
      integer ny
      integer mode
      double precision val
      double precision p
      double precision eps
      integer nc
      double precision el
      double precision y(ny,k)
      double precision wx(n)
      double precision wy(k)
      double precision c(nc,k)
      double precision stat(6)
      double precision b(2 * m - 1,n)
      double precision we(2 * m + 1,n)
      double precision bwe(2 * m + 1,n)
      double precision splcf2pywrap, splc
      splcf2pywrap = splc(m, n, k, y, ny, wx, wy, mode, val, p, ep
     &s, c, nc, stat, b, we, el, bwe)
      end


      subroutine f2pywraptrinv (trinvf2pywrap, b, e, m, n)
      external trinv
      integer m
      integer n
      double precision b(2 * m + 1,n)
      double precision e(2 * m + 1,n)
      double precision trinvf2pywrap, trinv
      trinvf2pywrap = trinv(b, e, m, n)
      end


      subroutine f2pywrapsplder (splderf2pywrap, ider, m, n, t, x,
     & c, l, q)
      external splder
      integer ider
      integer m
      integer n
      double precision t
      integer l
      double precision x(n)
      double precision c(n)
      double precision q(2 * m)
      double precision splderf2pywrap, splder
      splderf2pywrap = splder(ider, m, n, t, x, c, l, q)
      end

