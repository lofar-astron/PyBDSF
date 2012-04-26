c  main program
      common /cstak/ ds
      double precision ds(2000)
      common /param/ c
      real c
      external handle, bc, af, postd
      integer ndx, nxc, nxx, i, k, is(1000)
      integer nu, nv, nx, i1mach
      real ewe(1000), err, errpar(2), u(100), v(1), x(100)
      real errr, dt, xc(100), uc(100), eebsf, rs(1000)
      real ws(1000), xx(1000), tstop, r1mach
      logical ls(1000)
      complex cs(500)
      integer temp
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test  post on automatic, static mesh refinement.
c      u sub t = u sub xx + c * u sub x      on (0,1)
c the solution is
c      u(x,t) = exp(-c*x).
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(2000, 4)
      c = 50
      nu = 1
      nv = 0
      errpar(1) = 1e-1
      errpar(2) = 1e-1
      k = 4
      ndx = 8
      call umb(0e0, 1e0, ndx, k, xc, nxc)
c initial conditions for uc.
      call setr(nxc-k, 0e0, uc)
c infinity.
      err = r1mach(2)
   1  if (err .le. 1e-2) goto  6
c halve the crude x.
         call lumb(xc, nxc, 3, k, x, nx)
c fitting points for refinement.
         call lumd(x, nx, k, xx, nxx)
c uc on xx.
         call splne(k, xc, nxc, uc, xx, nxx, ewe)
c fit u to uc on mesh.
         call dl2sf(xx, ewe, nxx, k, x, nx, u)
         tstop = 1./r1mach(4)
         dt = 1e-6
         i = nx-2*(k-1)
         temp = i1mach(2)
         write (temp,  2) i
   2     format (18h solving for ndx =, i3)
         call post(u, nu, k, x, nx, v, nv, 0e0, tstop, dt, af, bc, 
     1      postd, errpar, handle)
c get run-time statistics.
         call postx
c error estimate for uc.
         err = eebsf(k, xc, nxc, uc, x, nx, u)
c error estimate for u.
         errr = err/16.
         temp = i1mach(2)
         write (temp,  3) err, errr
   3     format (21h error estimates uc =, 1pe10.2, 9h  and u =, 1p
     1      e10.2)
         nxc = nx
         do  4 i = 1, nx
            xc(i) = x(i)
   4        continue
         temp = nx-k
         do  5 i = 1, temp
            uc(i) = u(i)
   5        continue
         goto  1
   6  stop 
      end
      subroutine af(t, x, nx, u, ux, ut, utx, nu, v, vt, nv, a, 
     1   au, aux, aut, autx, av, avt, f, fu, fux, fut, futx, fv, fvt)
      integer nu, nx
      integer nv
      real t, x(nx), u(nx, nu), ux(nx, nu), ut(nx, nu), utx(nx, nu)
      real v(1), vt(1), a(nx, nu), au(nx, nu, nu), aux(nx, nu, nu), aut(
     1   nx, nu, nu)
      real autx(nx, nu, nu), av(1), avt(1), f(nx, nu), fu(nx, nu, nu), 
     1   fux(nx, nu, nu)
      real fut(nx, nu, nu), futx(nx, nu, nu), fv(1), fvt(1)
      common /param/ c
      real c
      integer i
      do  1 i = 1, nx
         a(i, 1) = ux(i, 1)+c*u(i, 1)
         aux(i, 1, 1) = 1
         au(i, 1, 1) = c
         f(i, 1) = ut(i, 1)
         fut(i, 1, 1) = 1
   1     continue
      return
      end
      subroutine bc(t, l, r, u, ux, ut, utx, nu, v, vt, nv, b, bu,
     1   bux, but, butx, bv, bvt)
      integer nu
      integer nv
      real t, l, r, u(nu, 2), ux(nu, 2), ut(nu, 2)
      real utx(nu, 2), v(1), vt(1), b(nu, 2), bu(nu, nu, 2), bux(nu, nu,
     1   2)
      real but(nu, nu, 2), butx(nu, nu, 2), bv(1), bvt(1)
      common /param/ c
      real c
      real exp
      b(1, 1) = u(1, 1)-1.
      b(1, 2) = u(1, 2)-exp(-c)
      bu(1, 1, 1) = 1
      bu(1, 1, 2) = 1
      return
      end
      subroutine handle(t0, u0, v0, t, u, v, nu, nxmk, nv, k, x, 
     1   nx, dt, tstop)
      integer nxmk, nu, nx
      integer nv, k
      real t0, u0(nxmk, nu), v0(1), t, u(nxmk, nu), v(1)
      real x(nx), dt, tstop
      common /time/ tt
      real tt
      external uofx
      integer i1mach
      real eu, eesff
      integer temp
c output and checking routine.
      if (t0 .ne. t) goto 2
         temp = i1mach(2)
         write (temp,  1) t
   1     format (16h restart for t =, 1pe10.2)
         return
   2  tt = t
      eu = eesff(k, x, nx, u, uofx)
      temp = i1mach(2)
      write (temp,  3) t, eu
   3  format (15h error in u(x, , 1pe10.2, 4h ) =, 1pe10.2)
      return
      end
      subroutine uofx(x, nx, u, w)
      integer nx
      real x(nx), u(nx), w(nx)
      common /param/ c
      real c
      common /time/ t
      real t
      integer i
      real exp
      do  1 i = 1, nx
         u(i) = exp((-c)*x(i))
   1     continue
      return
      end
