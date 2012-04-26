C$TEST  DPOST10
c  main program
      common /cstak/ ds
      double precision ds(2000)
      external handle, dpostd, bc, af
      integer ndx, nxh, i, k, is(1000), nu
      integer nv, nx, i1mach
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision deebsf, err, dabs, u(100), v(1), x(100)
      double precision dmax1, dt, ue(100), uh(100), xh(100), ws(500)
      double precision tstop
      integer temp
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to estimate x and t error as sum.
c      u sub t = u sub xx + f      on (0,1)
c where f is chosen so that the solution is
c      u(x,t) = exp(xt).
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(2000, 4)
      nu = 1
      nv = 0
      errpar(1) = 0
      errpar(2) = 1e-2
      k = 4
      ndx = 4
      tstop = 1
      dt = 1d-2
c crude mesh.
      call dumb(0d0, 1d0, ndx, k, x, nx)
c initial conditions for u.
      call setd(nx-k, 1d0, u)
      temp = i1mach(2)
      write (temp,  1) 
   1  format (36h solving on crude mesh using errpar.)
      call dpost(u, nu, k, x, nx, v, nv, 0d0, tstop, dt, af, bc, dpostd,
     1   errpar, handle)
c get run-time statistics.
      call dpostx
c halve the mesh spacing.
      call dumb(0d0, 1d0, 2*ndx-1, k, xh, nxh)
c initial conditions for uh.
      call setd(nxh-k, 1d0, uh)
      dt = 1d-2
      temp = i1mach(2)
      write (temp,  2) 
   2  format (38h solving on refined mesh using errpar.)
      call dpost(uh, nu, k, xh, nxh, v, nv, 0d0, tstop, dt, af, bc, 
     1   dpostd, errpar, handle)
c get run-time statistics.
      call dpostx
c estimate u error.
      err = deebsf(k, x, nx, u, xh, nxh, uh)
      write (6,  3) err
   3  format (24h u error from u and uh =, 1pe10.2)
c initial conditions for ue.
      call setd(nx-k, 1d0, ue)
      dt = 1d-2
      errpar(1) = errpar(1)/10.
      errpar(2) = errpar(2)/10.
      temp = i1mach(2)
      write (temp,  4) 
   4  format (39h solving on crude mesh using errpar/10.)
      call dpost(ue, nu, k, x, nx, v, nv, 0d0, tstop, dt, af, bc, 
     1   dpostd, errpar, handle)
c get run-time statistics.
      call dpostx
      err = 0
      temp = nx-k
      do  5 i = 1, temp
         err = dmax1(err, dabs(u(i)-ue(i)))
   5     continue
      write (6,  6) err
   6  format (24h u error from u and ue =, 1pe10.2)
      stop 
      end
      subroutine af(t, x, nx, u, ux, ut, utx, nu, v, vt, nv, a, 
     1   au, aux, aut, autx, av, avt, f, fu, fux, fut, futx, fv, fvt)
      integer nu, nx
      integer nv
      double precision t, x(nx), u(nx, nu), ux(nx, nu), ut(nx, nu), utx(
     1   nx, nu)
      double precision v(1), vt(1), a(nx, nu), au(nx, nu, nu), aux(nx, 
     1   nu, nu), aut(nx, nu, nu)
      double precision autx(nx, nu, nu), av(1), avt(1), f(nx, nu), fu(
     1   nx, nu, nu), fux(nx, nu, nu)
      double precision fut(nx, nu, nu), futx(nx, nu, nu), fv(1), fvt(1)
      integer i
      double precision dexp
      do  1 i = 1, nx
         a(i, 1) = -ux(i, 1)
         aux(i, 1, 1) = -1
         f(i, 1) = (x(i)-t**2)*dexp(x(i)*t)-ut(i, 1)
         fut(i, 1, 1) = -1
   1     continue
      return
      end
      subroutine bc(t, l, r, u, ux, ut, utx, nu, v, vt, nv, b, bu,
     1   bux, but, butx, bv, bvt)
      integer nu
      integer nv
      double precision t, l, r, u(nu, 2), ux(nu, 2), ut(nu, 2)
      double precision utx(nu, 2), v(1), vt(1), b(nu, 2), bu(nu, nu, 2),
     1   bux(nu, nu, 2)
      double precision but(nu, nu, 2), butx(nu, nu, 2), bv(1), bvt(1)
      double precision dexp
      b(1, 1) = u(1, 1)-1d0
      b(1, 2) = u(1, 2)-dexp(t)
      bu(1, 1, 1) = 1
      bu(1, 1, 2) = 1
      return
      end
      subroutine handle(t0, u0, v0, t, u, v, nu, nxmk, nv, k, x, 
     1   nx, dt, tstop)
      integer nxmk, nu, nx
      integer nv, k
      double precision t0, u0(nxmk, nu), v0(1), t, u(nxmk, nu), v(1)
      double precision x(nx), dt, tstop
      common /time/ tt
      double precision tt
      external uofx
      integer i1mach
      double precision deesff, eu
      integer temp
c output and checking routine.
      if (t0 .ne. t) goto 2
         temp = i1mach(2)
         write (temp,  1) t
   1     format (16h restart for t =, 1pe10.2)
         return
   2  tt = t
      eu = deesff(k, x, nx, u, uofx)
      temp = i1mach(2)
      write (temp,  3) t, eu
   3  format (14h error in u(x,, 1pe10.2, 4h ) =, 1pe10.2)
      return
      end
      subroutine uofx(x, nx, u, w)
      integer nx
      double precision x(nx), u(nx), w(nx)
      common /time/ t
      double precision t
      integer i
      double precision dexp
      do  1 i = 1, nx
         u(i) = dexp(x(i)*t)
   1     continue
      return
      end
