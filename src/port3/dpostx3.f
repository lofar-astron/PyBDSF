C$TEST  DPOST3
c  main program
      common /cstak/ ds
      double precision ds(2000)
      external dee, handle, bc, af
      integer ndx, k, is(1000), nu, nv, nmesh
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision u(100), v(1), mesh(100), dt, ws(500), tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test dpost on
c      u sub t = u sub xx + v + f      on (0,1)
c        v sub t = u( 1/2, t )
c where f is chosen so that the solution is
c      u(x,t) = cos(xt)   and    v(t) = 2 sin(t/2).
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(2000, 4)
      nu = 1
      nv = 1
      errpar(1) = 1e-2
c essentially relative error.
      errpar(2) = 1e-6
      tstop = 1
      dt = 1d-6
      k = 4
      ndx = 4
c ndx uniform mesh points on (0,1).
      call dumb(0d0, 1d0, ndx, k, mesh, nmesh)
c initial conditions for u.
      call setd(nmesh-k, 1d0, u)
c initial value for v.
      v(1) = 0
      call dpost(u, nu, k, mesh, nmesh, v, nv, 0d0, tstop, dt, af, bc, 
     1   dee, errpar, handle)
c check for errors and stack usage statistics.
      call wrapup
      stop 
      end
      subroutine af(t, x, nx, u, ux, ut, utx, nu, v, vt, nv, a, 
     1   au, aux, aut, autx, av, avt, f, fu, fux, fut, futx, fv, fvt)
      integer nu, nv, nx
      double precision t, x(nx), u(nx, nu), ux(nx, nu), ut(nx, nu), utx(
     1   nx, nu)
      double precision v(nv), vt(nv), a(nx, nu), au(nx, nu, nu), aux(nx,
     1   nu, nu), aut(nx, nu, nu)
      double precision autx(nx, nu, nu), av(nx, nu, nv), avt(nx, nu, nv)
     1   , f(nx, nu), fu(nx, nu, nu), fux(nx, nu, nu)
      double precision fut(nx, nu, nu), futx(nx, nu, nu), fv(nx, nu, nv)
     1   , fvt(nx, nu, nv)
      integer i
      double precision dcos, dsin
      do  1 i = 1, nx
         a(i, 1) = -ux(i, 1)
         aux(i, 1, 1) = -1
         f(i, 1) = v(1)-ut(i, 1)-x(i)*dsin(x(i)*t)+t**2*dcos(x(i)*t)-
     1      2d0*dsin(t/2d0)
         fut(i, 1, 1) = -1
         fv(i, 1, 1) = 1
   1     continue
      return
      end
      subroutine bc(t, l, r, u, ux, ut, utx, nu, v, vt, nv, b, bu,
     1   bux, but, butx, bv, bvt)
      integer nu, nv
      double precision t, l, r, u(nu, 2), ux(nu, 2), ut(nu, 2)
      double precision utx(nu, 2), v(nv), vt(nv), b(nu, 2), bu(nu, nu, 2
     1   ), bux(nu, nu, 2)
      double precision but(nu, nu, 2), butx(nu, nu, 2), bv(nu, nv, 2), 
     1   bvt(nu, nv, 2)
      double precision dcos
      b(1, 1) = u(1, 1)-1d0
      b(1, 2) = u(1, 2)-dcos(t)
      bu(1, 1, 1) = 1
      bu(1, 1, 2) = 1
      return
      end
      subroutine dee(t, k, x, nx, u, ut, nu, nxmk, v, vt, nv, d, 
     1   du, dut, dv, dvt)
      integer nxmk, nu, nv, nx
      integer k
      double precision t, x(nx), u(nxmk, nu), ut(nxmk, nu), v(nv), vt(
     1   nv)
      double precision d(nv), du(nv, nxmk, nu), dut(nv, nxmk, nu), dv(
     1   nv, nv), dvt(nv, nv)
      integer intrvd, i, ileft
      double precision xi(1), basis(10)
      integer temp
      xi(1) = 0.5d0
c find 0.5 in mesh.
      ileft = intrvd(nx, x, xi(1))
      if (k .gt. 10) call seterr(
     1   41hdee - k .gt. 10, need more space in basis, 41, 1, 2)
c b-spline basis at xi(1).
      call dbspln(k, x, nx, xi, 1, ileft, basis)
      d(1) = vt(1)
      dvt(1, 1) = 1
c vt(1) - u(0.5,t) = 0.
      do  1 i = 1, k
         temp = ileft+i-k
         d(1) = d(1)-u(temp, 1)*basis(i)
         temp = ileft+i-k
         du(1, temp, 1) = du(1, temp, 1)-basis(i)
   1     continue
      return
      end
      subroutine handle(t0, u0, v0, t, u, v, nu, nxmk, nv, k, x, 
     1   nx, dt, tstop)
      integer nxmk, nu, nv, nx
      integer k
      double precision t0, u0(nxmk, nu), v0(nv), t, u(nxmk, nu), v(nv)
      double precision x(nx), dt, tstop
      common /time/ tt
      double precision tt
      external uofx
      integer i1mach
      double precision deesff, dabs, dsin, eu, ev
      integer temp
c output and checking routine.
      if (t0 .eq. t) return
c uofx needs time.
      tt = t
      eu = deesff(k, x, nx, u, uofx)
      ev = dabs(v(1)-2d0*dsin(t/2d0))
      temp = i1mach(2)
      write (temp,  1) t, eu, ev
   1  format (14h error in u(x,, 1pe10.2, 4h ) =, 1pe10.2, 6h   v =, 1p
     1   e10.2)
      return
      end
      subroutine uofx(x, nx, u, w)
      integer nx
      double precision x(nx), u(nx), w(nx)
      common /time/ t
      double precision t
      integer i
      double precision dcos
      do  1 i = 1, nx
         u(i) = dcos(x(i)*t)
   1     continue
      return
      end
