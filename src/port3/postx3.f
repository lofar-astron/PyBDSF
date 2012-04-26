c  main program
      common /cstak/ ds
      double precision ds(2000)
      external dee, handle, bc, af
      integer ndx, k, is(1000), nu, nv, nmesh
      real errpar(2), u(100), v(1), mesh(100), dt, rs(1000)
      real ws(1000), tstop
      logical ls(1000)
      complex cs(500)
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test  post on
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
      dt = 1e-6
      k = 4
      ndx = 4
c ndx uniform mesh points on (0,1).
      call umb(0e0, 1e0, ndx, k, mesh, nmesh)
c initial conditions for u.
      call setr(nmesh-k, 1e0, u)
c initial value for v.
      v(1) = 0
      call post(u, nu, k, mesh, nmesh, v, nv, 0e0, tstop, dt, af, bc, 
     1   dee, errpar, handle)
c check for errors and stack usage statistics.
      call wrapup
      stop 
      end
      subroutine af(t, x, nx, u, ux, ut, utx, nu, v, vt, nv, a, 
     1   au, aux, aut, autx, av, avt, f, fu, fux, fut, futx, fv, fvt)
      integer nu, nv, nx
      real t, x(nx), u(nx, nu), ux(nx, nu), ut(nx, nu), utx(nx, nu)
      real v(nv), vt(nv), a(nx, nu), au(nx, nu, nu), aux(nx, nu, nu), 
     1   aut(nx, nu, nu)
      real autx(nx, nu, nu), av(nx, nu, nv), avt(nx, nu, nv), f(nx, nu),
     1   fu(nx, nu, nu), fux(nx, nu, nu)
      real fut(nx, nu, nu), futx(nx, nu, nu), fv(nx, nu, nv), fvt(nx, 
     1   nu, nv)
      integer i
      real cos, sin
      do  1 i = 1, nx
         a(i, 1) = -ux(i, 1)
         aux(i, 1, 1) = -1
         f(i, 1) = v(1)-ut(i, 1)-x(i)*sin(x(i)*t)+t**2*cos(x(i)*t)-2.*
     1      sin(t/2.)
         fut(i, 1, 1) = -1
         fv(i, 1, 1) = 1
   1     continue
      return
      end
      subroutine bc(t, l, r, u, ux, ut, utx, nu, v, vt, nv, b, bu,
     1   bux, but, butx, bv, bvt)
      integer nu, nv
      real t, l, r, u(nu, 2), ux(nu, 2), ut(nu, 2)
      real utx(nu, 2), v(nv), vt(nv), b(nu, 2), bu(nu, nu, 2), bux(nu, 
     1   nu, 2)
      real but(nu, nu, 2), butx(nu, nu, 2), bv(nu, nv, 2), bvt(nu, nv, 2
     1   )
      real cos
      b(1, 1) = u(1, 1)-1.
      b(1, 2) = u(1, 2)-cos(t)
      bu(1, 1, 1) = 1
      bu(1, 1, 2) = 1
      return
      end
      subroutine dee(t, k, x, nx, u, ut, nu, nxmk, v, vt, nv, d, 
     1   du, dut, dv, dvt)
      integer nxmk, nu, nv, nx
      integer k
      real t, x(nx), u(nxmk, nu), ut(nxmk, nu), v(nv), vt(nv)
      real d(nv), du(nv, nxmk, nu), dut(nv, nxmk, nu), dv(nv, nv), dvt(
     1   nv, nv)
      integer intrvr, i, ileft
      real xi(1), basis(10)
      integer temp
      xi(1) = 0.5e0
c find 0.5 in mesh.
      ileft = intrvr(nx, x, xi(1))
      if (k .gt. 10) call seterr(
     1   41hdee - k .gt. 10, need more space in basis, 41, 1, 2)
c b-spline basis at xi(1).
      call bspln(k, x, nx, xi, 1, ileft, basis)
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
      real t0, u0(nxmk, nu), v0(nv), t, u(nxmk, nu), v(nv)
      real x(nx), dt, tstop
      common /time/ tt
      real tt
      external uofx
      integer i1mach
      real abs, sin, eu, ev, eesff
      integer temp
c output and checking routine.
      if (t0 .eq. t) return
c uofx needs time.
      tt = t
      eu = eesff(k, x, nx, u, uofx)
      ev = abs(v(1)-2.*sin(t/2.))
      temp = i1mach(2)
      write (temp,  1) t, eu, ev
   1  format (14h error in u(x,, 1pe10.2, 4h ) =, 1pe10.2, 6h   v =, 1p
     1   e10.2)
      return
      end
      subroutine uofx(x, nx, u, w)
      integer nx
      real x(nx), u(nx), w(nx)
      common /time/ t
      real t
      integer i
      real cos
      do  1 i = 1, nx
         u(i) = cos(x(i)*t)
   1     continue
      return
      end
