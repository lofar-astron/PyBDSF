c  main program
      common /cstak/ ds
      double precision ds(2000)
      external dee, handle, bc, af
      integer ndx, k, is(1000), nu, nv, nmesh
      real errpar(2), u(100), v(1), atan, mesh(100), dt
      real rs(1000), ws(1000), tstop
      logical ls(1000)
      complex cs(500)
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test  post on
c      u sub t = u sub xx - u**3 + f      on (-pi,+pi)
c subject to periodic boundary conditions,
c where f is chosen so that the solution is
c      u(x,t) = cos(x)*sin(t).
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(2000, 4)
      nu = 1
      nv = 1
      errpar(1) = 0
c absolute error.
      errpar(2) = 1e-2
      tstop = 8.*atan(1e0)
      dt = 0.4
c make a mesh of ndx uniform points on (-pi,+pi).
      k = 4
      ndx = 7
      call umb((-4.)*atan(1e0), 4.*atan(1e0), ndx, k, mesh, nmesh)
c initial conditions for u.
      call setr(nmesh-k, 0e0, u)
c initial conditions for v.
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
         f(i, 1) = (-ut(i, 1))-u(i, 1)**3+cos(x(i))*(cos(t)+sin(t)+cos(x
     1      (i))**2*sin(t)**3)
         fut(i, 1, 1) = -1
         fu(i, 1, 1) = (-3.)*u(i, 1)**2
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
      b(1, 1) = ux(1, 1)-v(1)
      b(1, 2) = ux(1, 2)-v(1)
      bux(1, 1, 1) = 1
      bv(1, 1, 1) = -1
      bux(1, 1, 2) = 1
      bv(1, 1, 2) = -1
      return
      end
      subroutine dee(t, k, x, nx, u, ut, nu, nxmk, v, vt, nv, d, 
     1   du, dut, dv, dvt)
      integer nxmk, nu, nv, nx
      integer k
      real t, x(nx), u(nxmk, nu), ut(nxmk, nu), v(nv), vt(nv)
      real d(nv), du(nv, nxmk, nu), dut(nv, nxmk, nu), dv(nv, nv), dvt(
     1   nv, nv)
      integer temp
c u(-pi,t) - u(+pi,t) = 0.
      temp = nx-k
      d(1) = u(1, 1)-u(temp, 1)
      du(1, 1, 1) = 1
      temp = nx-k
      du(1, temp, 1) = -1
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
      real eu, eesff, ev
      integer temp
c output and checking routine.
      if (t0 .eq. t) return
c uofx needs time.
      tt = t
      eu = eesff(k, x, nx, u, uofx)
      ev = v(1)
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
      real cos, sin
      do  1 i = 1, nx
         u(i) = cos(x(i))*sin(t)
   1     continue
      return
      end
