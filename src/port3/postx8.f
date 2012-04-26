c  main program
      common /cstak/ ds
      double precision ds(5000)
      common /time/ t
      real t
      common /kmesh/ k, nmesh
      integer k, nmesh
      common /cmesh/ mesh
      real mesh(100)
      external dee, handle, uofx, bc, af
      integer ndx, i, is(1000), nu, nv
      real errpar(2), u(100), v(100), dt, rs(1000), ws(1000)
      real tstop
      logical ls(1000)
      complex cs(500)
      integer temp
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test  post on the integro-pde
c      u sub t = 2 * u sub xx - int(0,1) exp(x-y)*u(y) dy      on (0,1)
c subject to given dirichlet bcs, chosen so that the solution is
c      u(x,t) = exp(t+x).
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(5000, 4)
      nu = 1
      errpar(1) = 0
c absolute error.
      errpar(2) = 1e-2
      tstop = 1
      dt = 1e-2
      k = 4
c ndx uniform mesh points on (0,1).
      ndx = 7
      call umb(0e0, 1e0, ndx, k, mesh, nmesh)
      nv = nmesh-k
c uofx needs t.
      t = 0
c ics for u.
      call l2sff(uofx, k, mesh, nmesh, u)
      temp = nmesh-k
      do  1 i = 1, temp
         v(i) = u(i)
   1     continue
c ics for v.
      call post(u, nu, k, mesh, nmesh, v, nv, 0e0, tstop, dt, af, bc, 
     1   dee, errpar, handle)
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
      common /kmesh/ k, nmesh
      integer k, nmesh
      common /cmesh/ mesh
      real mesh(100)
      integer i
      do  1 i = 1, nx
         a(i, 1) = 2.*ux(i, 1)
         aux(i, 1, 1) = 2
         f(i, 1) = ut(i, 1)
         fut(i, 1, 1) = 1
   1     continue
c get the integral.
      call intgrl(k, mesh, nmesh, v, x, nx, f, fv)
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
      real exp
      b(1, 1) = u(1, 1)-exp(t)
      b(1, 2) = u(1, 2)-exp(t+1.)
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
      integer i
      do  1 i = 1, nxmk
         d(i) = u(i, 1)-v(i)
         du(i, i, 1) = 1
         dv(i, i) = -1
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
      real eu, eesff
      integer temp
c output and checking routine.
      if (t0 .ne. t) goto 2
         temp = i1mach(2)
         write (temp,  1) t0, dt
   1     format (16h restart for t =, 1pe10.2, 7h   dt =, 1pe10.2)
         return
   2  tt = t
      eu = eesff(k, x, nx, u, uofx)
      temp = i1mach(2)
      write (temp,  3) t, eu
   3  format (14h error in u(x,, 1pe10.2, 4h ) =, 1pe10.2)
      return
      end
      subroutine uofx(x, nx, u, w)
      integer nx
      real x(nx), u(nx), w(nx)
      common /time/ t
      real t
      integer i
      real exp
      do  1 i = 1, nx
         u(i) = exp(t+x(i))
   1     continue
      return
      end
      subroutine intgrl(k, mesh, nmesh, v, x, nx, f, fv)
      integer nx, nmesh
      integer k
      real mesh(nmesh), v(1), x(nx), f(nx), fv(nx, 1)
      integer mgq, i, j, l, ix
      real ewe, ker, wgq(3), xgq(3), b(3, 4, 200), keru
      real xx(3)
      logical first
      integer temp, temp1
      data first/.true./
c to compute
c    f = integral from mesh(1) to mesh(nmesh)
c       kernel(x,y,sum(i=1,...,nmesh-k) v(i)*b(i,y)) dy
c  and
c    fv = d(f)/d(v).
c assume that call kernel(x,y,u,ker,keru) returns
c     ker = kernel(x,y,u) and
c     keru = partial kernel / partial u.
c v(nmesh-k),fv(nx,nmesh-k)
c the following declaration is specific to k = 4 splines.
      if (nmesh-k .gt. 200) call seterr(27hintgrl - nmesh-k .gt. nxmax
     1   , 27, 1, 2)
c need more local space.
      if (k .ne. 4) call seterr(17hintgrl - k .ne. 4, 17, 2, 2)
c use k-1 point gaussian-quadrature rule on each interval.
      mgq = k-1
      if (first) call gqm11(mgq, xgq, wgq)
c only get gq rule once, its expensive.
c the gaussian quadrature rule.
c do integral interval by interval.
      temp = nmesh-k
      do  6 i = k, temp
c g.q. points on (mesh(i), mesh(i+1)).
         do  1 j = 1, mgq
            xx(j) = 0.5*(mesh(i+1)+mesh(i))+0.5*(mesh(i+1)-mesh(i))*xgq(
     1         j)
   1        continue
         if (first) call bspln(k, mesh, nmesh, xx, mgq, i, b(1, 1, i))
c only get b-spline basis once, its expensive.
         do  5 j = 1, mgq
c get sum() v()*b()(xx).
            ewe = 0
            do  2 l = 1, k
               temp1 = i+l-k
               ewe = ewe+v(temp1)*b(j, l, i)
   2           continue
            do  4 ix = 1, nx
c get kernel and partial.
               call kernel(x(ix), xx(j), ewe, ker, keru)
               f(ix) = f(ix)+0.5*ker*(mesh(i+1)-mesh(i))*wgq(j)
               do  3 l = 1, k
                  temp1 = i+l-k
                  fv(ix, temp1) = fv(ix, temp1)+0.5*b(j, l, i)*keru*(
     1               mesh(i+1)-mesh(i))*wgq(j)
   3              continue
   4           continue
   5        continue
   6     continue
      first = .false.
      return
      end
      subroutine kernel(x, y, u, ker, keru)
      real x, y, u, ker, keru
      real exp
c to evaluate the kernel exp(x-y)*u(y) and its partial wrt. u.
      keru = exp(x-y)
      ker = keru*u
      return
      end
