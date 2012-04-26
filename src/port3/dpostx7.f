C$TEST  DPOST7
c  main program
      common /cstak/ ds
      double precision ds(4000)
      common /time/ t
      double precision t
      common /param/ vc, x, xi0
      double precision vc(4), x(3), xi0
      external dee, handle, uofx, bc, af
      integer ndx, idlumb, istkgt, k, iu, is(1000)
      integer nu, nv, immmd, imesh, nmesh
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision tstart, d, v(4), dt, xb(3), ws(500)
      double precision tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to test dpost on
c      u sub t = u sub xx + f      on (20,10**6)
c where f is chosen so that the solution is
c      u(x,t) = exp(-x*t),
c and x(1,t) is chosen so that the boundary-layer is tracked
c implicitly by forcing u(x(1,t)/2.3/d,t) = 1/e.
c this is the same as requiring the exact solution to have
c u(x(1,t),t) = 10 ** -d.
c v(1,2,3) gives the moving mesh, v(4) is time.
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(4000, 4)
      call enter(1)
      nu = 1
      nv = 4
      errpar(1) = 1e-2
c mixed relative and absolute error.
      errpar(2) = 1e-2
      d = 3
c w(xi0,t) = 1/e.
      xi0 = 1./2.3/d
      tstart = 20
      tstop = 1d+6
      dt = 1d-2
      k = 4
c ndx uniform mesh points on each interval of xb.
      ndx = 6
      xb(1) = 0
      xb(2) = 1
      xb(3) = 2
c get mesh on port stack.
      imesh = idlumb(xb, 3, ndx, k, nmesh)
c make 1d0 of multiplicity k-1.
      imesh = immmd(imesh, nmesh, 1d0, k-1)
      x(1) = 0
      x(2) = 2.3*d/tstart
      x(3) = 1
c initial values for v.
      call dlplmg(3, x, vc)
c get u on port stack.
      iu = istkgt(nmesh-k, 4)
c uofx needs time.
      t = tstart
      vc(4) = tstart
c uofx needs v for mapping.
      call movefd(nv, vc, v)
c initial conditions for u.
      call dl2sff(uofx, k, ws(imesh), nmesh, ws(iu))
c output ics.
      call handle(t-1d0, ws(iu), v, t, ws(iu), v, nu, nmesh-k, nv, k, 
     1   ws(imesh), nmesh, dt, tstop)
      call dpost(ws(iu), nu, k, ws(imesh), nmesh, v, nv, tstart, tstop
     1   , dt, af, bc, dee, errpar, handle)
      call leave
      call wrapup
      stop 
      end
      subroutine af(t, xi, nx, u, ux, ut, utx, nu, v, vt, nv, a, 
     1   au, aux, aut, autx, av, avt, f, fu, fux, fut, futx, fv, fvt)
      integer nu, nv, nx
      double precision t, xi(nx), u(nx, nu), ux(nx, nu), ut(nx, nu), 
     1   utx(nx, nu)
      double precision v(nv), vt(nv), a(nx, nu), au(nx, nu, nu), aux(nx,
     1   nu, nu), aut(nx, nu, nu)
      double precision autx(nx, nu, nu), av(nx, nu, nv), avt(nx, nu, nv)
     1   , f(nx, nu), fu(nx, nu, nu), fux(nx, nu, nu)
      double precision fut(nx, nu, nu), futx(nx, nu, nu), fv(nx, nu, nv)
     1   , fvt(nx, nu, nv)
      common /dpostf/ failed
      logical failed
      integer i
      double precision xxi(99), xtv(99), xvv(99), x(99), xxiv(99), ax(
     1   99)
      double precision fx(99), xt(99), xv(99), dexpl
      logical temp
      temp = v(2) .le. v(1)
      if (.not. temp) temp = v(2) .ge. v(3)
      if (.not. temp) goto 1
         failed = .true.
         return
c map xi into x.
   1  call dlplm(xi, nx, v, 3, x, xxi, xxiv, xv, xvv, xt, xtv)
c map u into x system.
      call dpostu(xi, x, xt, xxi, xv, vt, nx, 3, ux, ut, nu, ax, fx)
      do  2 i = 1, nx
         a(i, 1) = -ux(i, 1)
         aux(i, 1, 1) = -1
         f(i, 1) = (-ut(i, 1))-dexpl((-x(i))*v(4))*(x(i)+v(4)**2)
         fut(i, 1, 1) = -1
         fv(i, 1, 4) = (-dexpl((-x(i))*v(4)))*(2d0*v(4)+(x(i)+v(4)**2)*(
     1      -x(i)))
         fx(i) = (-dexpl((-x(i))*v(4)))*(1d0-v(4)*x(i)-v(4)**3)
   2     continue
c map a and f into xi system.
      call dposti(xi, x, xt, xxi, xv, xtv, xxiv, xvv, nx, ux, ut, nu, v,
     1   vt, nv, 1, 3, a, ax, au, aux, aut, autx, av, avt, f, fx, fu, 
     2   fux, fut, futx, fv, fvt)
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
      double precision dexpl
c u(0,t) = 1
      b(1, 1) = u(1, 1)-1d0
c u(1,t) = exp(-t)
      b(1, 2) = u(1, 2)-dexpl(-v(4))
      bu(1, 1, 1) = 1
      bu(1, 1, 2) = 1
      bv(1, 4, 2) = dexpl(-v(4))
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
      common /param/ vc, xc, xi0
      double precision vc(4), xc(3), xi0
      integer intrvd, i, ileft
      double precision dexp, bx(10), xx(1)
      integer temp
      d(1) = v(1)
c x(0,v) = 0.
      dv(1, 1) = 1
      xx(1) = xi0
      ileft = intrvd(nx, x, xx(1))
c get the b-spline basis at xx.
      call dbspln(k, x, nx, xx, 1, ileft, bx)
      d(2) = -dexp(-1d0)
c d(2) = w(xi0,t) - exp(-1).
      do  1 i = 1, k
         temp = ileft+i-k
         d(2) = d(2)+u(temp, 1)*bx(i)
         temp = ileft+i-k
         du(2, temp, 1) = bx(i)
   1     continue
      d(3) = v(3)-1d0
c x(2,v) = 1.
      dv(3, 3) = 1
      d(4) = vt(4)-1d0
      dvt(4, 4) = 1
      return
      end
      subroutine handle(t0, u0, v0, t, u, v, nu, nxmk, nv, k, x, 
     1   nx, dt, tstop)
      integer nxmk, nu, nv, nx
      integer k
      double precision t0, u0(nxmk, nu), v0(nv), t, u(nxmk, nu), v(nv)
      double precision x(nx), dt, tstop
      common /param/ vc, xx, xi0
      double precision vc(4), xx(3), xi0
      common /time/ tt
      double precision tt
      external uofx
      integer i1mach
      double precision deesff, dlplmt, eu, ev
      integer temp
c output and checking routine.
      if (t0 .ne. t) goto 2
         temp = i1mach(2)
         write (temp,  1) t, dt
   1     format (16h restart for t =, 1pe10.2, 7h   dt =, 1pe10.2)
         return
c let dt carry v(2) down by no more than a factor of 10.
   2  dt = dlplmt(t, v, nv, t0, v0, 1d-1, dt)
      tt = t
c uofx needs v for mapping.
      call movefd(nv, v, vc)
      eu = deesff(k, x, nx, u, uofx)
c error in position of boundary layer.
      ev = v(2)-1d0/xi0/t
      temp = i1mach(2)
      write (temp,  3) t, eu, ev
   3  format (14h error in u(x,, 1pe10.2, 4h ) =, 1pe10.2, 6h   v =, 1p
     1   e10.2)
      return
      end
      subroutine uofx(xi, nx, u, w)
      integer nx
      double precision xi(nx), u(nx), w(nx)
      common /cstak/ ds
      double precision ds(500)
      common /param/ vc, x, xi0
      double precision vc(4), x(3), xi0
      common /time/ t
      double precision t
      integer ixv, ixx, istkgt, i, is(1000)
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision ws(500), dexpl
      integer temp
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c the port library stack and its aliases.
      call enter(1)
      ixx = istkgt(nx, 4)
c space for x and xv.
      ixv = istkgt(3*nx, 4)
c map into user system.
      call dlplmx(xi, nx, vc, 3, ws(ixx), ws(ixv))
      do  1 i = 1, nx
         temp = ixx+i
         u(i) = dexpl((-ws(temp-1))*t)
   1     continue
      call leave
      return
      end
