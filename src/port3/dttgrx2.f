C$TEST  DTTGR2
c  main program
      common /cstak/ ds
      double precision ds(350000)
      external handle, bc, af
      integer ndx, ndy, istkgt, is(1000), iu, ix
      integer iy, nu, kx, nx, ky, ny
      integer idumb
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision tstart, dt, lx, ly, rx, ry
      double precision ws(500), tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to solve two coupled, nonlinear heat equations.
c   u1 sub t = div . ( u1x, u1y ) - u1*u2 + g1
c   u2 sub t = div . ( u2x, u2y ) - u1*u2 + g2
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 4)
      call enter(1)
      nu = 2
      lx = 0
      rx = 1
      ly = 0
      ry = 1
      kx = 4
      ky = 4
      ndx = 3
      ndy = 3
      tstart = 0
      tstop = 1
      dt = 1e-2
      errpar(1) = 1e-2
      errpar(2) = 1e-4
c uniform grid.
      ix = idumb(lx, rx, ndx, kx, nx)
c uniform grid.
      iy = idumb(ly, ry, ndy, ky, ny)
c space for the solution.
      iu = istkgt(nu*(nx-kx)*(ny-ky), 4)
      call setd(nu*(nx-kx)*(ny-ky), 1d0, ws(iu))
      call dttgr(ws(iu), nu, kx, ws(ix), nx, ky, ws(iy), ny, tstart, 
     1   tstop, dt, af, bc, errpar, handle)
      call leave
      call wrapup
      stop 
      end
      subroutine af(t, x, nx, y, ny, nu, u, ut, ux, uy, uxt, uyt
     1   , a, au, aut, aux, auy, auxt, auyt, f, fu, fut, fux, fuy, fuxt,
     2   fuyt)
      integer nu, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny, nu), ut(nx, ny, nu), 
     1   ux(nx, ny, nu)
      double precision uy(nx, ny, nu), uxt(nx, ny, nu), uyt(nx, ny, nu),
     1   a(nx, ny, nu, 2), au(nx, ny, nu, nu, 2), aut(nx, ny, nu, nu, 2)
      double precision aux(nx, ny, nu, nu, 2), auy(nx, ny, nu, nu, 2), 
     1   auxt(nx, ny, nu, nu, 2), auyt(nx, ny, nu, nu, 2), f(nx, ny, nu)
     2   , fu(nx, ny, nu, nu)
      double precision fut(nx, ny, nu, nu), fux(nx, ny, nu, nu), fuy(nx,
     1   ny, nu, nu), fuxt(nx, ny, nu, nu), fuyt(nx, ny, nu, nu)
      integer p, q
      double precision dexp
      do  2 q = 1, ny
         do  1 p = 1, nx
            a(p, q, 1, 1) = ux(p, q, 1)
            aux(p, q, 1, 1, 1) = 1
            a(p, q, 1, 2) = uy(p, q, 1)
            auy(p, q, 1, 1, 2) = 1
            f(p, q, 1) = ut(p, q, 1)+u(p, q, 1)*u(p, q, 2)
            fu(p, q, 1, 1) = u(p, q, 2)
            fu(p, q, 1, 2) = u(p, q, 1)
            fut(p, q, 1, 1) = 1
            a(p, q, 2, 1) = ux(p, q, 2)
            aux(p, q, 2, 2, 1) = 1
            a(p, q, 2, 2) = uy(p, q, 2)
            auy(p, q, 2, 2, 2) = 1
            f(p, q, 2) = ut(p, q, 2)+u(p, q, 1)*u(p, q, 2)
            fu(p, q, 2, 1) = u(p, q, 2)
            fu(p, q, 2, 2) = u(p, q, 1)
            fut(p, q, 2, 2) = 1
            f(p, q, 1) = f(p, q, 1)-(dexp(t*(x(p)-y(q)))*(x(p)-y(q)-2d0*
     1         t*t)+1d0)
            f(p, q, 2) = f(p, q, 2)-(dexp(t*(y(q)-x(p)))*(y(q)-x(p)-2d0*
     1         t*t)+1d0)
   1        continue
   2     continue
      return
      end
      subroutine bc(t, x, nx, y, ny, lx, rx, ly, ry, u, ut, ux, 
     1   uy, uxt, uyt, nu, b, bu, but, bux, buy, buxt, buyt)
      integer nu, nx, ny
      double precision t, x(nx), y(ny), lx, rx, ly
      double precision ry, u(nx, ny, nu), ut(nx, ny, nu), ux(nx, ny, nu)
     1   , uy(nx, ny, nu), uxt(nx, ny, nu)
      double precision uyt(nx, ny, nu), b(nx, ny, nu), bu(nx, ny, nu, 
     1   nu), but(nx, ny, nu, nu), bux(nx, ny, nu, nu), buy(nx, ny, nu
     2   , nu)
      double precision buxt(nx, ny, nu, nu), buyt(nx, ny, nu, nu)
      integer i, j
      double precision dexp
      do  2 j = 1, ny
         do  1 i = 1, nx
            bu(i, j, 1, 1) = 1
            b(i, j, 1) = u(i, j, 1)-dexp(t*(x(i)-y(j)))
            bu(i, j, 2, 2) = 1
            b(i, j, 2) = u(i, j, 2)-dexp(t*(y(j)-x(i)))
   1        continue
   2     continue
      return
      end
      subroutine handle(t0, u0, t, u, nv, dt, tstop)
      integer nv
      double precision t0, u0(nv), t, u(nv), dt, tstop
      common /cstak/ ds
      double precision ds(500)
      common /d7tgrp/ errpar, nu, mxq, myq
      integer nu, mxq, myq
      real errpar(2)
      common /d7tgrm/ kx, ix, nx, ky, iy, ny
      integer kx, ix, nx, ky, iy, ny
      integer ifa, ita(2), ixa(2), nta(2), nxa(2), idlumd
      integer ixs, iys, nxs, nys, istkgt, i
      integer j, iewe, ka(2), ma(2), is(1000), i1mach
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision dabs, erru, dmax1, ws(500)
      integer temp, temp1, temp2
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c the port library stack and its aliases.
      if (t0 .ne. t) goto 2
         write (6,  1) t
   1     format (16h restart for t =, 1pe10.2)
         return
   2  call enter(1)
c find the error in the solution at 2*kx * 2*ky points / mesh rectangle.
c x search grid.
      ixs = idlumd(ws(ix), nx, 2*kx, nxs)
c y search grid.
      iys = idlumd(ws(iy), ny, 2*ky, nys)
c u search grid values.
      iewe = istkgt(nu*nxs*nys, 4)
c the exact solution.
      call ewe(t, ws(ixs), nxs, ws(iys), nys, ws(iewe), nu)
      ka(1) = kx
      ka(2) = ky
      ita(1) = ix
      ita(2) = iy
      nta(1) = nx
      nta(2) = ny
      ixa(1) = ixs
      ixa(2) = iys
      nxa(1) = nxs
      nxa(2) = nys
      ma(1) = 0
c get solution.
      ma(2) = 0
c approximate solution values.
      ifa = istkgt(nxs*nys, 4)
      do  5 j = 1, nu
c evaluate them.
         temp = (j-1)*(nx-kx)*(ny-ky)
         call dtsd1(2, ka, ws, ita, nta, u(temp+1), ws, ixa, nxa, ma, 
     1      ws(ifa))
c error in solution values.
         erru = 0
         temp = nxs*nys
         do  3 i = 1, temp
            temp2 = iewe+i-1+(j-1)*nxs*nys
            temp1 = ifa+i
            erru = dmax1(erru, dabs(ws(temp2)-ws(temp1-1)))
   3        continue
         temp = i1mach(2)
         write (temp,  4) t, j, erru
   4     format (14h error in u(.,, 1pe10.2, 1h,, i2, 3h) =, 1pe10.2)
   5     continue
      call leave
      return
      end
      subroutine ewe(t, x, nx, y, ny, u, nu)
      integer nu, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny, nu)
      integer i, j, p
      real float
      double precision dble, dexp
c the exact solution.
      do  3 p = 1, nu
         do  2 i = 1, nx
            do  1 j = 1, ny
               u(i, j, p) = dexp(dble(float((-1)**(p+1)))*t*(x(i)-y(j)))
   1           continue
   2        continue
   3     continue
      return
      end
