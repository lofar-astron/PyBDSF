C$TEST  DTTGR1P
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
c to solve the heat equation with solution u == t*x*y,
c   grad . ( u + ux + .1 * uy, u + uy + .1 * ux ) = ut + ux + uy +g(x,t)
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 4)
      call enter(1)
      nu = 1
      lx = 0
      rx = 1
      ly = 0
      ry = 1
      kx = 2
      ky = 2
      ndx = 3
      ndy = 3
      tstart = 0
      tstop = 1
      dt = 1
      errpar(1) = 1e-2
      errpar(2) = 1e-4
c uniform grid.
      ix = idumb(lx, rx, ndx, kx, nx)
c uniform grid.
      iy = idumb(ly, ry, ndy, ky, ny)
c space for the solution.
      iu = istkgt(nu*(nx-kx)*(ny-ky), 4)
c initial conditions for u.
      call setd(nu*(nx-kx)*(ny-ky), 0d0, ws(iu))
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
      integer i, p, q
      do  3 i = 1, nu
         do  2 q = 1, ny
            do  1 p = 1, nx
               a(p, q, i, 1) = ux(p, q, i)+.1*uy(p, q, i)+u(p, q, i)
               a(p, q, i, 2) = uy(p, q, i)+.1*ux(p, q, i)+u(p, q, i)
               aux(p, q, i, i, 1) = 1
               auy(p, q, i, i, 2) = 1
               auy(p, q, i, i, 1) = .1
               aux(p, q, i, i, 2) = .1
               au(p, q, i, i, 1) = 1
               au(p, q, i, i, 2) = 1
               f(p, q, i) = ut(p, q, i)+ux(p, q, i)+uy(p, q, i)
               fut(p, q, i, i) = 1
               fux(p, q, i, i) = 1
               fuy(p, q, i, i) = 1
               f(p, q, i) = f(p, q, i)+.2*t-x(p)*y(q)
   1           continue
   2        continue
   3     continue
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
      do  2 j = 1, ny
         do  1 i = 1, nx
            bu(i, j, 1, 1) = 1
            b(i, j, 1) = u(i, j, 1)-t*x(i)*y(j)
   1        continue
   2     continue
      return
      end
      subroutine handle(t0, u0, t, u, nv, dt, tstop)
      integer nv
      double precision t0, u0(nv), t, u(nv), dt, tstop
      common /d7tgrp/ errpar, nu, mxq, myq
      integer nu, mxq, myq
      real errpar(2)
      common /d7tgrm/ kx, ix, nx, ky, iy, ny
      integer kx, ix, nx, ky, iy, ny
      if (t0 .ne. t) goto 2
         write (6,  1) t
   1     format (16h restart for t =, 1pe10.2)
         return
c print results.
   2  call gerr(kx, ix, nx, ky, iy, ny, u, nu, t)
      return
      end
      subroutine gerr(kx, ix, nx, ky, iy, ny, u, nu, t)
      integer kx, ix, nx, ky, iy, ny
      integer nu
      double precision u(1), t
      common /cstak/ ds
      double precision ds(500)
      integer ifa, ita(2), ixa(2), nta(2), nxa(2), idlumd
      integer ixs, iys, nxs, nys, istkgt, i
      integer ka(2), ma(2), is(1000), i1mach
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision ws(500)
      integer temp, temp1
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to print the solution at each time-step.
c u(nx-kx,ny,ky,nu).
c the port library stack and its aliases.
      call enter(1)
c find the solution at 2 * 2 points / mesh rectangle.
c x search grid.
      ixs = idlumd(ws(ix), nx, 2, nxs)
c y search grid.
      iys = idlumd(ws(iy), ny, 2, nys)
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
c evaluate them.
      call dtsd1(2, ka, ws, ita, nta, u, ws, ixa, nxa, ma, ws(ifa))
      temp1 = ifa+nxs*nys-1
      temp = i1mach(2)
      write (temp,  1) t, (ws(i), i = ifa, temp1)
   1  format (3h u(, 1pe10.2, 7h,.,.) =, (1p5e10.2/20x,1p4e10.2))
      call leave
      return
      end
      subroutine ewe(t, x, nx, y, ny, u, nu)
      integer nu, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny, nu)
      integer i, j, p
c the exact solution.
      do  3 p = 1, nu
         do  2 i = 1, nx
            do  1 j = 1, ny
               u(i, j, p) = t*x(i)*y(j)
   1           continue
   2        continue
   3     continue
      return
      end
