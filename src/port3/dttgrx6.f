C$TEST  DTTGR6
c  main program
      common /cstak/ ds
      double precision ds(350000)
      external handle, bc, af
      integer iue, ndx, ndy, iur, ixr, iyr
      integer nxr, nyr, istkgt, i, is(1000), iu
      integer ix, iy, nu, kx, nx, ky
      integer ny, i1mach
      real errpar(2), rs(1000), float
      logical ls(1000)
      complex cs(500)
      double precision tstart, dble, dabs, eerr, erre, errr
      double precision dmax1, dt, lx, ly, rx, ry
      double precision ws(500), tstop
      integer temp, temp1, temp2
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to get error estimates for laplaces equation with real ( z*log(z) ) as
c solution.
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 4)
      call enter(1)
      nu = 1
      lx = 0
      rx = 1
      ly = 0
      ry = 1
      kx = 4
      ky = 4
      ndx = 2
      ndy = 2
      tstart = 0
      tstop = 1
      dt = 1
      errpar(1) = 1e-2
      errpar(2) = 1e-4
      nx = ndx+2*(kx-1)
c space for x mesh.
      ix = istkgt(nx, 4)
      do  1 i = 1, kx
         temp = ix+i
         ws(temp-1) = 0
         temp = ix+nx-i
         ws(temp) = rx
   1     continue
c 0 and rx mult = kx.
      temp = ndx-1
      do  2 i = 1, temp
         temp2 = ix+kx-2+i
         ws(temp2) = rx*(dble(float(i-1))/(dble(float(ndx))-1d0))**kx
   2     continue
      ny = ndy+2*(ky-1)
c space for y mesh.
      iy = istkgt(ny, 4)
      do  3 i = 1, ky
         temp = iy+i
         ws(temp-1) = 0
         temp = iy+ny-i
         ws(temp) = ry
   3     continue
c 0 and ry mult = ky.
      temp = ndy-1
      do  4 i = 1, temp
         temp2 = iy+ky-2+i
         ws(temp2) = ry*(dble(float(i-1))/(dble(float(ndy))-1d0))**ky
   4     continue
c space for the solution.
      iu = istkgt(nu*(nx-kx)*(ny-ky), 4)
      call setd(nu*(nx-kx)*(ny-ky), 0d0, ws(iu))
      temp = i1mach(2)
      write (temp,  5) 
   5  format (23h solving on crude mesh.)
      call dttgr(ws(iu), nu, kx, ws(ix), nx, ky, ws(iy), ny, tstart, 
     1   tstop, dt, af, bc, errpar, handle)
      dt = 1
      ndx = 2*ndx-1
c refine mesh.
      ndy = 2*ndy-1
      nxr = ndx+2*(kx-1)
c space for x mesh.
      ixr = istkgt(nxr, 4)
      do  6 i = 1, kx
         temp = ixr+i
         ws(temp-1) = 0
         temp = ixr+nxr-i
         ws(temp) = rx
   6     continue
c 0 and rx mult = kx.
      temp = ndx-1
      do  7 i = 1, temp
         temp2 = ixr+kx-2+i
         ws(temp2) = rx*(dble(float(i-1))/(dble(float(ndx))-1d0))**kx
   7     continue
      nyr = ndy+2*(ky-1)
c space for y mesh.
      iyr = istkgt(nyr, 4)
      do  8 i = 1, ky
         temp = iyr+i
         ws(temp-1) = 0
         temp = iyr+nyr-i
         ws(temp) = ry
   8     continue
c 0 and ry mult = ky.
      temp = ndy-1
      do  9 i = 1, temp
         temp2 = iyr+ky-2+i
         ws(temp2) = ry*(dble(float(i-1))/(dble(float(ndy))-1d0))**ky
   9     continue
c space for the solution.
      iur = istkgt(nu*(nxr-kx)*(nyr-ky), 4)
      call setd(nu*(nxr-kx)*(nyr-ky), 0d0, ws(iur))
      temp = i1mach(2)
      write (temp,  10) 
  10  format (25h solving on refined mesh.)
      call dttgr(ws(iur), nu, kx, ws(ixr), nxr, ky, ws(iyr), nyr, 
     1   tstart, tstop, dt, af, bc, errpar, handle)
      dt = 1
      errpar(1) = errpar(1)/10.
      errpar(2) = errpar(2)/10.
c space for the solution.
      iue = istkgt(nu*(nx-kx)*(ny-ky), 4)
      call setd(nu*(nx-kx)*(ny-ky), 0d0, ws(iue))
      temp = i1mach(2)
      write (temp,  11) 
  11  format (24h solving with errpar/10.)
      call dttgr(ws(iue), nu, kx, ws(ix), nx, ky, ws(iy), ny, tstart, 
     1   tstop, dt, af, bc, errpar, handle)
      errr = eerr(kx, ix, nx, ky, iy, ny, ws(iu), nu, ixr, nxr, iyr, 
     1   nyr, ws(iur), tstop)
      erre = 0
      temp = nu*(nx-kx)*(ny-ky)
      do  12 i = 1, temp
         temp2 = iu+i
         temp1 = iue+i
         erre = dmax1(erre, dabs(ws(temp2-1)-ws(temp1-1)))
  12     continue
      temp = i1mach(2)
      write (temp,  13) erre
  13  format (24h u error from u and ue =, 1pe10.2)
      temp = i1mach(2)
      write (temp,  14) errr
  14  format (24h u error from u and ur =, 1pe10.2)
      call leave
      call wrapup
      stop 
      end
      subroutine af(t, xi, nx, yi, ny, nu, u, ut, ux, uy, uxt, 
     1   uyt, a, au, aut, aux, auy, auxt, auyt, f, fu, fut, fux, fuy, 
     2   fuxt, fuyt)
      integer nu, nx, ny
      double precision t, xi(nx), yi(ny), u(nx, ny, nu), ut(nx, ny, nu),
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
               a(p, q, i, 1) = ux(p, q, i)
               a(p, q, i, 2) = uy(p, q, i)
               aux(p, q, i, i, 1) = 1
               auy(p, q, i, i, 2) = 1
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
      double precision r, dcos, dlog, dsin, datan, theta
      double precision dsqrt
      do  6 j = 1, ny
         do  5 i = 1, nx
            if (y(j) .ne. ly) goto 1
               b(i, j, 1) = uy(i, j, 1)
c neumann data on bottom.
               buy(i, j, 1, 1) = 1
               goto  4
   1           r = dsqrt(x(i)**2+y(j)**2)
c dirichlet data.
               if (x(i) .le. 0d0) goto 2
                  theta = datan(y(j)/x(i))
                  goto  3
   2              theta = 2d0*datan(1d0)
   3           b(i, j, 1) = u(i, j, 1)-r*(dcos(theta)*dlog(r)-theta*
     1            dsin(theta))
               bu(i, j, 1, 1) = 1
   4        continue
   5        continue
   6     continue
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
      integer iewe, ka(2), ma(2), is(1000), i1mach
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision dabs, erru, dmax1, ws(500)
      integer temp, temp1, temp2
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to get and print the error at each time-step.
c u(nx-kx,ny-ky,nu).
c the port library stack and its aliases.
      call enter(1)
c find the error in the solution at 2*kx * 2*ky points / mesh rectangle.
c x search grid.
      ixs = idlumd(ws(ix), nx, 2*kx, nxs)
c y search grid.
      iys = idlumd(ws(iy), ny, 2*ky, nys)
c u search grid values.
      iewe = istkgt(nxs*nys, 4)
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
c evaluate them.
      call dtsd1(2, ka, ws, ita, nta, u, ws, ixa, nxa, ma, ws(ifa))
c error in solution values.
      erru = 0
      temp = nxs*nys
      do  1 i = 1, temp
         temp2 = iewe+i
         temp1 = ifa+i
         erru = dmax1(erru, dabs(ws(temp2-1)-ws(temp1-1)))
   1     continue
      temp = i1mach(2)
      write (temp,  2) t, erru
   2  format (14h error in u(.,, 1pe10.2, 3h) =, 1pe10.2)
      call leave
      return
      end
      double precision function eerr(kx, ix, nx, ky, iy, ny, u, 
     1   nu, ixr, nxr, iyr, nyr, ur, t)
      integer kx, ix, nx, ky, iy, ny
      integer nu, ixr, nxr, iyr, nyr
      double precision u(1), ur(1), t
      common /cstak/ ds
      double precision ds(500)
      integer ifa, ita(2), ixa(2), nta(2), nxa(2), idlumd
      integer ixs, iys, nxs, nys, istkgt, i
      integer ifar, ka(2), ma(2), is(1000), i1mach
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision dabs, erru, dmax1, ws(500)
      integer temp, temp1, temp2
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to get and print the error estimate at each time-step.
c u(nx-kx,ny-ky,nu), ur(nxr-kx,nyr-ky,nu).
c the port library stack and its aliases.
      call enter(1)
c find the error in the solution at 2*kx * 2*ky points / fine mesh recta
cngle.
c x search grid.
      ixs = idlumd(ws(ixr), nxr, 2*kx, nxs)
c y search grid.
      iys = idlumd(ws(iyr), nyr, 2*ky, nys)
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
      ka(1) = kx
      ka(2) = ky
      ita(1) = ixr
      ita(2) = iyr
      nta(1) = nxr
      nta(2) = nyr
      ixa(1) = ixs
      ixa(2) = iys
      nxa(1) = nxs
      nxa(2) = nys
      ma(1) = 0
c get solution.
      ma(2) = 0
c approximate solution values.
      ifar = istkgt(nxs*nys, 4)
c evaluate them.
      call dtsd1(2, ka, ws, ita, nta, ur, ws, ixa, nxa, ma, ws(ifar))
c error in solution values.
      erru = 0
      temp = nxs*nys
      do  1 i = 1, temp
         temp2 = ifar+i
         temp1 = ifa+i
         erru = dmax1(erru, dabs(ws(temp2-1)-ws(temp1-1)))
   1     continue
      call leave
      eerr = erru
      return
      end
      subroutine ewe(t, x, nx, y, ny, u, nu)
      integer nu, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny, nu)
      integer i, j, p
      double precision r, dcos, dlog, dsin, datan, theta
      double precision dsqrt
c the exact solution.
      do  7 p = 1, nu
         do  6 i = 1, nx
            do  5 j = 1, ny
               r = dsqrt(x(i)**2+y(j)**2)
               if (x(i) .le. 0d0) goto 1
                  theta = datan(y(j)/x(i))
                  goto  2
   1              theta = 2d0*datan(1d0)
   2           if (r .le. 0d0) goto 3
                  u(i, j, p) = r*(dcos(theta)*dlog(r)-theta*dsin(theta))
                  goto  4
   3              u(i, j, p) = 0
   4           continue
   5           continue
   6        continue
   7     continue
      return
      end
