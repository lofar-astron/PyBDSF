C$TEST  DTTGU5
c  main program5
      common /cstak/ ds
      double precision ds(350000)
      external handlu, bc, af
      integer ndx, ndy, istkgt, is(1000), iu, ix, temp, temp1
      integer nu, nr, iyb(5), ixb(5), kx, ky
      integer nxr(5), nyr(5), kxr(5), kyr(5)
      integer idumb
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision tstart, dt, rx
      double precision ws(500), tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to solve laplaces equation with real ( z*log(z) ) as solution.
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 4)
      call enter(1)
      nu = 1
      kx = 4
      ky = 4
      ndx = 3
      ndy = 3
      nr = 5
      tstart = 0
      dt = 1.d0
      tstop =1.d0
      errpar(1) = 1e-2
      errpar(2) = 1e-4
      nx = ndx+2*(kx-1)
      rx=1.0d0
c space for x mesh for rectangle 1
      ix = istkgt(nx, 4)
c 0 and rx mult = kx.
      ixb(1)=ix
      do  1 i = 1, kx
         temp = ix+i
         ws(temp-1) = 0
         temp = ix+nx-i
         ws(temp) = rx
   1     continue
      temp = ndx-1
      do  2 i = 1, temp
         temp1 = ix+kx-2+i
         ws(temp1) = rx*(dble(float(i-1))/(dble(float(ndx))-1d0))**kx
   2     continue
c rectangle 2 has same grid in x direction as rectanlge 1
      ixb(2)=istkgt(nx, 4)
      call dcopy(nx, ws(ix), 1, ws(ixb(2)), 1)
c uniform grid for rectanlges 3,4, and 5 in x direction
      ixb(3) = idumb(1.0d0, 2.0d0, ndx, kx, nxr(3))
      ixb(4) = idumb(2.0d0, 3.0d0, ndx, kx, nxr(4))
      ixb(5) = idumb(2.0d0, 3.0d0, ndx, kx, nxr(5))
      ny = ndy+2*(ky-1)
c rectangles 1,3, and 4 use the same grid in the y direction as
c is used for the x direction in rectangle 1
c space for y mesh.
      iyb(1) = istkgt(ny, 4)
      call dcopy( nx, ws(ix), 1, ws(iyb(1)), 1)
      iyb(3) =istkgt(ny, 4)
      call dcopy( nx, ws(ix), 1, ws(iyb(3)), 1)
      iyb(4) =istkgt(ny, 4)
      call dcopy( nx, ws(ix), 1, ws(iyb(4)), 1)
c rectangles 2 and 5 use uniform mesh in y direction
      iyb(2) = idumb(1.0d0, 2.0d0, ndy, ky, nyr(2))
      iyb(5) = idumb(1.0d0, 2.0d0, ndy, ky, nyr(5))
c space for the solution.
      nnu=0
      do 3 i=1,nr
         nxr(i)=nx
         nyr(i)=ny
         nnu=nnu+nu*((nxr(i)-kx)*(nyr(i)-ky))
 3    continue
      iu = istkgt(nnu, 4)
      do 4 i=1,nr
         kxr(i)=kx
         kyr(i)=ky
 4    continue
      call setd(nnu, 0.0d0,ws(iu))
      call dttgu(ws(iu),nu,nr,kxr,ws,nxr,ixb,kyr,ws,nyr,iyb,tstart,
     1   tstop, dt, af, bc, errpar, handlu)
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
      subroutine ewe2(t, x, nx, y, ny, u, inu, ir)
      integer inu,ir, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny)
      integer i, j
      double precision r, dcos, dlog, dsin, datan, theta
      double precision dsqrt
c the exact solution.
         do  6 i = 1, nx
            do  5 j = 1, ny
               r = dsqrt(x(i)**2+y(j)**2)
               if (x(i) .le. 0d0) goto 1
                  theta = datan(y(j)/x(i))
                  goto  2
   1              theta = 2d0*datan(1d0)
   2           if (r .le. 0d0) goto 3
                  u(i, j) = r*(dcos(theta)*dlog(r)-theta*dsin(theta))
                  goto  4
   3              u(i, j) = 0
   4           continue
   5           continue
   6        continue
      return
      end
      subroutine handlu(t0, u0, t, u, nv, dt, tstop)
      integer nv
      double precision t0, u0(nv), t, u(nv), dt, tstop
      common /d7tgup/ errpar, nu, mxp, myp
      integer nu
      real errpar(2)
      common /d7tgum/  kxp,ix,nxp,kyp,iy,nyp,nxnyt,nr,iup
      integer kx, ix, nx, ky, iy, ny
      common /cstak/is
      integer is(1000)
      iwrite=i1mach(2)
      if (t0 .ne. t) goto 2
         write (iwrite,  1) t
   1     format (16h restart for t =, 1pe10.2)
         return
c get and print the error.
   2    continue
        write(iwrite, 3)t
   3    format(6h at t=,1pe10.2)
        ius=1
        do 5 inu = 1, nu
           iyr=iy
           ixr=ix
           do 4 ir=1,nr
              ir1=ir-1
              nx=is(nxp+ir1)
              ny=is(nyp+ir1)
              kx=is(kxp+ir1)
              ky=is(kyp+ir1)
              call gerr(kx, ixr, nx, ky, iyr, ny, u(ius), inu, t, ir)
              ixr=ixr+nx
              iyr=iyr+ny
              ius=ius+(nx-kx)*(ny-ky)
   4        continue
   5  continue
      return
      end
      subroutine gerr(kx, ix, nx, ky, iy, ny, u, inu, t, ir)
      integer kx, ix, nx, ky, iy, ny, inu, ir
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
c for variable inu for rectangle ir
c u(nx-kx,ny-ky).
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
      call ewe2(t, ws(ixs), nxs, ws(iys), nys, ws(iewe), inu, ir)
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
      write (temp,  2) ir, inu, erru
   2  format(9h for rect,i3,14h error in u(.,,  i2,
     1       3h) =, 1pe10.2)
      call leave
      return
      end
