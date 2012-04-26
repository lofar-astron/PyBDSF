C$TEST  TTGU5
c  main program
      common /cstak/ ds
      real ds(350000)
      external handlu, bc, af
      integer ndx, ndy, istkgt, is(1000), iu, ix, temp, temp1
      integer nu, nr, iyb(5), ixb(5), kx, ky
      integer nxr(5), nyr(5), kxr(5), kyr(5)
      integer  IUMB
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      real tstart, dt, rx
      real ws(500), tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to solve laplaces equation with real ( z*log(z) ) as solution.
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 3)
      call enter(1)
      nu = 1
      kx = 4
      ky = 4
      ndx = 3
      ndy = 3
      nr = 5
      tstart = 0
      dt = 1.e0
      tstop =1.e0
      errpar(1) = 1e-2
      errpar(2) = 1e-4
      nx = ndx+2*(kx-1)
      rx=1.0e0
c space for x mesh for rectangle 1
      ix = istkgt(nx, 3)
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
         ws(temp1) = rx*((float(i-1))/((float(ndx))-1e0))**kx
   2     continue
c rectangle 2 has same grid in x direction as rectanlge 1
      ixb(2)=istkgt(nx, 3)
      call  SCOPY(nx, ws(ix), 1, ws(ixb(2)), 1)
c uniform grid for rectanlges 3,4, and 5 in x direction
      ixb(3) =  IUMB(1.0e0, 2.0e0, ndx, kx, nxr(3))
      ixb(4) =  IUMB(2.0e0, 3.0e0, ndx, kx, nxr(4))
      ixb(5) =  IUMB(2.0e0, 3.0e0, ndx, kx, nxr(5))
      ny = ndy+2*(ky-1)
c rectangles 1,3, and 4 use the same grid in the y direction as
c is used for the x direction in rectangle 1
c space for y mesh.
      iyb(1) = istkgt(ny, 3)
      call  SCOPY( nx, ws(ix), 1, ws(iyb(1)), 1)
      iyb(3) =istkgt(ny, 3)
      call  SCOPY( nx, ws(ix), 1, ws(iyb(3)), 1)
      iyb(4) =istkgt(ny, 3)
      call  SCOPY( nx, ws(ix), 1, ws(iyb(4)), 1)
c rectangles 2 and 5 use uniform mesh in y direction
      iyb(2) =  IUMB(1.0e0, 2.0e0, ndy, ky, nyr(2))
      iyb(5) =  IUMB(1.0e0, 2.0e0, ndy, ky, nyr(5))
c space for the solution.
      nnu=0
      do 3 i=1,nr
         nxr(i)=nx
         nyr(i)=ny
         nnu=nnu+nu*((nxr(i)-kx)*(nyr(i)-ky))
 3    continue
      iu = istkgt(nnu, 3)
      do 4 i=1,nr
         kxr(i)=kx
         kyr(i)=ky
 4    continue
      call SETR(nnu, 0.0e0,ws(iu))
      call ttgu(ws(iu),nu,nr,kxr,ws,nxr,ixb,kyr,ws,nyr,iyb,tstart,
     1   tstop, dt, af, bc, errpar, handlu)
      call leave
      call wrapup
      stop 
      end
      subroutine af(t, xi, nx, yi, ny, nu, u, ut, ux, uy, uxt, 
     1   uyt, a, au, aut, aux, auy, auxt, auyt, f, fu, fut, fux, fuy, 
     2   fuxt, fuyt)
      integer nu, nx, ny
      real t, xi(nx), yi(ny), u(nx, ny, nu), ut(nx, ny, nu),
     1   ux(nx, ny, nu)
      real uy(nx, ny, nu), uxt(nx, ny, nu), uyt(nx, ny, nu),
     1   a(nx, ny, nu, 2), au(nx, ny, nu, nu, 2), aut(nx, ny, nu, nu, 2)
      real aux(nx, ny, nu, nu, 2), auy(nx, ny, nu, nu, 2), 
     1   auxt(nx, ny, nu, nu, 2), auyt(nx, ny, nu, nu, 2), f(nx, ny, nu)
     2   , fu(nx, ny, nu, nu)
      real fut(nx, ny, nu, nu), fux(nx, ny, nu, nu), fuy(nx,
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
      real t, x(nx), y(ny), lx, rx, ly
      real ry, u(nx, ny, nu), ut(nx, ny, nu), ux(nx, ny, nu)
     1   , uy(nx, ny, nu), uxt(nx, ny, nu)
      real uyt(nx, ny, nu), b(nx, ny, nu), bu(nx, ny, nu, 
     1   nu), but(nx, ny, nu, nu), bux(nx, ny, nu, nu), buy(nx, ny, nu
     2   , nu)
      real buxt(nx, ny, nu, nu), buyt(nx, ny, nu, nu)
      integer i, j
      real r,  COS, ALOG,  SIN,  ATAN, theta
      real  SQRT
      do  6 j = 1, ny
         do  5 i = 1, nx
            if (y(j) .ne. ly) goto 1
               b(i, j, 1) = uy(i, j, 1)
c neumann data on bottom.
               buy(i, j, 1, 1) = 1
               goto  4
   1           r =  SQRT(x(i)**2+y(j)**2)
c dirichlet data.
               if (x(i) .le. 0e0) goto 2
                  theta =  ATAN(y(j)/x(i))
                  goto  3
   2              theta = 2e0* ATAN(1e0)
   3           b(i, j, 1) = u(i, j, 1)-r*( COS(theta)*ALOG(r)-theta*
     1             SIN(theta))
               bu(i, j, 1, 1) = 1
   4        continue
   5        continue
   6     continue
      return
      end
      subroutine ewe2(t, x, nx, y, ny, u, inu, ir)
      integer inu,ir, nx, ny
      real t, x(nx), y(ny), u(nx, ny)
      integer i, j
      real r,  COS, ALOG,  SIN,  ATAN, theta
      real  SQRT
c the exact solution.
         do  6 i = 1, nx
            do  5 j = 1, ny
               r =  SQRT(x(i)**2+y(j)**2)
               if (x(i) .le. 0e0) goto 1
                  theta =  ATAN(y(j)/x(i))
                  goto  2
   1              theta = 2e0* ATAN(1e0)
   2           if (r .le. 0e0) goto 3
                  u(i, j) = r*( COS(theta)*ALOG(r)-theta* SIN(theta))
                  goto  4
   3              u(i, j) = 0
   4           continue
   5           continue
   6        continue
      return
      end
      subroutine handlu(t0, u0, t, u, nv, dt, tstop)
      integer nv
      real t0, u0(nv), t, u(nv), dt, tstop
      common /a7tgup/ errpar, nu, mxp, myp
      integer nu
      real errpar(2)
      common /a7tgum/  kxp,ix,nxp,kyp,iy,nyp,nxnyt,nr,iup
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
      real u(1), t
      common /cstak/ ds
      real ds(500)
      integer ifa, ita(2), ixa(2), nta(2), nxa(2),  ILUMD
      integer ixs, iys, nxs, nys, istkgt, i
      integer iewe, ka(2), ma(2), is(1000), i1mach
      real rs(1000)
      logical ls(1000)
      complex cs(500)
      real  ABS, erru, AMAX1, ws(500)
      integer temp, temp1, temp2
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to get and print the error at each time-step.
c for variable inu for rectangle ir
c u(nx-kx,ny-ky).
c the port library stack and its aliases.
      call enter(1)
c find the error in the solution at 2*kx * 2*ky points / mesh rectangle.
c x search grid.
      ixs =  ILUMD(ws(ix), nx, 2*kx, nxs)
c y search grid.
      iys =  ILUMD(ws(iy), ny, 2*ky, nys)
c u search grid values.
      iewe = istkgt(nxs*nys, 3)
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
      ifa = istkgt(nxs*nys, 3)
c evaluate them.
      call tsd1(2, ka, ws, ita, nta, u, ws, ixa, nxa, ma, ws(ifa))
c error in solution values.
      erru = 0
      temp = nxs*nys
      do  1 i = 1, temp
         temp2 = iewe+i
         temp1 = ifa+i
         erru = AMAX1(erru,  ABS(ws(temp2-1)-ws(temp1-1)))
   1     continue
      temp = i1mach(2)
      write (temp,  2) ir, inu, erru
   2  format(9h for rect,i3,14h error in u(.,,  i2,
     1       3h) =, 1pe10.2)
      call leave
      return
      end
