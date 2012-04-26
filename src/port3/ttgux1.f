C$TEST  TTGU1
c  main program
c to solve the heat equation with solution u == t*x*y,
c   grad . ( u + ux + .1 * uy, u + uy + .1 * ux ) = ut + ux + uy +g(x,t)
      common /cstak/ ds
      real ds(350000)
      integer ixb(4), iyb(4), nxr(4), nyr(4), kxr(4), kyr(4)
      external handlu, bc, af
      integer ndx, ndy, istkgt, is(1000), iu, nu, kx, ky,  IUMB
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      real tstart, dt, tstop, ws(500)
c the port library stack and its aliases.
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c initialize the port library stack length.
      call istkin(350000, 3)
      call enter(1)
      nu = 1
      kx = 2
      ky = 2
      ndx = 3
      ndy = 3
      nr=4
      tstart = 0.e0
      tstop = 1.e0
      dt = 1
      errpar(1) = 1e-2
      errpar(2) = 1e-4
c uniform grid.
c
c make grid for t-shaped region
c
      ixb(1) =  IUMB(-1.0e0, 0.0e0, ndx, kx, nxr(1))
      ixb(2) =  IUMB(0.0e0, 1.0e0, ndx, kx, nxr(2))
      ixb(3) =  IUMB(1.0e0, 2.0e0, ndx, kx, nxr(3))
      ixb(4) =  IUMB(0.0e0, 1.0e0, ndx, kx, nxr(4))
      iyb(1) =  IUMB(0.0e0, 1.0e0, ndy, ky, nyr(1))
      iyb(2) =  IUMB(0.0e0, 1.0e0, ndy, ky, nyr(2))
      iyb(3) =  IUMB(0.0e0, 1.0e0, ndy, ky, nyr(3))
      iyb(4) =  IUMB(-1.0e0, 0.0e0, ndy, ky, nyr(4))
      nnu =nu*((nxr(1)-kx)*(nyr(1)+nyr(3)-2*ky)+
     1  (nxr(2)-kx)*(nyr(2)-ky)+
     4  (nxr(4)-kx)*(nyr(4)-ky))
      nr=4
c space for the solution.
      iu = istkgt(nnu, 3)
      do 1 i=1,nr
         kxr(i)=kx
         kyr(i)=ky
1     continue
c initial conditions for u.
      call SETR(nnu, 0.e0,ws(iu))
      call  ttgu(ws(iu),nu,nr,kxr,ws,nxr,ixb,kyr,ws,nyr,iyb,tstart,
     1   tstop, dt, af, bc, errpar, handlu)
      call leave
      call wrapup
      stop 
      end
      subroutine af(t, x, nx, y, ny, nu, u, ut, ux, uy, uxt, uyt
     1   , a, au, aut, aux, auy, auxt, auyt, f, fu, fut, fux, fuy, fuxt,
     2   fuyt)
      integer nu, nx, ny
      real t, x(nx), y(ny), u(nx, ny, nu), ut(nx, ny, nu), 
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
      real t, x(nx), y(ny), lx, rx, ly
      real ry, u(nx, ny, nu), ut(nx, ny, nu), ux(nx, ny, nu)
     1   , uy(nx, ny, nu), uxt(nx, ny, nu)
      real uyt(nx, ny, nu), b(nx, ny, nu), bu(nx, ny, nu, 
     1   nu), but(nx, ny, nu, nu), bux(nx, ny, nu, nu), buy(nx, ny, nu
     2   , nu)
      real buxt(nx, ny, nu, nu), buyt(nx, ny, nu, nu)
      integer i, j
      do  2 j = 1, ny
         do  1 i = 1, nx
            bu(i, j, 1, 1) = 1
            b(i, j, 1) = u(i, j, 1)-t*x(i)*y(j)
   1        continue
   2     continue
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
      subroutine ewe2(t, x, nx, y, ny, u, inu, ir)
      integer inu, ir, nx, ny
      real t, x(nx), y(ny), u(nx, ny)
      integer i, j
c the exact solution.
         do  2 i = 1, nx
            do  1 j = 1, ny
               u(i, j) = t*x(i)*y(j)
   1           continue
   2        continue
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
