C$TEST  DTTGU3
c  main program
      common /cstak/ ds
      double precision ds(350000)
      external handlu, bc, af
      integer ndx, ndy, istkgt, is(1000), iu
      integer nu, nr, iyb(3), ixb(3), kx, ky
      integer nxr(3), nyr(3), kxr(3), kyr(3)
      integer idumb
      real errpar(2), rs(1000)
      logical ls(1000)
      complex cs(500)
      double precision tstart, dt, ws(500)
      double precision tstop
      equivalence (ds(1), cs(1), ws(1), rs(1), is(1), ls(1))
c to solve the layered heat equation, with kappa = 1, 1/2, 1/3,
c   div . ( kappa(x,y) * grad u ) = ut + g
c the port library stack and its aliases.
c initialize the port library stack length.
      call istkin(350000, 4)
      call enter(1)
      nu = 1
      nr = 3
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
      ixb(1) = idumb(0.0d0, 1.0d0, ndx, kx, nxr(1))
      ixb(2) = idumb(0.0d0, 1.0d0, ndx, kx, nxr(2))
      ixb(3) = idumb(0.0d0, 1.0d0, ndx, kx, nxr(3))
      iyb(1) = idumb(0.0d0, 1.0d0, ndy, ky, nyr(1))
      iyb(2) = idumb(1.0d0, 2.0d0, ndy, ky, nyr(2))
      iyb(3) = idumb(2.0d0, 3.0d0, ndy, ky, nyr(3))
c space for the solution.
      nnu=0
      do 1 i=1,nr
         nnu=nnu+nu*((nxr(i)-kx)*(nyr(i)-ky))
 1    continue
      iu = istkgt(nnu, 4)
      do 2 i=1,nr
         kxr(i)=kx
         kyr(i)=ky
 2    continue
      call setd(nnu, 0.d0,ws(iu))
      call dttgu(ws(iu),nu,nr,kxr,ws,nxr,ixb,kyr,ws,nyr,iyb,tstart,
     1   tstop, dt, af, bc, errpar, handlu)
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
      double precision kappa
      logical temp
      do  7 i = 1, nu
         do  6 q = 1, ny
            do  5 p = 1, nx
               if (y(q) .ge. 1d0) goto 1
                  kappa = 1
                  goto  4
   1              if (y(q) .ge. 2d0) goto 2
                     kappa = 0.5
                     goto  3
   2                 kappa = 1d0/3d0
   3           continue
   4           a(p, q, i, 1) = kappa*ux(p, q, i)
               aux(p, q, i, i, 1) = kappa
               a(p, q, i, 2) = kappa*uy(p, q, i)
               auy(p, q, i, i, 2) = kappa
               f(p, q, i) = ut(p, q, i)
               fut(p, q, i, i) = 1
               f(p, q, i) = f(p, q, i)-y(q)/kappa
               temp = 1d0 .lt. y(q)
               if (temp) temp = y(q) .lt. 2d0
               if (temp) f(p, q, i) = f(p, q, i)+1d0
               temp = 2d0 .lt. y(q)
               if (temp) temp = y(q) .lt. 3d0
               if (temp) f(p, q, i) = f(p, q, i)+3d0
   5           continue
   6        continue
   7     continue
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
      logical temp
      do  6 j = 1, ny
         do  5 i = 1, nx
            temp = x(i) .eq. lx
            if (.not. temp) temp = x(i) .eq. rx
            if (.not. temp) goto 1
               bux(i, j, 1, 1) = 1
c left or right.
c neumann bcs.
               b(i, j, 1) = ux(i, j, 1)
               goto  4
   1           if (y(j) .ne. ly) goto 2
                  b(i, j, 1) = u(i, j, 1)
c bottom.
                  bu(i, j, 1, 1) = 1
                  goto  3
   2              b(i, j, 1) = u(i, j, 1)-6d0*t
c top.
                  bu(i, j, 1, 1) = 1
   3        continue
   4        continue
   5        continue
   6     continue
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
      subroutine ewe2(t, x, nx, y, ny, u, inu, ir)
      integer inu, ir, nx, ny
      double precision t, x(nx), y(ny), u(nx, ny), dble
      integer i, j
c the exact solution.
         do  6 i = 1, nx
            do  5 j = 1, ny
              u(i, j) = dble(float(ir))*t*y(j)-dble(float(ir-1))*t
              if(ir.eq.3) u(i,j)=u(i,j)-t
   5           continue
   6        continue
      return
      end
