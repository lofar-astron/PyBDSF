c! roundness modified for python

        subroutine pytess_roundness(n,m,ngens,xgens,ygens,
     /             snrgens,eps,code,volrank)
        implicit none
        integer n,m,ngens,i,areavec(ngens)
        integer roundfacold(ngens),niter
        real*8 volrank(n,m),xgens(ngens),ygens(ngens),snrgens(ngens)
        real*8 eps,roundfac(ngens)
        real*8 roundpix(ngens),x(ngens),y(ngens)
        character code*1

cf2py   intent(in) n,m,ngens,xgens,ygens,snrgens,code,eps
cf2py   intent(out) volrank

        do i=1,ngens
         roundfac(i)=2.d0/3.d0
         x(i)=0.d0
         y(i)=0.d0
         roundpix(i)=0.d0
         areavec(i)=1
        end do

        niter=0
333     continue
        call tess_bin_complicated(n,m,ngens,xgens,ygens,
     /        snrgens,volrank,roundpix,x,y,niter,code,eps)
        niter=niter+1
        do i=1,ngens
         roundfacold(i)=roundfac(i)
        end do
        call tile_roundness(volrank,n,m,ngens,xgens,
     /             ygens,roundfac,roundpix,x,y)
        call calc_area_tess(volrank,n,m,ngens,areavec)
        if (niter.eq.1.or.(niter.lt.2)) then
         goto 333
        end if
        i=int(snrgens(1))

        return
        end

c! same as simple but weights are fn of each pixel now rather than bin
c! if code='s' then each pt belongs to one bin. If not then fuzzy tesselation
        subroutine tess_bin_complicated(n,m,ngens,xgens,ygens,
     /             snrgens,volrank,roundpix,x,y,niter,code,eps)
        implicit none
        integer n,m,ngens,i,j,k,minind(n,m),l,niter
        real*8 volrank(n,m),xgens(ngens),ygens(ngens),dist,dist1
        real*8 dumr,snrgens(ngens),eps,roundpix(ngens)
        real*8 x(ngens),y(ngens),dumr1,wts
        character code*1

        do j=1,m
         do i=1,n
          volrank(i,j)=0.d0
          dumr=1.d90
          do k=1,ngens
           if (niter.eq.0) then
            wts=1.d0
           else
            dumr1=sqrt((x(k)-i)*(x(k)-i)+(y(k)-j)*(y(k)-j))
            dumr1=dumr1*roundpix(k)
            wts=1.d0/dumr1
           end if
           dist=sqrt((i-xgens(k))*(i-xgens(k))+
     /                 (j-ygens(k))*(j-ygens(k)))/wts
           if (dist.lt.dumr) then
            dumr=dist
            minind(i,j)=k       ! minind(i,j) is number of nearest generator
           end if
          end do
         end do
        end do
c!
        if (code.eq.'s') then
         do j=1,m
          do i=1,n
           volrank(i,j)=1.d0*minind(i,j)
          end do
         end do
        else
         do j=1,m
          do i=1,n
           do k=1,ngens
            l=minind(i,j)
            if (k.ne.l) then
             if (niter.eq.0) then
              wts=1.d0
             else
              dumr1=sqrt((x(k)-i)*(x(k)-i)+(y(k)-j)*(y(k)-j))
              dumr1=dumr1*roundpix(k)
              wts=1.d0/dumr1
             end if
             dist=sqrt((i-xgens(k))*(i-xgens(k))+
     /            (j-ygens(k))*(j-ygens(k)))/wts
             dist1=sqrt((i-xgens(minind(i,j)))*(i-xgens(minind(i,j)))+
     /            (j-ygens(minind(i,j)))*(j-ygens(minind(i,j))))/wts
             if (dist.le.(1.d0+eps)*dist1)
     /           volrank(i,j)=volrank(i,j)+1.d0*(minind(i,j)+k)
            end if
           end do
          end do
         end do
        end if
        i=int(snrgens(1))

        return
        end

        subroutine calc_area_tess(volrank,n,m,x,areavec)
        implicit none
        integer n,m,x,areavec(x),i,j
        real*8 volrank(n,m)

        do i=1,x
         areavec(i)=0
        end do
        do j=1,m
         do i=1,n
          areavec(int(volrank(i,j)))=areavec(int(volrank(i,j)))+1
         end do
        end do

        return
        end
c!
c!
c! check roundness.
c! modify to make roundpix not include dist so that u dont have to
c! define huge 3d arrays which crash.
        subroutine tile_roundness(volrank,n,m,ngens,xgens,
     /             ygens,roundfac,roundpix,x,y)
        implicit none
        include "constants.h"
        integer n,m,i,j,ngens,ind,npix(ngens),k
        real*8 volrank(n,m),area(ngens),sumrad(ngens),dist
        real*8 xgens(ngens),ygens(ngens),roundfac(ngens)
        real*8 x(ngens),y(ngens),roundpix(ngens)

        do i=1,ngens
         area(i)=0.d0
         sumrad(i)=0.d0
         npix(i)=0
         x(i)=0.d0
         y(i)=0.d0
        end do

        do j=1,m
         do i=1,n
          ind=int(volrank(i,j))
          x(ind)=x(ind)+i
          y(ind)=y(ind)+j
          npix(ind)=npix(ind)+1
         end do
        end do
        do i=1,ngens
         x(i)=x(i)/npix(i)
         y(i)=y(i)/npix(i)
        end do

        do i=1,ngens
         npix(i)=0
        end do
        do j=1,m
         do i=1,n
          ind=int(volrank(i,j))
          dist=sqrt((xgens(ind)-i)*(xgens(ind)-i)+
     /              (ygens(ind)-j)*(ygens(ind)-j))
          dist=sqrt((x(ind)-i)*(x(ind)-i)+
     /              (y(ind)-j)*(y(ind)-j))
          area(ind)=area(ind)+1
          sumrad(ind)=sumrad(ind)+dist
          npix(ind)=npix(ind)+1
         end do
        end do

        do i=1,ngens
         roundfac(i)=(sumrad(i)/npix(i))/(sqrt(area(i)/pi))
        end do

        do k=1,ngens
         roundpix(k)=1.d0/(sumrad(k)/npix(k))
        end do

        return
        end





