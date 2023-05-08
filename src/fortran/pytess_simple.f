c! if code='s' then each pt belongs to one bin. If not then fuzzy tesselation
c! cut out just this part for pythn cos it takes forever to do this part
c! in python

        subroutine pytess_simple(n,m,ngens,xgens,ygens,
     /             snrgens,wts,eps,code,volrank)
        implicit none
        integer n,m,ngens,i,j,k,minind(n,m),l
        real*8 volrank(n,m),xgens(ngens),ygens(ngens),dist
        real*8 dumr,snrgens(ngens),wts(ngens),eps,distmin
        character code*1

cf2py   intent(in) n,m,ngens,xgens,ygens,snrgens,wts,code,eps
cf2py   intent(out) volrank

        do j=1,m
         do i=1,n
          volrank(i,j)=0.d0
          dumr=1.d90
          do k=1,ngens
           dist=sqrt((i-xgens(k))*(i-xgens(k))+
     /               (j-ygens(k))*(j-ygens(k)))/wts(k)
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
             dist=sqrt((i-xgens(k))*(i-xgens(k))+
     /                 (j-ygens(k))*(j-ygens(k)))/wts(k)
             distmin=sqrt((i-xgens(l))*(i-xgens(l))+
     /                    (j-ygens(l))*(j-ygens(l)))/wts(l)
             if (dist.le.(1.d0+eps)*distmin)
     /           volrank(i,j)=volrank(i,j)+1.d0*(l+k)
            end if
           end do
          end do
         end do
        end if

        i = int(snrgens(1))  ! so there is no warning while compiling

        return
        end
c!
