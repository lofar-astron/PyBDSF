
        subroutine gaul2gaulbin(fname)
        implicit none
        character fname*500,word*10,f2*500
        integer gaulid,islid,flag,blc1,blc2,trc1,trc2,srcid,nchar
        real*8 tot,dtot,peak,epeak,ra,era,dec,edec,xpix,expix,ypix
        real*8 eypix,bmaj,ebmaj,bmin,ebmin,bpa,ebpa,dbmaj,edbmaj,dbmin
        real*8 edbmin,dbpa,edbpa,sstd,sav,rstd,rav,chisq,q,dumr1,dumr2
        real*8 dumr3,dumr4,dumr5,dumr6

cf2py   fname

        f2 = fname(1:nchar(fname))//'.bin'
        open(unit=21, file=fname)
        open(unit=22, file=f2, form='unformatted')
        word = ''

        do while (word.ne.'fmt')
         read (21,*) word
        end do
200     read (21,*,END=100) gaulid,islid,flag,tot,dtot,peak,epeak,ra,
     /   era,dec,edec,xpix,expix,ypix,eypix,bmaj,ebmaj,bmin,ebmin,bpa,
     /   ebpa,dbmaj,edbmaj,dbmin,edbmin,dbpa,edbpa,sstd,sav,rstd,
     /   rav,chisq,q,srcid,blc1,blc2,trc1,trc2,dumr1,dumr2,
     /   dumr3,dumr4,dumr5,dumr6
        write (22) gaulid,islid,flag,tot,dtot,peak,epeak,ra,
     /   era,dec,edec,xpix,expix,ypix,eypix,bmaj,ebmaj,bmin,ebmin,bpa,
     /   ebpa,dbmaj,edbmaj,dbmin,edbmin,dbpa,edbpa,sstd,sav,rstd,
     /   rav,chisq,q,srcid,blc1,blc2,trc1,trc2,dumr1,dumr2,
     /   dumr3,dumr4,dumr5,dumr6
        goto 200

100     close(21)
        close(22)

        return
        end



