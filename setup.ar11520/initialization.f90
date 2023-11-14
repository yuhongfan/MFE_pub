Subroutine initialize_problemsize
       Use problemsize
       Use input_namelists
       implicit none
!
       inmax=nr+5
       jnmax=nth+5
       knmax=nph+5
       in=inmax
       jn=jnmax
       kn=knmax
!
       return
End Subroutine initialize_problemsize

Subroutine initialize_grid
       Use input_namelists
       Use problemsize
       Use grid
       implicit none
       integer :: i, j, k
       integer :: in_in, jn_in, kn_in
!
       is=3
       js=3
       ks=3
       ie=in-3
       je=jn-3
       ke=kn-3
       ione=1
       jone=1
       kone=1
!
       ism1=is-ione
       jsm1=js-jone
       ksm1=ks-kone
       ism2=is-ione-ione
       jsm2=js-jone-jone
       ksm2=ks-kone-kone
       isp1=is+ione
       jsp1=js+jone
       ksp1=ks+kone
       isp2=is+ione+ione
       jsp2=js+jone+jone
       ksp2=ks+kone+kone
!
       iem1=ie-ione
       jem1=je-jone
       kem1=ke-kone
       iem2=ie-ione-ione
       jem2=je-jone-jone
       kem2=ke-kone-kone
       iep1=ie+ione
       jep1=je+jone
       kep1=ke+kone
       iep2=ie+ione+ione
       jep2=je+jone+jone
       kep2=ke+kone+kone
       iep3=ie+ione+ione+ione
       jep3=je+jone+jone+jone
       kep3=ke+kone+kone+kone
!
       allocate(dx1a(in))
       allocate(dx2a(jn))
       allocate(dx3a(kn))

       allocate(dx1b(in))
       allocate(dx2b(jn))
       allocate(dx3b(kn))

       allocate(dx1ai(in))
       allocate(dx2ai(jn))
       allocate(dx3ai(kn))

       allocate(dx1bi(in))
       allocate(dx2bi(jn))
       allocate(dx3bi(kn))

       allocate(x1a(in))
       allocate(x2a(jn))
       allocate(x3a(kn))

       allocate(x1b(in))
       allocate(x2b(jn))
       allocate(x3b(kn))

       allocate(g2a(in))
       allocate(dg2ad1(in))
       allocate(g2ai(in))

       allocate(g2b(in))
       allocate(dg2bd1(in))
       allocate(g2bi(in))

       allocate(g31a(in))
       allocate(dg31ad1(in))
       allocate(g31ai(in))

       allocate(g31b(in))
       allocate(dg31bd1(in))
       allocate(g31bi(in))

       allocate(g32a(jn))
       allocate(dg32ad2(jn))
       allocate(g32ai(jn))

       allocate(g32b(jn))
       allocate(dg32bd2(jn))
       allocate(g32bi(jn))

       allocate(dvl1a(in))
       allocate(dvl1ai(in))
       allocate(dvl1b(in))
       allocate(dvl1bi(in))

       allocate(dvl2a(jn))
       allocate(dvl2ai(jn))
       allocate(dvl2b(jn))
       allocate(dvl2bi(jn))

       allocate(dvl3a(kn))
       allocate(dvl3ai(kn))
       allocate(dvl3b(kn))
       allocate(dvl3bi(kn))

       open(unit=13, file=trim(gridfile_in), &
        form='unformatted',access='stream')
       read(13) in_in,jn_in,kn_in
       if(in_in .ne. in .or. jn_in .ne. jn &
         .or. kn_in .ne. kn) then
         write(6,*) 'in_in jn_in kn_in= ', &
         in_in, jn_in, kn_in, 'but in jn kn = ', in,jn,kn
         stop
       endif
       read(13) (x1a(i),i=ism2,iep3)
       read(13) (x1b(i),i=ism2,iep2)
       read(13) (x2a(j),j=jsm2,jep3)
       read(13) (x2b(j),j=jsm2,jep2)
       read(13) (x3a(k),k=ksm2,kep3)
       read(13) (x3b(k),k=ksm2,kep2)
       
       read(13) (dx1a(i),i=ism2,iep2)
       read(13) (dx1b(i),i=ism2,iep3)
       read(13) (dx2a(j),j=jsm2,jep2)
       read(13) (dx2b(j),j=jsm2,jep3)
       read(13) (dx3a(k),k=ksm2,kep2)
       read(13) (dx3b(k),k=ksm2,kep3)
       
       read(13) (g2a(i),i=ism2,iep3)
       read(13) (g2b(i),i=ism2,iep2)
       read(13) (g31a(i),i=ism2,iep3)
       read(13) (g31b(i),i=ism2,iep2)
       read(13) (g32a(j),j=jsm2,jep3)
       read(13) (g32b(j),j=jsm2,jep2)

       read(13) (dg2bd1(i),i=ism2,iep3)
       read(13) (dg2ad1(i),i=ism2,iep2)
       read(13) (dg31bd1(i),i=ism2,iep3)
       read(13) (dg31ad1(i),i=ism2,iep2)
       read(13) (dg32bd2(j),j=jsm2,jep3)
       read(13) (dg32ad2(j),j=jsm2,jep2)

       read(13) (dvl1a(i),i=ism2,iep2)
       read(13) (dvl1b(i),i=ism2,iep3)
       read(13) (dvl2a(j),j=jsm2,jep2)
       read(13) (dvl2b(j),j=jsm2,jep3)
       read(13) (dvl3a(k),k=ksm2,kep2)
       read(13) (dvl3b(k),k=ksm2,kep3)

       read(13) (dx1ai(i),i=ism2,iep2)
       read(13) (dx1bi(i),i=ism2,iep3)
       read(13) (dx2ai(j),j=jsm2,jep2)
       read(13) (dx2bi(j),j=jsm2,jep3)
       read(13) (dx3ai(k),k=ksm2,kep2)
       read(13) (dx3bi(k),k=ksm2,kep3)

       read(13) (g2ai(i),i=ism2,iep3)
       read(13) (g2bi(i),i=ism2,iep2)
       read(13) (g31ai(i),i=ism2,iep3)
       read(13) (g31bi(i),i=ism2,iep2)
       read(13) (g32ai(j),j=jsm2,jep3)
       read(13) (g32bi(j),j=jsm2,jep2)

       read(13) (dvl1ai(i),i=ism2,iep2)
       read(13) (dvl1bi(i),i=ism2,iep3)
       read(13) (dvl2ai(j),j=jsm2,jep2)
       read(13) (dvl2bi(j),j=jsm2,jep3)
       read(13) (dvl3ai(k),k=ksm2,kep2)
       read(13) (dvl3bi(k),k=ksm2,kep3)
       close(13)
!
       return
End Subroutine initialize_grid
