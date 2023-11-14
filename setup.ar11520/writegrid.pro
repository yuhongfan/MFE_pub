PRO writegrid,in,jn,kn,file, $
x1a,x1b,x2a,x2b,x3a,x3b, $
dx1a,dx1b,dx2a,dx2b,dx3a,dx3b, $
g2a,g2b,g31a,g31b,g32a,g32b, $
dg2bd1,dg2ad1,dg31bd1,dg31ad1,dg32bd2,dg32ad2, $
dvl1a,dvl1b,dvl2a,dvl2b,dvl3a,dvl3b, $
dx1ai,dx1bi,dx2ai,dx2bi,dx3ai,dx3bi, $
g2ai,g2bi,g31ai,g31bi,g32ai,g32bi, $
dvl1ai,dvl1bi,dvl2ai,dvl2bi,dvl3ai,dvl3bi

is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

openw,/get_lun,unit,file
writeu,unit,in,jn,kn
writeu,unit,x1a(ism2:iep3)
writeu,unit,x1b(ism2:iep2)
writeu,unit,x2a(jsm2:jep3)
writeu,unit,x2b(jsm2:jep2)
writeu,unit,x3a(ksm2:kep3)
writeu,unit,x3b(ksm2:kep2)

writeu,unit,dx1a(ism2:iep2)
writeu,unit,dx1b(ism2:iep3)
writeu,unit,dx2a(jsm2:jep2)
writeu,unit,dx2b(jsm2:jep3)
writeu,unit,dx3a(ksm2:kep2)
writeu,unit,dx3b(ksm2:kep3)

writeu,unit,g2a(ism2:iep3)
writeu,unit,g2b(ism2:iep2)
writeu,unit,g31a(ism2:iep3)
writeu,unit,g31b(ism2:iep2)
writeu,unit,g32a(jsm2:jep3)
writeu,unit,g32b(jsm2:jep2)

writeu,unit,dg2bd1(ism2:iep3)
writeu,unit,dg2ad1(ism2:iep2)
writeu,unit,dg31bd1(ism2:iep3)
writeu,unit,dg31ad1(ism2:iep2)
writeu,unit,dg32bd2(jsm2:jep3)
writeu,unit,dg32ad2(jsm2:jep2)

writeu,unit,dvl1a(ism2:iep2)
writeu,unit,dvl1b(ism2:iep3)
writeu,unit,dvl2a(jsm2:jep2)
writeu,unit,dvl2b(jsm2:jep3)
writeu,unit,dvl3a(ksm2:kep2)
writeu,unit,dvl3b(ksm2:kep3)

writeu,unit,dx1ai(ism2:iep2)
writeu,unit,dx1bi(ism2:iep3)
writeu,unit,dx2ai(jsm2:jep2)
writeu,unit,dx2bi(jsm2:jep3)
writeu,unit,dx3ai(ksm2:kep2)
writeu,unit,dx3bi(ksm2:kep3)

writeu,unit,g2ai(ism2:iep3)
writeu,unit,g2bi(ism2:iep2)
writeu,unit,g31ai(ism2:iep3)
writeu,unit,g31bi(ism2:iep2)
writeu,unit,g32ai(jsm2:jep3)
writeu,unit,g32bi(jsm2:jep2)

writeu,unit,dvl1ai(ism2:iep2)
writeu,unit,dvl1bi(ism2:iep3)
writeu,unit,dvl2ai(jsm2:jep2)
writeu,unit,dvl2bi(jsm2:jep3)
writeu,unit,dvl3ai(ksm2:kep2)
writeu,unit,dvl3bi(ksm2:kep3)

free_lun,unit

end
