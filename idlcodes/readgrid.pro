PRO readgrid,file,in,jn,kn, $
x1a,x1b,x2a,x2b,x3a,x3b, $
dx1a,dx1b,dx2a,dx2b,dx3a,dx3b, $
g2a,g2b,g31a,g31b,g32a,g32b, $
dg2bd1,dg2ad1,dg31bd1,dg31ad1,dg32bd2,dg32ad2, $
dvl1a,dvl1b,dvl2a,dvl2b,dvl3a,dvl3b, $
dx1ai,dx1bi,dx2ai,dx2bi,dx3ai,dx3bi, $
g2ai,g2bi,g31ai,g31bi,g32ai,g32bi, $
dvl1ai,dvl1bi,dvl2ai,dvl2bi,dvl3ai,dvl3bi

;Purpose: read the grid.dat file to get all the grid variables for the MFE code
;
;Usage: readgrid,file,in,jn,kn,x1a,x1b,x2a,x2b,x3a,x3b,dx1a,dx1b,dx2a,dx2b,dx3a,dx3b,g2a,g2b,g31a,g31b,g32a,g32b,dg2bd1,dg2ad1,dg31bd1,dg31ad1,dg32bd2,dg32ad2,dvl1a,dvl1b,dvl2a,dvl2b,dvl3a,dvl3b,dx1ai,dx1bi,dx2ai,dx2bi,dx3ai,dx3bi,g2ai,g2bi,g31ai,g31bi,g32ai,g32bi,dvl1ai,dvl1bi,dvl2ai,dvl2bi,dvl3ai,dvl3bi
;file (input): file name (full path file name) for the grid.dat file
;other arguments (return): the various grid variables for the MFE code

in=long(1)
jn=long(1)
kn=long(1)
openr,/get_lun,unit,file
readu,unit,in,jn,kn

is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

x1a=dblarr(iep3-ism2+1)
x1b=dblarr(iep2-ism2+1)
x2a=dblarr(jep3-jsm2+1)
x2b=dblarr(jep2-jsm2+1)
x3a=dblarr(kep3-ksm2+1)
x3b=dblarr(kep2-ksm2+1)

dx1a=dblarr(iep2-ism2+1)
dx1b=dblarr(iep3-ism2+1)
dx2a=dblarr(jep2-jsm2+1)
dx2b=dblarr(jep3-jsm2+1)
dx3a=dblarr(kep2-ksm2+1)
dx3b=dblarr(kep3-ksm2+1)

g2a=dblarr(iep3-ism2+1)
g2b=dblarr(iep2-ism2+1)
g31a=dblarr(iep3-ism2+1)
g31b=dblarr(iep2-ism2+1)
g32a=dblarr(jep3-jsm2+1)
g32b=dblarr(jep2-jsm2+1)

dg2bd1=dblarr(iep3-ism2+1)
dg2ad1=dblarr(iep2-ism2+1)
dg31bd1=dblarr(iep3-ism2+1)
dg31ad1=dblarr(iep2-ism2+1)
dg32bd2=dblarr(jep3-ism2+1)
dg32ad2=dblarr(jep2-ism2+1)

dvl1a=dblarr(iep2-ism2+1)
dvl1b=dblarr(iep3-ism2+1)
dvl2a=dblarr(jep2-jsm2+1)
dvl2b=dblarr(jep3-jsm2+1)
dvl3a=dblarr(kep2-ksm2+1)
dvl3b=dblarr(kep3-ksm2+1)

dx1ai=dblarr(iep2-ism2+1)
dx1bi=dblarr(iep3-ism2+1)
dx2ai=dblarr(jep2-jsm2+1)
dx2bi=dblarr(jep3-jsm2+1)
dx3ai=dblarr(kep2-ksm2+1)
dx3bi=dblarr(kep3-ksm2+1)

g2ai=dblarr(iep3-ism2+1)
g2bi=dblarr(iep2-ism2+1)
g31ai=dblarr(iep3-ism2+1)
g31bi=dblarr(iep2-ism2+1)
g32ai=dblarr(jep3-jsm2+1)
g32bi=dblarr(jep2-jsm2+1)

dvl1ai=dblarr(iep2-ism2+1)
dvl1bi=dblarr(iep3-ism2+1)
dvl2ai=dblarr(jep2-jsm2+1)
dvl2bi=dblarr(jep3-jsm2+1)
dvl3ai=dblarr(kep2-ksm2+1)
dvl3bi=dblarr(kep3-ksm2+1)

readu,unit,x1a
readu,unit,x1b
readu,unit,x2a
readu,unit,x2b
readu,unit,x3a
readu,unit,x3b

readu,unit,dx1a
readu,unit,dx1b
readu,unit,dx2a
readu,unit,dx2b
readu,unit,dx3a
readu,unit,dx3b

readu,unit,g2a
readu,unit,g2b
readu,unit,g31a
readu,unit,g31b
readu,unit,g32a
readu,unit,g32b

readu,unit,dg2bd1
readu,unit,dg2ad1
readu,unit,dg31bd1
readu,unit,dg31ad1
readu,unit,dg32bd2
readu,unit,dg32ad2

readu,unit,dvl1a
readu,unit,dvl1b
readu,unit,dvl2a
readu,unit,dvl2b
readu,unit,dvl3a
readu,unit,dvl3b

readu,unit,dx1ai
readu,unit,dx1bi
readu,unit,dx2ai
readu,unit,dx2bi
readu,unit,dx3ai
readu,unit,dx3bi

readu,unit,g2ai
readu,unit,g2bi
readu,unit,g31ai
readu,unit,g31bi
readu,unit,g32ai
readu,unit,g32bi

readu,unit,dvl1ai
readu,unit,dvl1bi
readu,unit,dvl2ai
readu,unit,dvl2bi
readu,unit,dvl3ai
readu,unit,dvl3bi

free_lun,unit

end
