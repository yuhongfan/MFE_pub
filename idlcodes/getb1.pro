PRO getb1,datadir,it,b1,x1_b1,x2_b1,x3_b1,time_b1

;Purpose: Read one b1 data array inside domain (no ghost zones) at time step it
;
;Usage: getb1,datadir,it,b1,x1_b1,x2_b1,x3_b1,time_b1
;datadir (input): full path of the mfe output data directory
;it (input): time step number
;b1 (return): array for b1
;x1_b1 (return): dimension-1 grid position of the b1 array
;x2_b1 (return): dimension-2 grid position of the b1 array
;x3_b1 (return): dimension-3 grid position of the b1 array
;time_b1 (return): time
;
;Use: readphysparams.pro readgrid.pro

filephysparams=datadir+'/physparams.dat'
readphysparams,filephysparams,g_const,rgas,kboltz,mproton,pi,gamma,muconst,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,msol,rsol

filegrid=datadir+'/grid.dat'
readgrid,filegrid,in,jn,kn, $
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

x1_b1=x1a(is:iep1)
x2_b1=x2b(js:je)
x3_b1=x3b(ks:ke)

fileb1=datadir+'/b1_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
readdata,fileb1,data,time,in,jn,kn
b1=data(is:iep1,js:je,ks:ke)
time_b1=time

end
