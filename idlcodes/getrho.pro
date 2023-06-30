PRO getrho,datadir,it,rho,x1_rho,x2_rho,x3_rho,time_rho

;Purpose: Read one density data array inside domain (no ghost zones) at time step it
;
;Usage: getrho,datadir,it,rho,x1_rho,x2_rho,x3_rho,time_rho
;datadir (input): full path of the mfe output data directory
;it (input): time step number
;rho (return): array for density
;x1_rho (return): dimension-1 grid position of the rho
;x2_rho (return): dimension-2 grid position of the rho
;x3_rho (return): dimension-3 grid position of the rho
;time_rho (return): time
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

x1_rho=x1b(is:ie)
x2_rho=x2b(js:je)
x3_rho=x3b(ks:ke)

filed=datadir+'/d_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
readdata,filed,data,time,in,jn,kn
rho=data(is:ie,js:je,ks:ke)
time_rho=time

end
