;read in grid and physparams.dat
header='header'
print,'header=?'
read,header
endian='little'

;read grid.dat
in=long(85)
jn=long(85)
kn=long(85)
file=header+'grid_in.dat'

if endian eq 'big' then begin
openr,/get_lun,unit,file,/swap_if_little_endian
endif else begin
openr,/get_lun,unit,file
endelse
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

;read physparams.dat
g_const=0.D0
rgas=0.D0
kboltz=0.D0
mproton=0.D0
pi=0.D0
gamma=0.D0
unit_rho=0.D0
unit_len=0.D0
unit_b=0.D0
unit_temp=0.D0
unit_v=0.D0
unit_time=0.D0
msol=0.D0
rsol=0.D0
temp_c=0.D0
rho_c=0.D0
muconst=0.D0

file=header+'physparams.dat'
print,'read ',file
if endian eq 'big' then begin
openr,/get_lun,unit,file,/swap_if_little_endian
endif else begin
openr,/get_lun,unit,file
endelse

readu,unit,g_const,rgas,kboltz,mproton,pi,gamma,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,$
msol,rsol,temp_c,rho_c,muconst

free_lun,unit

help,g_const,rgas,kboltz,mproton,pi,gamma,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,$
msol,rsol,temp_c,rho_c,muconst

rstar=kboltz/(muconst*mproton)
gacc_c=g_const*msol/rsol/rsol
as2=rstar*temp_c/unit_v/unit_v
ovhp0=unit_len/(rstar*temp_c/gacc_c)
temp_hpres=dx1a(is)*unit_len*gacc_c/rstar
hp_prom=2.5e4*rstar/gacc_c

unit_p=unit_rho*unit_v*unit_v
unit_eng=unit_rho*unit_v*unit_v*unit_len^3

help,in-5,jn-5,kn-5

end
