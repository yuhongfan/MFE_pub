;Purpose: Read in a snapshot from the PDFI data to initialize and write out a MFE grid (rebinned), grid_in.dat, physparams.dat file, and the lower boundary normal magnetic field lbdata.dat
;Uses: writegrid.pro

;directory for pdfi_run initial state
header='/glade/derecho/scratch/yfan/ar11520/initdata/'

;input idl save file name for the PDFI electric field data
;file='/glade/campaign/hao/cmemodel/yfan/PDFI/abdullah_2023/save_files/efield_pdfi_ss_11520_20120712_145400.sav'
file='/glade/campaign/hao/cmemodel/yfan/PDFI/abdullah_2023/save_files/efield_pdfi_ss_11520_20120712_133000.sav'
restore,file

tiny=1.D-99
pi=3.1415926535898D0

;PDFI units
clight=2.99792458D10
unit_efield=1.D0/2.99792458D2  ;[stat-Volt]=2.99792458D2[Volt]
dt=12.D0*60.D0

;Physical parameters and units for physparams.dat
g_const = 6.672D-8
rgas = 8.3D7
kboltz = 1.381D-16
mproton = 1.673D-24
mpovme = 1836.15267245D0
muconst = 0.5D0
gamma=1.666666667D0
unit_b=0.20000000000000D+02
unit_rho=0.83650000000000D-15
unit_len=0.69600000000000D+11
unit_temp=0.10000000000000D+07
unit_v=unit_b/sqrt(4.D0*pi*unit_rho)
unit_time=unit_len/unit_v
msol=1.99D33
rsol=6.96D10
gacc_s=g_const*msol/rsol/rsol
rho_c=1.D12*mproton
temp_c=2.D4
rstar=kboltz/(muconst*mproton)

filephys=header+'/physparams.dat'
openw,/get_lun,unit,filephys
writeu, unit, g_const
writeu, unit, rgas
writeu, unit, kboltz
writeu, unit, mproton
writeu, unit, pi
writeu, unit, gamma
writeu, unit, unit_rho
writeu, unit, unit_len
writeu, unit, unit_b
writeu, unit, unit_temp
writeu, unit, unit_v
writeu, unit, unit_time
writeu, unit, msol
writeu, unit, rsol
writeu, unit, temp_c
writeu, unit, rho_c
writeu, unit, muconst
free_lun,unit

;generating grid
jn_orig=n_elements(co_lat_pc_a)+4
kn_orig=n_elements(lon_pc_a)+4
jn=n_elements(co_lat_pc_a)+4
kn=n_elements(lon_pc_a)+4
in=long(432-72+5)
del2=(co_lat_pc_a(-1)-co_lat_pc_a(0))/double(jn-5)
del3=(lon_pc_a(-1)-lon_pc_a(0))/double(kn-5)

del1=3.D7/unit_len
ratr=1.005D0
rmin=1.D0
ic=long(0)
igeo=(in-5)
factorr=double(ic)+(1.D0-ratr^(igeo-1.d0))/(1.D0-ratr)+ratr^(igeo-2.D0)
rmax=rmin+del1*factorr

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

for i=0,is+ic+1 do begin
x1a(i)=double(i-is)*del1+rmin
endfor
for i=is+ic+2,ie do begin
x1a(i)=x1a(i-1)+del1*ratr^(i-is-ic-2)
endfor
x1a(ie+1)=x1a(ie)+del1*ratr^(ie-is-ic-2)
x1a(ie+2)=x1a(ie+1)+del1*ratr^(ie-is-ic-2)
x1a(ie+3)=x1a(ie+2)+del1*ratr^(ie-is-ic-2)
for i=0,in-2 do begin
dx1a(i)=x1a(i+1)-x1a(i)
endfor
for i=0,in-2 do begin
x1b(i)=x1a(i)+0.5D0*dx1a(i)
endfor
for i=1,in-2 do begin
dx1b(i)=x1b(i)-x1b(i-1)
endfor
dx1b(0)=dx1b(1)
dx1b(in-1)=dx1b(in-2)

x2a(js:jep1)=co_lat_pc_a
x2b(js:je)=x2a(js:je)+del2/2.D0
dx2a=replicate(del2,jn-1)
dx2b=replicate(del2,jn)
x2a(jsm1)=x2a(js)-dx2a(jsm1)
x2a(jsm2)=x2a(jsm1)-dx2a(jsm2)
x2a(jep2)=x2a(jep1)+dx2a(jep1)
x2a(jep3)=x2a(jep2)+dx2a(jep2)
x2b(jsm1)=x2a(jsm1)+0.5D0*dx2a(jsm1)
x2b(jsm2)=x2a(jsm2)+0.5D0*dx2a(jsm2)
x2b(jep1)=x2a(jep1)+0.5D0*dx2a(jep1)
x2b(jep2)=x2a(jep2)+0.5D0*dx2a(jep2)

x3a(ks:kep1)=lon_pc_a
x3b(ks:ke)=x3a(ks:ke)+del3/2.D0
dx3a=replicate(del3,kn-1)
dx3b=replicate(del3,kn)
x3a(ksm1)=x3a(ks)-dx3a(ksm1)
x3a(ksm2)=x3a(ksm1)-dx3a(ksm2)
x3a(kep2)=x3a(kep1)+dx3a(kep1)
x3a(kep3)=x3a(kep2)+dx3a(kep2)
x3b(ksm1)=x3a(ksm1)+0.5D0*dx3a(ksm1)
x3b(ksm2)=x3a(ksm2)+0.5D0*dx3a(ksm2)
x3b(kep1)=x3a(kep1)+0.5D0*dx3a(kep1)
x3b(kep2)=x3a(kep2)+0.5D0*dx3a(kep2)

;compute other grid variables
g2a=x1a
g31a=x1a
dg2bd1=replicate(1.D0,in)
dg31bd1=replicate(1.D0,in)
g2b=x1b
g31b=x1b
dg2ad1=replicate(1.D0,in-1)
dg31ad1=replicate(1.D0,in-1)

g32a=abs(sin(x2a))
g32b=abs(sin(x2b))
dg32bd2=cos(x2a)
dg32ad2=cos(x2b)

dvl1a=x1b*x1b*dx1a
dvl1b=x1a*x1a*dx1b
dvl2a = sin(x2b)*dx2a
dvl2b = sin(x2a)*dx2b
dvl3a = dx3a
dvl3b = dx3b

dx1bi = 1.D0 / ( dx1b + tiny )
g2ai = 1.D0 / ( g2a  + tiny )
g31ai = 1.D0 / ( g31a + tiny )
dvl1bi = 1.D0 / ( dvl1b + tiny )

dx1ai = 1.D0 / ( dx1a + tiny )
g2bi  = 1.D0 / ( g2b  + tiny )
g31bi = 1.D0 / ( g31b + tiny )
dvl1ai = 1.D0 / ( dvl1a + tiny )

dx2bi = 1.D0 / ( dx2b + tiny )
g32ai = 1.D0 / ( g32a + tiny )
dvl2bi = 1.D0 / ( dvl2b + tiny )

dx2ai = 1.D0 / ( dx2a + tiny )
g32bi = 1.D0 / ( g32b + tiny )
dvl2ai = 1.D0 / ( dvl2a + tiny )

dx3bi = 1.D0 / ( dx3b + tiny )
dvl3bi = 1.D0 / ( dvl3b + tiny )

dx3ai = 1.D0 / ( dx3a + tiny )
dvl3ai= 1.D0 / ( dvl3a + tiny )

h3_jrorig=g31a(is)*(g32b(js:je)#dx3b(ks:kep1))
h2_jrorig=g2a(is)*(dx2b(js:jep1)#replicate(1.D0,ke-ks+1))
area_jrorig=g2a(is)*g31a(is)*((dx2b(js:jep1)*g32a(js:jep1))#dx3b(ks:kep1))

;trim horizontal grid to jn-5=960 kn-5=1800
jtrm=42
ktrm=36
x2a=x2a(jsm2+jtrm:jep3-jtrm)
x2b=x2b(jsm2+jtrm:jep2-jtrm)
x3a=x3a(ksm2+ktrm:kep3-ktrm)
x3b=x3b(ksm2+ktrm:kep2-ktrm)
dx2a=dx2a(jsm2+jtrm:jep2-jtrm)
dx2b=dx2b(jsm2+jtrm:jep3-jtrm)
dx3a=dx3a(ksm2+ktrm:kep2-ktrm)
dx3b=dx3b(ksm2+ktrm:kep3-ktrm)
g32a=g32a(jsm2+jtrm:jep3-jtrm)
g32b=g32b(jsm2+jtrm:jep2-jtrm)
dg32bd2=dg32bd2(jsm2+jtrm:jep3-jtrm)
dg32ad2=dg32ad2(jsm2+jtrm:jep2-jtrm)
dvl2a=dvl2a(jsm2+jtrm:jep2-jtrm)
dvl2b=dvl2b(jsm2+jtrm:jep3-jtrm)
dvl3a=dvl3a(ksm2+ktrm:kep2-ktrm)
dvl3b=dvl3b(ksm2+ktrm:kep3-ktrm)
dx2ai=dx2ai(jsm2+jtrm:jep2-jtrm)
dx2bi=dx2bi(jsm2+jtrm:jep3-jtrm)
dx3ai=dx3ai(ksm2+ktrm:kep2-ktrm)
dx3bi=dx3bi(ksm2+ktrm:kep3-ktrm)
g32ai=g32ai(jsm2+jtrm:jep3-jtrm)
g32bi=g32bi(jsm2+jtrm:jep2-jtrm)
dvl2ai=dvl2ai(jsm2+jtrm:jep2-jtrm)
dvl2bi=dvl2bi(jsm2+jtrm:jep3-jtrm)
dvl3ai=dvl3ai(ksm2+ktrm:kep2-ktrm)
dvl3bi=dvl3bi(ksm2+ktrm:kep3-ktrm)

;compute jr(jsm2:jep3,ksm2:kep3) where the indexes correspond to the ones after trimming
;note before trimming, x2a(js:jep1)=co_lat_pc_a x3a(ks:kep1)=lon_pc_a
jr0=(bp0(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
*h3_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
-bp0(jtrm-2-1:jn_orig-5-jtrm+2-1,ktrm-2:kn_orig-5-ktrm+2) $
*h3_jrorig(jtrm-2-1:jn_orig-5-jtrm+2-1,ktrm-2:kn_orig-5-ktrm+2) $
-bt0(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
*h2_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
+bt0(jtrm-2:jn_orig-5-jtrm+2,ktrm-2-1:kn_orig-5-ktrm+2-1) $
*h2_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2-1:kn_orig-5-ktrm+2-1)) $
/area_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2)/unit_len

jr1=(bp1(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
*h3_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
-bp1(jtrm-2-1:jn_orig-5-jtrm+2-1,ktrm-2:kn_orig-5-ktrm+2) $
*h3_jrorig(jtrm-2-1:jn_orig-5-jtrm+2-1,ktrm-2:kn_orig-5-ktrm+2) $
-bt1(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
*h2_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2) $
+bt1(jtrm-2:jn_orig-5-jtrm+2,ktrm-2-1:kn_orig-5-ktrm+2-1) $
*h2_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2-1:kn_orig-5-ktrm+2-1)) $
/area_jrorig(jtrm-2:jn_orig-5-jtrm+2,ktrm-2:kn_orig-5-ktrm+2)/unit_len

djrdt=(jr1-jr0)/dt

;compute er(jsm2:jep3,ksm2:kep3) where the indexes correspond to the ones after trimming
er=er(jtrm-2:jn_orig-5+2-jtrm,ktrm-2:kn_orig-5+2-ktrm)

br0=br0(js-js+jtrm:je-js-jtrm,ks-ks+ktrm:ke-ks-ktrm)
bt0=bt0(js-js+jtrm:jep1-js-jtrm,ks-ks+ktrm:ke-ks-ktrm)
bp0=bp0(js-js+jtrm:je-js-jtrm,ks-ks+ktrm:kep1-ks-ktrm)
br1=br1(js-js+jtrm:je-js-jtrm,ks-ks+ktrm:ke-ks-ktrm)
bt1=bt1(js-js+jtrm:jep1-js-jtrm,ks-ks+ktrm:ke-ks-ktrm)
bp1=bp1(js-js+jtrm:je-js-jtrm,ks-ks+ktrm:kep1-ks-ktrm)
et=et(js-js+jtrm:je-js-jtrm,ks-ks+ktrm:kep1-ks-ktrm)
ep=ep(js-js+jtrm:jep1-js-jtrm,ks-ks+ktrm:ke-ks-ktrm)

jn=jn-jtrm*2
kn=kn-ktrm*2
is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

h3=g31a(is)*(g32a(js:jep1)#dx3a(ks:ke))
h2=g2a(is)*(dx2a(js:je)#replicate(1.D0,kep1-ks+1))
area=g2a(is)*g31a(is)*((dx2a(js:je)*g32b(js:je))#dx3a(ks:ke))
area_jr=g2a(is)*g31a(is)*((dx2b(jsm2:jep3)*g32a(jsm2:jep3))#dx3b(ksm2:kep3))

;rebin horizontal grid 2x2
nb=long(2)
in_rb=in
jn_rb=(jn-5)/nb+5
kn_rb=(kn-5)/nb+5
del2_rb=del2*double(nb)
del3_rb=del3*double(nb)

is_rb=3-1 & js_rb=3-1 & ks_rb=3-1
ie_rb=is_rb+(in_rb-6) & je_rb=js_rb+(jn_rb-6) & ke_rb=ks_rb+(kn_rb-6)
ism1_rb=is_rb-1 & jsm1_rb=js_rb-1 & ksm1_rb=ks_rb-1
ism2_rb=is_rb-2 & jsm2_rb=js_rb-2 & ksm2_rb=ks_rb-2
iep1_rb=ie_rb+1 & jep1_rb=je_rb+1 & kep1_rb=ke_rb+1
iep2_rb=ie_rb+2 & jep2_rb=je_rb+2 & kep2_rb=ke_rb+2
iep3_rb=ie_rb+3 & jep3_rb=je_rb+3 & kep3_rb=ke_rb+3
x1a_rb=dblarr(iep3_rb-ism2_rb+1)
x1b_rb=dblarr(iep2_rb-ism2_rb+1)
x2a_rb=dblarr(jep3_rb-jsm2_rb+1)
x2b_rb=dblarr(jep2_rb-jsm2_rb+1)
x3a_rb=dblarr(kep3_rb-ksm2_rb+1)
x3b_rb=dblarr(kep2_rb-ksm2_rb+1)
dx1a_rb=dblarr(iep2_rb-ism2_rb+1)
dx1b_rb=dblarr(iep3_rb-ism2_rb+1)
dx2a_rb=dblarr(jep2_rb-jsm2_rb+1)
dx2b_rb=dblarr(jep3_rb-jsm2_rb+1)
dx3a_rb=dblarr(kep2_rb-ksm2_rb+1)
dx3b_rb=dblarr(kep3_rb-ksm2_rb+1)

x1a_rb=x1a
x1b_rb=x1b
dx1a_rb=dx1a
dx1b_rb=dx1b

for j=js_rb,jep1_rb do begin
x2a_rb(j)=x2a(js+(j-js_rb)*nb)
end
dx2a_rb=replicate(del2_rb,jn_rb-1)
dx2b_rb=replicate(del2_rb,jn_rb)
x2a_rb(jsm1_rb)=x2a_rb(js_rb)-del2_rb
x2a_rb(jsm2_rb)=x2a_rb(jsm1_rb)-del2_rb
x2a_rb(jep2_rb)=x2a_rb(jep1_rb)+del2_rb
x2a_rb(jep3_rb)=x2a_rb(jep2_rb)+del2_rb
x2b_rb(jsm2_rb:jep2_rb)=x2a_rb(jsm2_rb:jep2_rb)+del2_rb/2.D0

for k=ks_rb,kep1_rb do begin
x3a_rb(k)=x3a(ks+(k-ks_rb)*nb)
end
dx3a_rb=replicate(del3_rb,kn_rb-1)
dx3b_rb=replicate(del3_rb,kn_rb)
x3a_rb(ksm1_rb)=x3a_rb(ks_rb)-del3_rb
x3a_rb(ksm2_rb)=x3a_rb(ksm1_rb)-del3_rb
x3a_rb(kep2_rb)=x3a_rb(kep1_rb)+del3_rb
x3a_rb(kep3_rb)=x3a_rb(kep2_rb)+del3_rb
x3b_rb(ksm2_rb:kep2_rb)=x3a_rb(ksm2_rb:kep2_rb)+del3_rb/2.D0

;compute other grid variables
g2a_rb=x1a_rb
g31a_rb=x1a_rb
dg2bd1_rb=replicate(1.D0,in_rb)
dg31bd1_rb=replicate(1.D0,in_rb)
g2b_rb=x1b_rb
g31b_rb=x1b_rb
dg2ad1_rb=replicate(1.D0,in_rb-1)
dg31ad1_rb=replicate(1.D0,in_rb-1)

g32a_rb=abs(sin(x2a_rb))
g32b_rb=abs(sin(x2b_rb))
dg32bd2_rb=cos(x2a_rb)
dg32ad2_rb=cos(x2b_rb)

dvl1a_rb=x1b_rb*x1b_rb*dx1a_rb
dvl1b_rb=x1a_rb*x1a_rb*dx1b_rb
dvl2a_rb = sin(x2b_rb)*dx2a_rb
dvl2b_rb = sin(x2a_rb)*dx2b_rb
dvl3a_rb = dx3a_rb
dvl3b_rb = dx3b_rb

dx1bi_rb = 1.D0 / ( dx1b_rb + tiny )
g2ai_rb = 1.D0 / ( g2a_rb  + tiny )
g31ai_rb = 1.D0 / ( g31a_rb + tiny )
dvl1bi_rb = 1.D0 / ( dvl1b_rb + tiny )

dx1ai_rb = 1.D0 / ( dx1a_rb + tiny )
g2bi_rb  = 1.D0 / ( g2b_rb  + tiny )
g31bi_rb = 1.D0 / ( g31b_rb + tiny )
dvl1ai_rb = 1.D0 / ( dvl1a_rb + tiny )

dx2bi_rb = 1.D0 / ( dx2b_rb + tiny )
g32ai_rb = 1.D0 / ( g32a_rb + tiny )
dvl2bi_rb = 1.D0 / ( dvl2b_rb + tiny )

dx2ai_rb = 1.D0 / ( dx2a_rb + tiny )
g32bi_rb = 1.D0 / ( g32b_rb + tiny )
dvl2ai_rb = 1.D0 / ( dvl2a_rb + tiny )

dx3bi_rb = 1.D0 / ( dx3b_rb + tiny )
dvl3bi_rb = 1.D0 / ( dvl3b_rb + tiny )

dx3ai_rb = 1.D0 / ( dx3a_rb + tiny )
dvl3ai_rb= 1.D0 / ( dvl3a_rb + tiny )

h3_rb=g31a_rb(is_rb)*(g32a_rb(js_rb:jep1_rb)#dx3a_rb(ks_rb:ke_rb))
h2_rb=g2a_rb(is_rb)*(dx2a_rb(js_rb:je_rb)#replicate(1.D0,kep1_rb-ks_rb+1))
area_rb=g2a_rb(is_rb)*g31a_rb(is_rb)*((dx2a_rb(js_rb:je_rb)*g32b_rb(js_rb:je_rb))#dx3a_rb(ks_rb:ke_rb))

br0_rb=dblarr(je_rb-js_rb+1,ke_rb-ks_rb+1)
br1_rb=dblarr(je_rb-js_rb+1,ke_rb-ks_rb+1)
et_rb=dblarr(je_rb-js_rb+1,kep1_rb-ks_rb+1)
ep_rb=dblarr(jep1_rb-js_rb+1,ke_rb-ks_rb+1)
er_rb=dblarr(jep1_rb-js_rb+1,kep1_rb-ks_rb+1)
jr0_rb=dblarr(jep1_rb-js_rb+1,kep1_rb-ks_rb+1)
jr1_rb=dblarr(jep1_rb-js_rb+1,kep1_rb-ks_rb+1)

for k=ks_rb,ke_rb do begin
for j=js_rb,je_rb do begin
br0_rb(j-js_rb,k-ks_rb)=total(br0((j-js_rb)*nb:(j-js_rb)*nb+nb-1,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1) $
*area((j-js_rb)*nb:(j-js_rb)*nb+nb-1,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1)) $
/area_rb(j-js_rb,k-ks_rb)
br1_rb(j-js_rb,k-ks_rb)=total(br1((j-js_rb)*nb:(j-js_rb)*nb+nb-1,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1) $
*area((j-js_rb)*nb:(j-js_rb)*nb+nb-1,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1)) $
/area_rb(j-js_rb,k-ks_rb)
endfor
endfor

for k=ks_rb,kep1_rb do begin
for j=js_rb,jep1_rb do begin
jr0_rb(j-js_rb,k-ks_rb)=total(jr0(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2) $
*area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2)) $
/total(area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2))
jr1_rb(j-js_rb,k-ks_rb)=total(jr1(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2) $
*area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2)) $
/total(area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2))
er_rb(j-js_rb,k-ks_rb)=total(er(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2) $
*area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2)) $
/total(area_jr(js+(j-js_rb)*nb-nb/2:js+(j-js_rb)*nb+nb/2,ks+(k-ks_rb)*nb-nb/2:ks+(k-ks_rb)*nb+nb/2))
endfor
endfor
djrdt_rb=(jr1_rb-jr0_rb)/dt

print,total(br0*area),total(br0_rb*area_rb)
print,total(br1*area),total(br1_rb*area_rb)

for k=ks_rb,kep1_rb do begin
for j=js_rb,je_rb do begin
et_rb(j-js_rb,k-ks_rb)=total(et((j-js_rb)*nb:(j-js_rb)*nb+nb-1,ks-ks+(k-ks_rb)*nb) $
*h2((j-js_rb)*nb:(j-js_rb)*nb+nb-1,ks-ks+(k-ks_rb)*nb)) $
/h2_rb(j-js_rb,k-ks_rb)
endfor
endfor

for k=ks_rb,ke_rb do begin
for j=js_rb,jep1_rb do begin
ep_rb(j-js_rb,k-ks_rb)=total(ep(js-js+(j-js_rb)*nb,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1) $
*h3(js-js+(j-js_rb)*nb,(k-ks_rb)*nb:(k-ks_rb)*nb+nb-1)) $
/h3_rb(j-js_rb,k-ks_rb)
endfor
endfor

;write grid and lbdata
filegrid=header+'/grid_in.dat'
writegrid,in_rb,jn_rb,kn_rb,filegrid, $
x1a_rb,x1b_rb,x2a_rb,x2b_rb,x3a_rb,x3b_rb, $
dx1a_rb,dx1b_rb,dx2a_rb,dx2b_rb,dx3a_rb,dx3b_rb, $
g2a_rb,g2b_rb,g31a_rb,g31b_rb,g32a_rb,g32b_rb, $
dg2bd1_rb,dg2ad1_rb,dg31bd1_rb,dg31ad1_rb,dg32bd2_rb,dg32ad2_rb, $
dvl1a_rb,dvl1b_rb,dvl2a_rb,dvl2b_rb,dvl3a_rb,dvl3b_rb, $
dx1ai_rb,dx1bi_rb,dx2ai_rb,dx2bi_rb,dx3ai_rb,dx3bi_rb, $
g2ai_rb,g2bi_rb,g31ai_rb,g31bi_rb,g32ai_rb,g32bi_rb, $
dvl1ai_rb,dvl1bi_rb,dvl2ai_rb,dvl2bi_rb,dvl3ai_rb,dvl3bi_rb

filelb=header+'/lbdata.dat'
openw,/get_lun,unit,filelb
writeu,unit,jn_rb,kn_rb
writeu,unit,br0_rb
free_lun,unit

end
