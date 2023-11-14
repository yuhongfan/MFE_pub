PRO read_data,data,time,file,in,jn,kn,endian=endian

if n_elements(endian) eq 0 then endian='little'
in=long(85)
jn=long(85)
kn=long(85)
time=0.D0

print,'read ',file
if endian eq 'little' then begin
openr,/get_lun,unit,file
endif else begin
openr,/get_lun,unit,file,/swap_if_little_endian
endelse
readu,unit,time
readu,unit,in,jn,kn

is=3-1 & js=3-1 & ks=3-1
ie=is+(in-6) & je=js+(jn-6) & ke=ks+(kn-6)
ism1=is-1 & jsm1=js-1 & ksm1=ks-1
ism2=is-2 & jsm2=js-2 & ksm2=ks-2
iep1=ie+1 & jep1=je+1 & kep1=ke+1
iep2=ie+2 & jep2=je+2 & kep2=ke+2
iep3=ie+3 & jep3=je+3 & kep3=ke+3

data=dblarr(iep2-ism2+1,jep2-jsm2+1,kep2-ksm2+1)
readu,unit,data
close,unit
free_lun,unit

return
end
