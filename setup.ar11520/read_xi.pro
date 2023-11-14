PRO read_xi,file,data,time,jn,kn,endian=endian

if n_elements(endian) eq 0 then endian='little'
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
readu,unit,jn,kn

data=dblarr(jn,kn)
readu,unit,data
free_lun,unit

return
end
