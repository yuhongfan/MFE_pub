PRO write_23qty,data,ntsub,jn,kn,file,endian=endian

if n_elements(endian) eq 0 then endian='little'

print,'write ',file
print,'ntsub,jn,kn=',ntsub,jn,kn

if endian eq 'little' then begin
openw,/get_lun,unit,file
endif else begin
openw,/get_lun,unit,file,/swap_if_little_endian
endelse

writeu,unit,ntsub,jn,kn
writeu,unit,data
free_lun,unit

return
end
