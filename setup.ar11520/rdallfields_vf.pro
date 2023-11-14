;.run getphysparams1.pro

it=0
print,'it=?'
read,it

file=header+'b1_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_b1,b1,time,file,in,jn,kn,endian=endian

file=header+'b2_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_b2,b2,time,file,in,jn,kn,endian=endian

file=header+'b3_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_b3,b3,time,file,in,jn,kn,endian=endian

file=header+'v1_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,v1,time,file,in,jn,kn,endian=endian

file=header+'v2_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,v2,time,file,in,jn,kn,endian=endian

file=header+'v3_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,v3,time,file,in,jn,kn,endian=endian

file=header+'d_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,d,time,file,in,jn,kn,endian=endian

file=header+'e_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,e,time,file,in,jn,kn,endian=endian

tfield=e*(gamma-1.D0)/d*unit_v*unit_v/(rstar*unit_temp)
p=e*(gamma-1.D0)

;if using hpbcond1
iq=0
print,'input q? yes=1'
read,iq
if iq eq 1 then begin
file=header+'q_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,q,time,file,in,jn,kn,endian=endian
endif

;if using hpbcond
iqflx=0
print,'input qflux? yes=1 '
read,iqflx
if iqflx eq 1 then begin
file=header+'q1_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,q1,time,file,in,jn,kn,endian=endian
file=header+'q2_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,q2,time,file,in,jn,kn,endian=endian
file=header+'q3_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
read_data,q3,time,file,in,jn,kn,endian=endian
endif

end
