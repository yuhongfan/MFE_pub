;specify the <XIDATA_DIR> path
header='/glade/derecho/scratch/yfan/ar11520/xifield/'

nt1=long(1)
nt2=long(1000)
nstep=long(1)

;read in the first xi field file to get the horizontal grid size
fileslc=header+'/xi_'+strtrim(string(nt1,format='(i4.4)'),1)+'.dat'
print,'read ',fileslc
read_xi,fileslc,xi,time,jn_in,kn_in
jn=long(jn_in)
kn=long(kn_in)
ntsub=long(50)
xislc1=dblarr(ntsub,jn,kn)
print,'jn,kn=',jn,kn

;open file for time table
filetimetbl=header+'timetable'
openw,/get_lun,unit_tbl,filetimetbl

tsub=0.D0
itsub=long(0)
itcount=long(0)
ifilesub=long(1)

for it=nt1,nt2,nstep do begin
fileslc=header+'/xi_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
print,'read ',fileslc
read_xi,fileslc,xi,time,jn_in,kn_in

tsub=time
itsub=itsub+1
itcount=itcount+1

xislc1(itcount-1,*,*)=xi

printf,unit_tbl,tsub,itsub,itcount,ifilesub,format='(e22.14,3i8)'
print,tsub,itsub,itcount,ifilesub,format='(e22.14,3i8)'

if(itcount eq ntsub) then begin

filexislc1=header+'xislc1_23.'+strtrim(string(ifilesub,format='(i4.4)'),1)+'.dat'
write_23qty,xislc1,ntsub,jn,kn,filexislc1

ifilesub=ifilesub+1
itcount=1
xislc1(itcount-1,*,*)=xi

endif

endfor

free_lun,unit_tbl

filexislc1=header+'xislc1_23.'+strtrim(string(ifilesub,format='(i4.4)'),1)+'.dat'
write_23qty,xislc1,ntsub,jn,kn,filexislc1

end
