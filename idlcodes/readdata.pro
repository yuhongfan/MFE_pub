PRO readdata,file,data,time,in,jn,kn
;read a field variable data file (including ghostzones)
;Usage: readdata,file,data,time,in,jn,kn
;file (input): file name of the field variable data file
;data (return): the field variable array
;time (return): time of the snapshot
;in (return): maximum grid points in the 1st-dimension
;jn (return): maximum grid points in the 2nd-dimension
;kn (return): maximum grid points in the 3rd-dimension

in=long(1)
jn=long(1)
kn=long(1)
time=0.D0

print,'read ',file
openr,/get_lun,unit,file
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
free_lun,unit

return
end
