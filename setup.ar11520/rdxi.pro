header='/export/data1/yfan/stitch_hinode_vg/'
print,'header=?'
read,header

it=long(1)
print,'it=?'
read,it

filexi=header+'/xi_'+strtrim(string(it,format='(i4.4)'),1)+'.dat'
print,'read ',filexi
read_xi,filexi,xi,time,jn,kn

end
