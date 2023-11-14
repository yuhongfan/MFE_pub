Pro contfill_cart,y,x1v,x2v,pos=pos,xst=xst,yst=yst,offset=offset,high=high,$
low=low,nl=nl,xrange=xrange, yrange=yrange,xticks=xticks,yticks=yticks,$
color=color,thick=thick,charthick=charthick,charsize=charsize,$
xtitle=xtitle,ytitle=ytitle,title=title
if n_elements(pos) le 0 then pos=[0.12,0.12,1.-0.06,1.-0.06]
if n_elements(xst) le 0 then xst=1
if n_elements(yst) le 0 then yst=1
if n_elements(offset) le 0 then offset=0.
if n_elements(high) le 0 then high=max(y)
if n_elements(low) le 0 then low=min(y)
if n_elements(nl) le 0 then nl=20
if n_elements(xticks) le 0 then xticks=['']
if n_elements(yticks) le 0 then yticks=['']
if n_elements(color) le 0 then color=!p.color
if n_elements(thick) le 0 then thick=1
if n_elements(charthick) le 0 then charthick=1
if n_elements(charsize) le 0 then charsize=1
if n_elements(xtitle) le 0 then xtitle=' '
if n_elements(ytitle) le 0 then ytitle=' '
if n_elements(title) le 0 then title=' '
y=reform(y)
ndms=size(y)
if ndms(0) lt 2 then begin
  print, 'array must be two dimensions'
  return
endif
idms=ndms(1:ndms(0))
isub=where(idms gt 1)
if n_elements(isub) ne 2 then begin
  print, 'array must be two dimensions'
  return
endif
if n_elements(x1v) ne idms(0) or n_elements(x2v) ne idms(1) then begin
  print,'size of x1v must = ',idms(0)
  print,'size of x2v must = ',idms(1)
  return
endif
print,'plotted dimensions: ', isub(0),isub(1)
nx=idms(isub(0))
ny=idms(isub(1))
array=fltarr(nx,ny)
array(*,*)=y
x1arr=x1v#replicate(1.,n_elements(x2v))
x2arr=replicate(1,n_elements(x1v))#x2v

if n_elements(xrange) le 0 then xrange=[min(x1arr),max(x1arr)]
if n_elements(yrange) le 0 then yrange=[min(x2arr),max(x2arr)]

delta=(high-low)*offset
dy=(high-low-2*delta)/float(nl-1)
lv=low+delta+dy*findgen(nl)
print,lv(0),lv(nl-1), format='(2e24.14)'
contour,array,x1arr,x2arr,level=lv,c_linestyle=(lv lt 0.)*1, $
xrange=xrange,yr=yrange,$
xs=xst,ys=yst,pos=pos,xtickname=xticks,xtitle=xtitle,ytitle=ytitle,$
title=title,ytickname=yticks, color=color,thick=thick,charthick=charthick,charsize=charsize,/fill

return
end
