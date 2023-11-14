Pro mgram_cart,vz,vx,vy,xx,yy,pos=pos,$
offset=offset,high=high,low=low,nl=nl,ispx=ispx,ispy=ispy,length=length,maxvec=maxvec,$
charsize=charsize,thick=thick,$
xtit=xtit,ytit=ytit,tit=tit,winnum=winnum,xwinsz=xwinsz,ywinsz=ywinsz,$
vztit=vztit,file=file
tiny=1.e-30
if n_elements(pos) le 0 then pos=[0.12,0.12,1.-0.06,1.-0.06]
if n_elements(offset) le 0 then offset=0.
if n_elements(high) le 0 then high=max(abs(vz))*1.1
if n_elements(low) le 0 then low=-max(abs(vz))*1.1
if n_elements(nl) le 0 then nl=20
if n_elements(ispx) le 0 then ispx=40
if n_elements(ispy) le 0 then ispy=40
if n_elements(length) le 0 then length=1
if n_elements(maxvec) le 0 then maxvec=max(sqrt(vx^2+vy^2))+tiny
if n_elements(charsize) le 0 then charsize=1
if n_elements(thick) le 0 then thick=1
if n_elements(xtit) le 0 then xtit=' '
if n_elements(ytit) le 0 then ytit=' '
if n_elements(tit) le 0 then tit=' '
if n_elements(winnum) le 0 then winnum=0
if n_elements(xwinsz) le 0 then xwinsz=600
if n_elements(ywinsz) le 0 then ywinsz=600
if n_elements(vztit) le 0 then vztit=' '

vx=reform(vx)
vy=reform(vy)
vz=reform(vz)
ndms=size(vz)
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
if n_elements(xx) ne idms(0) or n_elements(yy) ne idms(1) then begin
  print,'size of xx must = ',idms(0)
  print,'size of yy must = ',idms(1)
  return
endif
print,'plotted dimensions: ', isub(0),isub(1)

nx=idms(isub(0))
ny=idms(isub(1))
x1arr=xx#replicate(1.,n_elements(yy))
x2arr=replicate(1,n_elements(xx))#yy
delta=(high-low)*offset
dy=(high-low-2*delta)/float(nl-1)
lv=low+delta+dy*findgen(nl)
print,lv(0),lv(nl-1), format='(2e24.14)'

;clipping vz for values above high and below low
isub=where(vz gt high/1.01)
vz(isub)=high/1.01
isub=where(vz lt low/1.01)
vz(isub)=low/1.01

maxvecloc=max(sqrt((rebin(vx,ispx,ispy))^2+(rebin(vy,ispx,ispy))^2))+tiny

window,winnum,xsize=xwinsz,ysize=ywinsz
filectable='redblue_ctable.dat'
restore,filectable
rr(0)=255
gg(0)=255
bb(0)=255
tvlct,rr,gg,bb
contour,vz,x1arr,x2arr,level=lv,c_linestyle=(lv lt 0.)*1, $
xr=[min(x1arr),max(x1arr)],yr=[min(x2arr),max(x2arr)],$
xs=5,ys=5,pos=pos,/fill
loadct,0
plot,xx,yy,pos=pos,xs=1,ys=1,/noer,/nodata,charsize=charsize, $
xtit=xtit,ytit=ytit,tit=tit,color=1
velovect,rebin(vx,ispx,ispy),rebin(vy,ispx,ispy),rebin(xx,ispx),rebin(yy,ispy),$
color=1,length=length*maxvecloc/maxvec,xs=5,ys=5,pos=pos,/noer

;color table
tvlct,rr,gg,bb

xo=xwinsz*0.05
yo=ywinsz*0.1
xsz=xwinsz*0.15
ysz=15.
table=low+findgen(xsz)/float(xsz-1)*(high-low)#replicate(1.,ysz)
postbl=[xo/xwinsz,yo/ywinsz,(xo+xsz)/xwinsz,(yo+ysz)/ywinsz]
tv,bytscl(table,max=high,min=low),xo,yo,xsize=xsz,ysize=ysz
loadct,0
plot,[low,high],[0,1],xs=1,ys=5,pos=postbl,/noer,/nodata,$
color=1,xticks=1,charsize=charsize
plots,[low,low],[0,1],color=1
plots,[high,high],[0,1],color=1
xyouts,(xo+xsz*1.1)/xwinsz,yo/ywinsz,/normal,vztit,color=1

;write a .eps file
if n_elements(file) gt 0 then begin
set_plot,'ps'
device,xsize=4,ysize=4*ywinsz/xwinsz,xoffset=1,yoffset=1,/in
device,/col
device,bits=8
device,/enc
device,file=file
tvlct,rr,gg,bb
font=0

contour,vz,x1arr,x2arr,level=lv,c_linestyle=(lv lt 0.)*1, $
xr=[min(x1arr),max(x1arr)],yr=[min(x2arr),max(x2arr)],$
xs=5,ys=5,pos=pos,/fill
loadct,0
plot,xx,yy,pos=pos,xs=1,ys=1,/noer,/nodata,charsize=1, $
xtit=xtit,ytit=ytit,tit=tit,color=1,font=font
velovect,rebin(vx,ispx,ispy),rebin(vy,ispx,ispy),rebin(xx,ispx),rebin(yy,ispy),$
color=1,length=length*maxvecloc/maxvec,xs=5,ys=5,pos=pos,/noer

;color table
tvlct,rr,gg,bb
xo=0.05
yo=0.1
xsz=0.15
ysz=0.025
table=low+findgen(100)/float(100-1)*(high-low)#replicate(1.,15)
postbl=[xo,yo,xo+xsz,yo+ysz]
tv,bytscl(table,max=high,min=low),xo,yo,xsize=xsz,ysize=ysz,/normal
loadct,0
plot,[low,high],[0,1],xs=1,ys=5,pos=postbl,/noer,/nodata,$
color=1,xticks=1,charsize=1,font=font
plots,[low,low],[0,1],color=1
plots,[high,high],[0,1],color=1
;xyouts,0.07,0.15,/normal,vztit,color=1,charsize=1
xyouts,(xo+xsz*1.1),yo,/normal,vztit,charsize=1,charthick=2,color=1

device,/close
set_plot,'x'
endif

return
end
