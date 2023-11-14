Pro remapdisk, array, thgrid, phgrid, radius, diskmap, ycrd, zcrd, ny=ny,nz=nz
;array (input): 2D array for the quantity to be remapped
;thgrid (input): polar angle (rad) coordinate of the first dimension of array
;phgrid (input): azimuthal angle (rad) coordinate of the 2nd dimension of array
;radius (input): radius position of surface where the array is located (currently not used)
;diskmap (output): remapped array of (ny,nz) uniform grid
;ycrd (output): remapped coordinates (deg) in lon
;zcrd (output): remapped coordinates (deg) in lat
;ny (input): number of grid point in lon in the remapped grid
;nz (input): number of grid point in lat in the remapped grid
array=reform(array)
sz=size(array)
nth=n_elements(thgrid)
nph=n_elements(phgrid)
print,'array dimension s(0),s(1),s(2)',sz(0),sz(1),sz(2)

if sz(1) ne nth or sz(2) ne nph then begin
print,'array dimensions', sz(1),sz(2)
print,'do not match sizes of'
print,'thgrid,phgrid ', nth,nph
return
endif

;determine range of the diskmap
;y1=radius*max(sin(thgrid))*min(sin(phgrid))
;y2=radius*max(sin(thgrid))*max(sin(phgrid))
;z1=radius*min(cos(thgrid))
;z2=radius*max(cos(thgrid))
y1=min(phgrid*180./!pi)
y2=max(phgrid*180./!pi)
z1=min(90.-thgrid*180./!pi)
z2=max(90.-thgrid*180./!pi)

;determine resolution and grid for the map
res=(z2-z1)/float(nth)/2.
if n_elements(ny) eq 0 then ny=fix((y2-y1)/res)
if n_elements(nz) eq 0 then nz=fix((z2-z1)/res)
dy=(y2-y1)/float(ny)
dz=(z2-z1)/float(nz)
ycrd=fltarr(ny)
zcrd=fltarr(nz)
diskmap=fltarr(ny,nz)

;re-map
for k=0,nz-1 do begin
print,'k=',k
zp=z1+(0.5+float(k))*dz
zcrd(k)=zp
for j=0,ny-1 do begin
yp=y1+(0.5+float(j))*dy
ycrd(j)=yp
thp=(90.-zp)/180.*!pi
php=yp/180.*!pi
if (thp lt thgrid(0)) or (thp ge thgrid(nth-1)) $
or (php lt phgrid(0)) or (php ge phgrid(nph-1)) then begin
vp=0.
endif else begin
interp2d,array, thgrid, phgrid, thp, php, vp
endelse
diskmap(j,k)=vp
endfor
endfor

return
end
