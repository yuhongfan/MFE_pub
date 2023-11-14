PRO interp2d, array,x1v,x2v,x1,x2,vp
sz=size(array)
n1=n_elements(x1v)
n2=n_elements(x2v)

if sz(1) ne n1 or sz(2) ne n2 then begin
print,'array dimensions', sz(1),sz(2)
print,'do not match sizes of'
print,'x1v,x2v', n1,n2
return
endif

if (x1 lt x1v(0)) or (x1 ge x1v(n1-1)) $
or (x2 lt x2v(0)) or (x2 ge x2v(n2-1)) then begin
print, 'point out of box'
return
endif

;find x position in the grid
dmin=min(abs(x1v-x1))
isub=where(abs(x1v-x1) eq dmin)
isub=isub(0)
if(x1 lt x1v(isub)) then begin
  i1=isub-1
endif else begin
  i1=isub
endelse
dx1=(x1-x1v(i1))/(x1v(i1+1)-x1v(i1))

;find y position in the grid
dmin=min(abs(x2v-x2))
isub=where(abs(x2v-x2) eq dmin)
isub=isub(0)
if(x2 lt x2v(isub)) then begin
  i2=isub-1
endif else begin
  i2=isub
endelse
dx2=(x2-x2v(i2))/(x2v(i2+1)-x2v(i2))

;checking position
;print,'i1,i2',i1,i2
;print,'x1',x1v(i1),x1,x1v(i1+1)
;print,'x2',x2v(i2),x2,x2v(i2+1)

vp=(1.-dx1)*(1.-dx2)*array(i1,i2) $
  +(1.-dx1)*dx2*array(i1,i2+1) $
  +dx1*(1.-dx2)*array(i1+1,i2) $
  +dx1*dx2*array(i1+1,i2+1)

return
end
