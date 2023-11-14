;First run pdfi2mhd_rebin_chromfix.pro in idl
;then run this to plot a vector magnetogram given the 2D arrays br1, bt1, bp1

high=max(abs(br1))  ;maximum br in G
low=-high
maxvec=max(abs(bp1)) ;maximum field strength in G for transverse field vetor plot

highj=1.e-5 ;maximum jr in G/cm
lowj=-highj

ny=800  ;remapped horizontal dimension for the magnetogram
nz=600  ;remapped vertical dimension for the magnetogram

ispx=100 ;vector arrow numbers plotted in the horizontal dimension
ispy=60 ;vector arrow numbers plotted in the vertical dimension

winnum=0 ;window number for the plot
winnumj=1 ;window number for the plot
xwinsz=1200. ;xsize for plotting window

radius=1.
remapdisk,br1,x2b(js:je),x3b(ks:ke),radius,$
brmap,loncrd,latcrd,ny=ny,nz=nz
remapdisk,-bt1,x2a(js:jep1),x3b(ks:ke),1.,$
blatmap,loncrd,latcrd,ny=ny,nz=nz
remapdisk,bp1,x2b(js:je),x3a(ks:kep1),1.,$
blonmap,loncrd,latcrd,ny=ny,nz=nz
remapdisk,jr1_rb,x2a_rb(js_rb:jep1_rb),x3a_rb(ks_rb:kep1_rb),1.,$
jrmap,loncrd,latcrd,ny=ny,nz=nz

ywinsz=long((max(latcrd)-min(latcrd))/(max(loncrd)-min(loncrd))*xwinsz)

print,'ploting magnetic vectors'

fileeps=header+'init_mgram.eps'
filejeps=header+'init_jr_mgram.eps'

filepng=header+'init_mgram.png'
filejpng=header+'init_jr_mgram.png'

mgram_cart,brmap,blonmap,blatmap,loncrd,latcrd,$
pos=[0.2,0.2,0.96,0.96],high=high,low=low,maxvec=maxvec,$
ispx=ispx,ispy=ispy,length=3,nl=256,$
charsize=2,xtit='longitude (deg)',ytit='latitude (deg)',$
winnum=winnum,xwinsz=xwinsz,ywinsz=ywinsz,file=fileeps
write_png,filepng,tvrd(/true)

mgram_cart,jrmap,blonmap,blatmap,loncrd,latcrd,$
pos=[0.2,0.2,0.96,0.96],high=highj,low=lowj,maxvec=maxvec,$
ispx=ispx,ispy=ispy,length=3,nl=256,$
charsize=2,xtit='longitude (deg)',ytit='latitude (deg)',$
winnum=winnumj,xwinsz=xwinsz,ywinsz=ywinsz,file=filejeps
write_png,filejpng,tvrd(/true)

set_plot,'x'
end
