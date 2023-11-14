;load a color table from a sav file of rr,gg,bb
filectable='redblue_ctable.dat'
;print,'file for ctable?',filectable
;read,filectable
restore,filectable
tvlct,rr,gg,bb
whitebg
end
