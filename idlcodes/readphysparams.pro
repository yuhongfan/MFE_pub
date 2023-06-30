PRO readphysparams,file,g_const,rgas,kboltz,mproton,pi,gamma,muconst,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,msol,rsol

;Purpose: read the physparams.dat output file from MFE
;
;Usage: readphysparams,file,g_const,rgas,kboltz,mproton,pi,gamma,muconst,unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,msol,rsol
;file (input): file name (full path file name) for the physparams.dat file
;other arguments (return): the various physical constants and normalizing units (in CGS)

g_const=0.D0
rgas=0.D0
kboltz=0.D0
mproton=0.D0
pi=0.D0
gamma=0.D0
unit_rho=0.D0
unit_len=0.D0
unit_b=0.D0
unit_temp=0.D0
unit_v=0.D0
unit_time=0.D0
msol=0.D0
rsol=0.D0
temp_c=0.D0
rho_c=0.D0
muconst=0.D0

print,'read ',file
openr,/get_lun,unit,file

readu,unit,g_const,rgas,kboltz,mproton,pi,gamma,$
unit_rho,unit_len,unit_b,unit_temp,unit_v,unit_time,$
msol,rsol,temp_c,rho_c,muconst

free_lun,unit

end
