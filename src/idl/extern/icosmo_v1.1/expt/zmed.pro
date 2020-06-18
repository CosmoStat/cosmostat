; function zmed,z,pz
;
; May 06 - Modified by AA to be double precision
; Jan 03 - Written by A. Refregier
;
; PURPOSE: find the median redshift zmed of a redshift distribution
; p(z)
; INPUT: z: vector of redshift values (not necessarily equally spaced
;        pz: corresponding values of p(z) with arbitrary normalisation
; OPTIONAL INPUT: plot: plot differential and cummulative distributions
; OUTPUT: zmed: median of the p(z) distribution defined as:
;                 int_0^zmed dz p(z) / int_0^infinity dz p(z) = 0.5

;--------------------------------------------------------------------------
function func, zz
; function called by root finding routine fx_root.pro

common func_com,z,pint
negz=where(zz le 0.0d)
ppint=interpol(pint,z,zz)
;print,'func: z:',zz
;print,'func: pint;',ppint
ppint = ppint > 0.0d
if (negz(0) ne -1) then ppint(negz)=zz(negz)
return,ppint-.5d
end

;--------------------------------------------------------------------------
function zmed,pz,plotit=plotit

common func_com,z,pint

; declarations
z=pz.z
p=pz.p
n_z=n_elements(z)

; compute integral
pint=dblarr(n_z) ;& pint(0)=p(0)
for i=1L,n_z -1 do pint(i)=pint(i-1)+p(i)
pint=pint/pint(n_z-1)
; compute quick z_med and search limits
dz=z(1)-z(0)

dpint_med=min(abs(pint-.5d),med)
z_m_quick=z(med)

dpint_low=min(abs(pint-.2d),low)
z_low=min([z(low),z_m_quick-2.d*dz])

dpint_high=min(abs(pint-.8d),high)
z_high=max([z(high),z_m_quick+2.d*dz])



; compute z_med accurately
z_m=newton(z_m_quick,'func',/double,stepmax=0.1,itmax=100.)
if keyword_set(plotit) then begin
  plot,z,p,xtitle='z',ytitle='p(z) and p(<z)'
  oplot,z,pint,lines=2
  oplot,[z_m,z_m],[0.,100.],lines=1
  print,'search range: z_low, z_high:',z_low,z_high
  print,'quick: z_m, p(<z_m):',z_m_quick,func(z_m_quick)+.5
  print,'final: z_m, p(<z_m):',z_m,func(z_m)+.5
endif

return,z_m
end

