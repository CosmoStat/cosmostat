;+
; NAME: 
;				PLT_SHEAR
;
; PURPOSE: 
;				plot the shear map of a mass map.
;
; CALLING:
; 				plt_shear, m, /kappa
;
; INPUT: 
; 				m --- mass or shears map (file structure) 
;				
; KEYWORD:
;				kappa --- if mass map input
;				gamma --- if shears map input
; OUTPUT: 
;				overplot the shear on the mass map.
;
; HISTORY:
;				Written by A. Refregier & S. Pires Nov 2005
;-
;-------------------------------------------------------------------------------

pro plt_shear, m , kappa = kappa, gamma = gamma

if keyword_set(kappa) then begin
  sz = size(m.kappa)
  n1= sz[1]
  n2= sz[2]
  kappa_to_gamma, m, g
  gamma1=g.gamma1
  gamma2=g.gamma2
  k = m.kappa
endif

if keyword_set(gamma) then begin
  sz = size(m.gamma1)
  n1= sz[1]
  n2= sz[2]
  gamma1=m.gamma1
  gamma2=m.gamma2
  gamma_to_kappa, m,mm
  k = mm.kappa
endif

window, 0, xsize=500, ysize=575, retain=2
loadct, 15
plt_image, k, /frame,/col

;plot the shear map
for j=0, n2-1, 5 do begin
  for i=0, n1-1, 5 do begin
    g1 = gamma1[i,j]
    g2 = gamma2[i,j]
    ampli = sqrt(((g1)^2) + ((g2)^2))
    alph=(atan(g2,g1))/2. 
    oplot, [i,i + ampli*cos(alph)*40], [j, j + ampli*sin(alph)*40]
    oplot, [i,i - ampli*cos(alph)*40], [j, j - ampli*sin(alph)*40]
  endfor
endfor

end
