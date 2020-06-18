



function anomalies_lowquad, Data, nside=nside, theory=theory, norm=norm
; Written by Anais Rassat Feb 2012. 
; Adapted from code "isw_anomalies.pro Written October 2010 by Anais Rassat
; PURPOSE: Estimates quadrupole of map and returns probability that it
; is low compared to theory
;---------------------------------
; INPUT: 
; map: must be unitless (not mK or muK) overdensity
; nside: nside of input map
; Theory: theoretical value of quadrupole in l(l+1)/(2pi) muK^2
;         default is WMAP3 theory (see paper p32)
; norm: set to 1 if map is in muK
;---------------------------------
; OUTPUT: 
; prob: probability that quadrupole of input map is as low
;---------------------------------
map = Data
if keyword_set(norm) then map = map / 1d6/2.725
if not keyword_set(theory) then theory = 1252.d0 ; which is WMAP3 theory p32
if not keyword_set(nside) then nside = 32

df = 5 ; number of degrees of freedom for quadrupole = 2*l+1
delta = theory/2.d0/(2.d0+1)*2.d0*!dpi
;mrs_almtrans, map, map_alm, /tab, /complex, lmax=lmax ; alm's for ILC map

spec = mrs_powspec(map, lmax = 20)
nuk2 = 1d12*(2.725d)^2        
quad = spec[2]*nuk2 ; units = muK^2
fac = 5.d0

;prob =  imsl_chisqcdf(quad*fac/delta,df,0.d0,/double)
prob = chisqr_pdf(quad*fac/delta, df)

return, prob
end
