;+
; NAME:
;        mrs_var_noise_mc
;
; PURPOSE:
;       Estimate the variance of a non-stationary noise coefficients
;       for various wavelet filters using Monte-Carlo simulation of
;       the noise.
;
; CALLING:
;     var_noise = mrs_var_noise_mc(sigmaNoise, bkj, nMC=nMC)
;
; INPUTS:
;     sigmaNoise -- Healpix map  = contains the standard deviation of
;                   the noise in the pixel domain (RMS map).
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;    
; OUTPUTS:
;     var_noise  -- IDL array(npix,Nscale) = maps of the variance of
;                   the noise coefficients for each wavelet filter.
;
; INPUT KEYWORD: 
;     nMC         : int = number of Monte-Carlo simulations used to
;                   estimate the noise coefficients variances.
; EXAMPLE:
;      var_Noise = mrs_var_noise_mc(sigmaNoise, bkj)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------

function mrs_var_noise_mc, sigmaNoise, bkj, nMC=nMC

var_noise = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: Res = mrs_var_noise_mc(sigmaNoise, bkj, nMC=nMC)'
   goto, DONE
endif

if not keyword_set(nMC) then nMC=20
npix = N_elements(sigmaNoise)
nside = npix2nside(npix)
Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1

print, '--- Estimate Noise variance over ', nMC, ' Monte-Carlo trials'
var_noise = dblarr(npix,Nscale)
mean_noise = dblarr(npix,Nscale)
for k=0, nMC-1 do begin
   n = randomn(seed,npix)
   mrs_wt_coeff, sigmaNoise*n, snnt, bkj
   mean_noise = mean_noise + snnt.coef
   var_noise = var_noise + snnt.coef^2.
   print, '- Monte-Carlo iteration' + strcompress(k+1) + ' out of' + strcompress(nMC)
endfor
var_noise = (var_noise - (mean_noise^2.)/nMC) / (nMC -1 )

DONE:
return, var_noise

end

