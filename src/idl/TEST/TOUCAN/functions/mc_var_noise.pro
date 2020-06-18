;+
; NAME:
;        mc_var_noise
;
; PURPOSE:
;       Estimate the variance of a non-stationary noise coefficients
;       for various wavelet filters using Monte-Carlo simulation of
;       the noise.
;
; CALLING:
;     var_noise = mc_var_noise(sigmaNoise, Nscale=Nscale,
;     Firstl=Firstl, lmax=lmax, nMC=nMC, lin=lin, filters=bkj)
;
; INPUTS:
;     sigmaNoise -- Healpix map  = contains the standard deviation of
;                   the noise in the pixel domain.
;    
; OUTPUTS:
;     var_noise  -- IDL array(npix,Nscale) = maps of the variance of
;                   the noise coefficients for each wavelet filter.
;
; INPUT KEYWORD: 
;     Nscale      : int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;     Firstl      : int = initial multipole value to be reconstructed
;                   in Cl.
;     lmax        : int = final multipole value to be reconstructed in
;                   Cl.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     nMC         : int = number of Monte-Carlo simulations used to
;                   estimate the noise coefficients variances.
;     lin         : bool = if set, the width and spacing of the filters in the
;                   multipole domain is linear.
;     filters     : array of size lmax+1 by Nscale = wavelet filters in the
;                   multipole domain to be used to compute the compressed
;                   measurements.
;
; EXAMPLE:
;      var_Noise = mc_var_noise(sigmaNoise, Nscale=Nscale,
;      Firstl=Firstl, lmax=lmax)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------

function mc_var_noise, sigmaNoise, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtwidth=filtwidth, nMC=nMC, lin=lin, filters=bkj

if not keyword_set(nMC) then nMC=20
npix = N_elements(sigmaNoise)
nside = npix2nside(npix)

if keyword_set(bkj) then begin
   Nscale = (size(bkj))[2]
   lmax = (size(bkj))[1]-1
endif
print, '--- Estimate Noise variance over ', nMC, ' Monte-Carlo trials'
var_noise = dblarr(npix,Nscale)
mean_noise = dblarr(npix,Nscale)
for k=0, nMC-1 do begin
   n = randomn(seed,npix)
   ntousi_needlet, sigmaNoise*n, snnt, NScale=Nscale, Firstl=2, lmax=lmax, filtwidth=filtwidth, lin=lin, filters=bkj
   mean_noise = mean_noise + snnt.coef
   var_noise = var_noise + snnt.coef^2.
   print, '- Monte-Carlo iteration' + strcompress(k+1) + ' out of' + strcompress(nMC)
endfor
var_noise = (var_noise - (mean_noise^2.)/nMC) / (nMC -1 )

return, var_noise

end

