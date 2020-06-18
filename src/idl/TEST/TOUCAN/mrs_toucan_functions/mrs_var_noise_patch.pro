;+
; NAME:
;        mrs_var_noise_patch
;
; PURPOSE:
;       Estimate the variance of a non-stationary noise coefficients
;       for various wavelet filters using a block average over local
;       patches.
;
; CALLING:
;     var_noise = mrs_var_noise_patch(NoiseReal, bkj, bin_size=bin_size)
;
; INPUTS:
;     NoiseReal  -- Healpix map  = contains one noise realization in
;                   the pixel domain.
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;    
; OUTPUTS:
;     var_noise  -- IDL array(npix,Nscale) = maps of the variance of
;                   the noise coefficients for each wavelet filter.
;
; INPUT KEYWORDS:
;     bin_size   : int = size of the patch used to compute the block
;                  average of the noise variances.
;
; EXAMPLE:
;      var_Noise = mrs_var_noise_patch(NoiseReal, bkj)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------


function mrs_var_noise_patch, NoiseReal, bkj, bin_size=bin_size

var_noise = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: Res = mrs_var_noise_patch(NoiseReal, bkj, bin_size=bin_size)'
   goto, DONE
endif

if not keyword_set(bin_size) then bin_size = 16
npix = N_elements(NoiseReal)
nside = npix2nside(npix)
nb_bin = double(nside)/bin_size
TotalNbin = float(nb_bin) * nb_bin * 12.
Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1

var_noise = fltarr(npix, Nscale)
print, 'Computing Needlet decomposition of Noise realization over', Nscale, ' wavelet scales...'
mrs_wt_coeff, NoiseReal, noiset, bkj
print, 'Needlet decomposition done'
coeff = noiset.coef

for k = 0, Nscale-1 do begin
   ScaleData =  H2F(coeff[*,k])
   ScaleRMRS = ScaleData
   ScaleRMRS[*] = 0
   
   Ind = 0L
   for f=0, 11 do begin
      for i=0,nb_bin-1 do begin
         for j=0,nb_bin-1 do begin    
            sig_im =  sigma(ScaleData[i*bin_size :(i+1)*bin_size-1, j*bin_size :(j+1)*bin_size-1, f])
            ScaleRMRS[i*bin_size :(i+1)*bin_size-1, j*bin_size :(j+1)*bin_size-1, f] = sig_im
            Ind = Ind + 1
         endfor
      endfor
   endfor 
   var_noise[*,k] = (F2H(ScaleRMRS))^2.
   
endfor

DONE:
return, var_noise

end

