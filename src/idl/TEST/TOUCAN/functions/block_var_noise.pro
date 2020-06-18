;+
; NAME:
;        block_var_noise
;
; PURPOSE:
;       Estimate the variance of a non-stationary noise coefficients
;       for various wavelet filters using a block average over local
;       patches.
;
; CALLING:
;     var_noise = block_var_noise(sigmaNoise, Nscale=Nscale,
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
;     bin_size    : int = size of the patches used to estimate the
;                   noise coefficients variances.
;     lin         : bool = if set, the width and spacing of the filters in the
;                   multipole domain is linear.
;     filters     : array of size lmax+1 by Nscale = wavelet filters in the
;                   multipole domain to be used to compute the compressed
;                   measurements.
;
; EXAMPLE:
;      var_Noise = bock_var_noise(sigmaNoise, Nscale=Nscale,
;      Firstl=Firstl, lmax=lmax)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------


function block_var_noise, NoiseReal, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtwidth=filtwidth, bin_size=bin_size, lin=lin, filters=bkj

if not keyword_set(bin_size) then bin_size = 16
npix = N_elements(NoiseReal)
nside = npix2nside(npix)
nb_bin = double(nside)/bin_size
TotalNbin = float(nb_bin) * nb_bin * 12.

if keyword_set(bkj) then begin
   Nscale = (size(bkj))[2]
   lmax = (size(bkj))[1]-1
endif
var_noise = fltarr(npix, Nscale)
print, 'Computing Needlet decomposition of Noise realization over', Nscale, ' wavelet scales...'
ntousi_needlet, NoiseReal, noiset, NScale=Nscale, Firstl=Firstl, lmax=lmax, filtwidth=filtwidth, lin=lin, filters=bkj
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

;; filename = 'block_varnoise_'+strcompress(nside,/remove_all)+'_'+strcompress(Nscale,/remove_all)+'_'+strcompress(nMC,/remove_all)+'.save'
;; filepresent = file_search(filename,count=c)

return, var_noise

end

