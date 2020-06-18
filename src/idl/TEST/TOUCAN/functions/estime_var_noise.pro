;+
; NAME:
;        estime_var_noise
;
; PURPOSE:
;     Estimate the variance of the coefficients for each filter of a
;     given wavelet transform from a map of the noise standard
;     deviation in the image domain.
;
; CALLING:
;     var_noise = estime_var_noise(sigmaNoise, bkj, Nscale=Nscale,
;                 npix=npix, lmax=lmax, MonteCarlo=MonteCarlo,
;                 MCiter=MCiter)
;
; INPUTS:
;     SigmaNoise -- Healpix map = map of the noise standard deviation
;                   in the image domain.
;     bkj        -- wavelet transform to be used.
;    
; OUTPUTS:
;     var_noise  -- IDL dblarr[npix, Nscale-1] = Map for each wavelet
;                   filter of the noise coefficient variances.
;
; INPUT KEYWORD: 
;     Nscale      : int = number of scales/filters used in the wavelet
;                   decomposition.
;     npix        : long = Number of pixels of the input map
;                   (12*nside*nside).
;     lmax        : int = maximum l value in the Spherical Harmonic
;                   Space (Healpix).
;     MonteCarlo  : bool = if set, the noise coefficient variances are
;                   computed from Monte-Carlo simulations.
;     MCiter      : int = number of Monte-Carlo iterations to be
;                   used. Default = 20
;
; EXAMPLE:
;      Variances of the wavelet coefficients of a noise map:
;      var_noise = estime_var_noise(sigmaNoise, bkj)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------


;========================== estime_var_noise =============================

function estime_var_noise, sigmaNoise, bkj, Nscale=Nscale, npix=npix, lmax=lmax, MonteCarlo=MonteCarlo, MCiter=MCiter

var_noise = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: Res =  estime_var_noise(sigmaNoise, bkj, Nscale=Nscale, npix=npix, lmax=lmax, MonteCarlo=MonteCarlo, MCiter=MCiter)'
   goto, DONE
endif

if not keyword_set(Nscale) then Nscale = (size(bkj))[2]
if not keyword_set(npix) then npix = N_elements(sigmaNoise)
if not keyword_set(lmax) then lmax = (size(bkj))[1] - 1
if (keyword_set(MonteCalor) && (not keyword_set(MCiter))) then MCiter = 20

ns = N_elements(sigmaNoise)
var_noise = dblarr(npix, Nscale)
ell = findgen(lmax+1)*2.+1.

;; ;------ If stationary noise detected
if (ns eq 1) then begin
   print, 'Stationary noise considered'
   for i=0, Nscale-1 do begin
      var_noise[*,i] = sigmaNoise^2*total(bkj[*,i]^2*ell)/npix
   endfor
endif else begin
   ;; ;------ If non-stationary noise detected
   print, 'Non-stationary noise detected'
   
   if keyword_set(MonteCarlo) then begin
      ;; ;------ If Monte-Carlo option set
      print, 'Estimate Noise variance over ', MCiter, ' Monte-Carlo trials'
      for iter=1,MCiter do begin
         n = randomn(seed,npix)
         ;; compute the needlet decomposition of the noise
         mrs_almtrans, n*sigmaNoise, ALM, lmax=lmax
         ALM_HighResolImag = ALM.alm
         snnt = dblarr(npix, NScale)
         for i=0, Nscale-1 do begin
            h = bkj[*,i]
            alm_product2, ALM_HighResolImag, h, alm_h
            ALM.alm = alm_h
            mrs_almrec, ALM, Hscale
            snnt[*,i] = double(Hscale)
         endfor
         ;;
         var_noise = var_noise + snnt.coef^2
         print, iter, ' Monte-Carlo iteration'
      endfor
      var_noise = var_noise/(MCiter-1)
      
   endif else begin
      ;; ;------ Otherwise, compute the variance directly <---- do a study by patch
      n = randomn(seed,npix)
      noise_spec = mrs_powspec(n*sigmaNoise)
      ellcln = ell*noise_spec[0:lmax]
      for i=0, Nscale-1 do begin
         var_noise[*,i] = total(bkj[*,i]^2*ellcln)/(4.*!dpi)
      endfor
   endelse
endelse

DONE:
return, var_noise

end
