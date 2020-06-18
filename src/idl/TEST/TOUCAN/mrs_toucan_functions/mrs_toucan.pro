;+
; NAME:
;        mrs_toucan
;
; PURPOSE:
;  Computes the CMB power spectrum from a CMB map. Takes into account
;  non-stationary noise given by the RMS map SigmaNoise of the noise
;  variance in the pixel domain, or by a noise realization NoiseReal
;  in the pixel domain or direcly by the variance of the noise
;  coefficient for each wavelet filter in var_noise. For low values of
;  the multipole l (values below tousil), mrs_toucan uses the result
;  given by mrs_tousi and applies a window to merge the two solutions.
;
; CALLING:
;     toucanCl = mrs_toucan( Map, filters=bkj,
;     var_noise=var_noise, NoiseReal=NoiseReal, sigmaNoise=sigmaNoise,
;     Mask=Mask, Niter=Niter, Nscale=Nscale, Firstl=Firstl, lmax=lmax,
;     Filtwidth=l, lin=lin, tousil=tousil, inpainting=inpainting,
;     scaledMask=scaledMask, rescaleMask=rescaleMask, WTNscale=WTNscale, SWT=SWT, HWT=HWT,
;     OptWT=OptWT, Pos=Pos)
;
; INPUTS:
;     Map        -- Healpix map: Input CMB image to estimate the CMB
;                   power spectrum from.
;    
; OUTPUTS:
;     toucanCl   -- IDL 1D array  = Estimated Cl array.
;
; INPUT KEYWORD: 
;     filters     : array of size lmax+1 by Nscale = wavelet filters
;                   in the multipole domain to be used to compute the
;                   compressed measurements.
;     var_noise   : array containing Nscale Healpix maps = each map
;                   contains the variance of the noise at that wavelet
;                   scale.
;     NoiseReal   : Healpix map  = contains one noise realization in
;                   the pixel domain.
;     sigmaNoise  : Healpix map  = contains the standard deviation of
;                   the noise in the pixel domain (RMS map).
;     Mask        : Healpix map = mask applied to the input Map.
;     Niter       : int = number of toucan iterations for the reconstruction.
;     Nscale      : int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;     Firstl      : int = initial multipole value to be reconstructed
;                   in Cl.
;     lmax        : int = final multipole value to be reconstructed in
;                   Cl.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     lin         : bool = if set, the width and spacing of the filters in the
;                   multipole domain is linear.
;     tousil      : int = value of the multipole l below which the
;                   solution is an average between Tousi and Toucan.
;     inpainting  : bool = if set, an impainting of the map is
;                   performed before estimation of the compressed
;                   measurements.
;     scaledMask  : array containing Nscale Healpix maps = mask to be
;                   used in the computation of the compressed
;                   measurements at each wavelet scale.
;     RescaleMask : bool = if set, a rescaling of the mask is
;                   performed at each wavelet scale during the
;                   estimation of the compressed measurements.
;     WTNscale    : int = number of wavelet scales used to impose
;                   sparsity during reconstruction.
;     SWT         : bool = if set, soft-thresholding of the wavelet
;                   coefficients of the power spectrum is performed
;                   during reconstruction.
;     HWT         : bool = if set, hard-thresholding of the wavelet
;                   coefficients of the power spectrum is performed
;                   during reconstruction.
;     OptWT       : options to be used in the wavelet decomposition of
;                   the power spectrum during reconstruction.
;     Pos         : bool = if set, a positivity constraint is applied
;                   to the power spectrum during reconstruction.
;
; OUTPUT KEYWORD: 
;
; EXAMPLE:
;      Estimate the CMB power spectrum from a CMB map with stationary
;      noise with standard deviation of 5:
;      toucan_Cl =  toucan(Map)
;
;      Estimate the CMB power spectrum from a map with non-stationary
;      instrumental noise and a Mask:
;      toucan_Cl =  toucan(Map, sigmaNoise=sigmaNoise, Mask=Mask)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------------------------------------

;;============================================== mrs_toucan ===================================================

function mrs_toucan, Map, filters=bkj, var_noise=var_noise, NoiseReal=NoiseReal, sigmaNoise=sigmaNoise, Mask=Mask, Niter=Niter, Nscale=Nscale, Firstl=Firstl, lmax=lmax, Filtwidth=l, lin=lin, tousil=tousil, inpainting=inpainting, scaledMask=scaledMask, rescaleMask=rescaleMask, WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos

COMMON C_PLANCK

Result = -1
if N_PARAMS() LT 1  then begin 
   print, 'CALLING SEQUENCE: Res =  mrs_toucan( Map, filters=bkj, var_noise=var_noise, NoiseReal=NoiseReal, sigmaNoise=sigmaNoise, Mask=Mask, Niter=Niter, Nscale=Nscale, Firstl=Firstl, lmax=lmax, Filtwidth=l, lin=lin, tousil=tousil, inpainting=inpainting, scaledMask=scaledMask, rescaleMask=rescaleMask, WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos)'
   goto, DONE
endif

npix = N_elements(Map)
nside = npix2nside(npix)

;;---------------------------------------- Set default parameters ----------------------------------------
if not keyword_set(Niter) then Niter = fix(10*exp(float(Nscale)/10)) ; Niter = 500 ;
if keyword_set(bkj) then begin
   Nscale = (size(bkj))[2]
   lmax = (size(bkj))[1]-1
endif
if not keyword_set(lmax) then begin
   lmax = long( nside )  * 2l
   if lmax GT P_Lmax then  lmax = P_Lmax  
endif
;; if we don't have many l, we don't want to use the DCT, and we consider more wavelet scales
if lmax lt 1000l then OWT = 1
if not keyword_set(Firstl) then Firstl=0
if not keyword_set(Nscale) then Nscale = round(4*alog(lmax))
if not keyword_set(tousil) then tousil = min([max([100,2*Firstl]), lmax/2.])
if keyword_set(scaledMask) then rescaleMask = 1
if (arg_present(Mask) and not keyword_set(rescaleMask) and not arg_present(inpainting)) then inpainting = 1
if not keyword_set(bkj) then bkj = mrs_toucan_wt_filters(Nscale, lmax, Firstl=Firstl, Filtwidth=l, lin=lin, Bmat=Bmat)

;;------------------------------------------ Inpainting if set ------------------------------------------------------
if keyword_set(inpainting) then begin
   print, '--- inpainting the mask...'
   rescaledMask = 0
   mrs_write,'in_mask.fits',mask
   mrs_write,'in_data.fits',Map
   spawn,'mrs_matmask in_mask.fits out_mat.fits'
   spawn,'mrs_alm_inpainting -m3 -M out_mat.fits  in_data.fits in_mask.fits  out_inp_data.fits'
   Map = rims('out_inp_data.fits')
   print, '--- inpainting done.'
endif

;;------------------- compute the variance of the noise over the wavelet filters if not given -----------------------
if (not keyword_set(var_noise)) then begin
   if keyword_set(sigmaNoise) then var_noise = mrs_var_noise_mc(sigmaNoise, bkj) $
   else if keyword_set(NoiseReal) then var_noise = mrs_var_noise_patch(NoiseReal, bkj) 
endif

;;------------------ Intergrated Power Spectrum over the scales of a wavelet transform -----------------------------
print, '--- computing the compressed power spectrum...'
intCl = mrs_intcl(Map, bkj, var_noise=var_noise, Mask=Mask, scaledMask=scaledMask, rescaleMask=rescaleMask, EdgeMask=maskEdge, cutoffreq=cutoffreq, niter=niter, Bmat=Bmat)
print, '--- compressed power spectrum computed.'

;;------------------------------------------ Reconstructed Cl with Toucan ------------------------------------------
print, 'Compute estimated Power Spectrum with Master to estimate HWT mask...'
if not keyword_set(mask) then Mcl = mrs_powspec(map, lmax=lmax) else Mcl = mrs_master_powspec(map,mask,lmax=lmax)
if ((not keyword_set(NoiseReal)) and keyword_set(SigmaNoise) ) then begin
   NoiseReal = SigmaNoise*randomn(seed, npix)
   pn = mrs_powspec(NoiseReal)
end else pn = 0
powspecInit = mrs_tousi(Mcl[0:2*tousil], NoisePs=pn)
initl = min([tousil, lmax-100])

print, '--- reconstructing the power spectrum...'
uresult = mrs_intcl_to_cl(Map, intCl, bkj, Mask=Mask, Firstl=Firstl, Bmat=Bmat, Niter=Niter, initl=initl, WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos, Mcl=Mcl, Cl=Cl, Proj_Resi=Proj_Resi, True_Resi=True_Resi, part_nMSE=part_nMSE, uresult=uresult)

result = uresult
weight = 1+(exp(5.)-1)*(tousil-findgen(tousil+1))/tousil
weight = alog(weight)
weight2 = (1.-weight/5.)
result[0:tousil] = (weight2*result[0:tousil] + weight*powspecInit[0:tousil])/(weight+weight2) 

DONE:
return, Result

end
