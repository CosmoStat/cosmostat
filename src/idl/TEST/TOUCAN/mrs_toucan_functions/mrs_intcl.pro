;+
; NAME:
;        mrs_intcl
;
; PURPOSE:
;        Computes the power spectrum integrated over filters
;        of a wavelet decomposition.
;
; CALLING:
;     intcl  =  mrs_intcl(Map, bkj, var_noise=var_noise, Mask=Mask,
;     scaledMask=scaledMask, rescaleMask=rescaleMask,
;     EdgeMask=maskEdge, cutoffreq=cutoffreq, niter=niter, Bmat=Bmat)
;
; INPUTS:
;     Map        -- Healpix map: Input CMB image to compute the
;                   compressed measurements from.
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters
;                   in the multipole domain used to compute the
;                   compressed measurements.
;    
; OUTPUTS:
;     intcl     -- IDL dblarr(Nscale) = Compressed measurements
;                  corresponding to the power spectrum Cl integrated
;                  over filters of a wavelet transform.
;
; INPUT/OUTPUT KEYWORD: 
;     var_noise   : array containing Nscale Healpix maps = each map
;                   contains the variance of the noise at that wavelet
;                   scale.
;     Mask        : Healpix map = mask applied to the input Map.
;     ScaledMask  : IDL fltarr(npix, maxaddpix) = contains the max
;                   rescaled to be used if the RescaleMask option is
;                   set.
;     RescaleMask : bool = if set, a rescaling of the mask is
;                   performed at each wavelet scale during the
;                   estimation of the compressed measurements.
;     EdgeMask    : IDL intarr(npix) : map with ones at the locations
;                   of the pixels corresponding to the edge of the
;                   mask to increase the speed of the rescaling.
;     cutoffreq   : IDL fltarr(Nscale+1) = values of the multipoles l
;                   that bound the wavelet filters.
;
; OUTPUT KEYWORDS:
;     Bmat        : IDL dblarr(lmax+1,Nscale) = measurement matrix
;                   constructed from the wavelet filters.
;
; EXAMPLE:
;      Estimate the compressed power spectrum from a CMB map with
;      stationary non-stationary noise and no mask:
;      intcl  =  mrs_intcl(Map, bkj, var_noise=var_noise)
;
;      Estimate the compressed power spectrum from a CMB map with
;      stationary non-stationary noise and mask:
;      intcl  =  mrs_intcl(Map, bkj, var_noise=var_noise, Mask=Mask)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------


;;========================== mrs_intcl =============================

function mrs_intcl, Map, bkj, var_noise=var_noise, Mask=Mask, scaledMask=scaledMask, rescaleMask=rescaleMask, EdgeMask=maskEdge, cutoffreq=cutoffreq, niter=niter, Bmat=Bmat, needt=needt

intcl = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: Res =  mrs_intcl(Map, bkj, var_noise=var_noise, Mask=Mask, scaledMask=scaledMask, rescaleMask=rescaleMask, EdgeMask=maskEdge, cutoffreq=cutoffreq, niter=niter, Bmat=Bmat)'
   goto, DONE
endif

;; ---------------- Needlet transform  
; print, 'Computing Needlet decomposition over', Nscale, ' scales...'
mrs_wt_coeff, Map, needt, bkj
print, 'Needlet decomposition done'

;; ---------------- Parameters of the needlet transform 
npix = needt.Npix
Nscale = needt.NScale
intcl = dblarr(Nscale)

nside = needt.Nside
lmax = needt.lmax ; lmax = 3*nside
ell = findgen(lmax+1)*2+1
bkj = needt.tabpsi

;; ---------------- Default parameters  
if not keyword_set(niter) then niter = 5 ; number of iterations to refine the angular spectrum
if not keyword_set(var_noise) then var_noise = fltarr(npix, Nscale)

;; ------ pre-compute rescaled masks
if (keyword_set(Mask) && keyword_set(rescaleMask) && not (keyword_set(scaledMask)) ) then begin
   print, 'Mask detected'
   maxaddpix = 15
   ;; ------ pre-compute mask edge
   if not (keyword_set(maskEdge)) then begin
      print, 'Computing mask edge...'
      maskEdge = mrs_mask_edge(mask)
      print, 'Mask edge computed'
   endif
   daddpix = indgen(maxaddpix+1)
   print, 'Resizing mask over ', maxaddpix+1, ' sizes...'
   scaledMask = mrs_resize_mask_vec(mask,daddpix,maskEdge=maskEdge)
   print, 'Mask rescaled'
   nmask = total(scaledMask,1)
endif

;; ---------------- Estimate integrated Power Spectrum 

;; ;------ If a mask is present
if keyword_set(Mask) then begin
   
   ;; ;------ If the rescaling of the mask option is present
   if keyword_set(rescaleMask) then begin
      
;;       addpix = round(17*alog(2.-needt.modfreq/lmax))
;;       ind = where(addpix gt maxaddpix, c)
;;       if c ne 0 then addpix(ind) = maxaddpix

;;       addpix = fix(17*alog(2.-needt.cutoffreq[0:Nscale-1]/lmax))
;;       ind = where(addpix gt 1,c)
;;       if c ne 0 then addpix(ind) = maxaddpix

      addpix = round(1000/needt.tabwidth)
      ind = where (cutoffreq/lmax gt 0.8, c)
      if c gt 0 then addpix(ind)=0
      ind = where(addpix gt maxaddpix, c)
      if c gt 0 then addpix(ind) = maxaddpix
      
      ;; ;------ Compute the integrated power spectrum
      for i=0, Nscale-1 do begin
         whichMask = addpix[i]
         intcl[i] = total(scaledMask[*,whichMask]*(needt.coef[*,i]^2-var_noise[*,i]))/nMask[whichMask]
      endfor

;;          ;; ; Refining of the angular Power Spectrum
;;          for iter=1,niter do begin
;;             ;; ; compute optimal weights
;;             for i=0,Nscale-1 do begin
;;                whichMask = ind2element[i]
;;                num = scaledMask[*,whichMask]*(intcl[i] + var_noise[*,i])^(-2)
;;                wkj[*,i] =  num/total(num)
;;                intcl[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
;;             endfor
;;          endfor
   
   endif else begin
      ;; ;------ If the rescaling of the mask option is not present
      nMask = total(total(Mask))
      
      for i=0, Nscale-1 do begin
         intcl[i] = total(Mask*(needt.coef[*,i]^2-var_noise[*,i]))/nMask
      endfor
      
;;          ;; ; Refining of the angular Power Spectrum
;;          for iter=1,niter do begin
;;             ;; ; compute optimal weights
;;             for i=0,Nscale-1 do begin
;;                num = Mask*(intcl[i] + var_noise[*,i])^(-2)
;;                wkj[*,i] =  num/total(num)
;;                intcl[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
;;             endfor
;;          endfor
      
      
   endelse
   
endif else begin
   
;; ;------ If no mask is present
   
   
   for i=0, Nscale-1 do begin
      intcl[i] = total(needt.coef[*,i]^2-var_noise[*,i])/npix
   endfor

;;       ;; ; Refining of the angular Power Spectrum
;;       for iter=1,niter do begin
;;          ;; ; compute optimal weights
;;          for i=0,Nscale-1 do begin
;;             num = (intcl[i] + var_noise[*,i])^(-2)
;;             wkj[*,i] =  num/total(num)
;;             intcl[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
;;          endfor
;;       endfor
   
endelse

if arg_present(Bmat) then begin
;   print, 'computing Bmat'
   Bmat = dblarr(lmax+1,Nscale)
   lcl = findgen(lmax+1)*2+1
   for i=0, Nscale-1 do begin
      Bmat(*,i) = bkj[*,i]^2*lcl
   endfor
   Bmat = Bmat / (4. * !dpi)
endif

DONE:
return, intcl

end

