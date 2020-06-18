;+
; NAME:
;        compressed_powspec
;
; PURPOSE:
;        Computes the compressed power spectrum. The measurements
;        correspond to the true power spectrum integrated over filters
;        of a wavelet decomposition.
;
; CALLING:
;     ang_ps  =  compressed_powspec(Map, sigmaNoise=sigmaNoise,
;     filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax,
;     scaledMask=scaledMask, Bmat=Bmat, rescaleMask=rescaleMask,
;     EdgeMask=maskEdge, niter=niter,
;     var_noise=var_noise, cutoffreq=cutoffreq)

;
; INPUTS:
;     Map        -- Healpix map: Input CMB image to compute the
;                   compressed measurements from.
;    
; OUTPUTS:
;     ang_ps     -- IDL dblarr(Nscale) = Compressed measurements.
;
; INPUT/OUTPUT KEYWORD: 
;     sigmaNoise  : Healpix map = map of the noise standard
;                   deviation in the image domain.
;     filters     : array of size lmax+1 by Nscale = wavelet filters
;                   in the multipole domain to be used to compute the
;                   compressed measurements.
;     Mask        : Healpix map = mask applied to the input Map.
;     Nscale      : int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;     Firstl      : int = initial multipole value to be reconstructed
;                   in Cl.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     lmax        : int = final multipole value to be reconstructed in
;                   Cl.
;     ScaledMask  : IDL fltarr(npix, maxaddpix) = contains the max
;                   rescaled to be used if the RescaleMask option is
;                   set.
;     Bmat        : IDL dblarr(lmax+1,Nscale) = measurement matrix
;                   constructed from the wavelet filters.
;     RescaleMask : bool = if set, a rescaling of the mask is
;                   performed at each wavelet scale during the
;                   estimation of the compressed measurements.
;     EdgeMask    : IDL intarr(npix) : map with ones at the locations
;                   of the pixels corresponding to the edge of the
;                   mask to increase the speed of the rescaling.
;     var_noise   : array containing Nscale Healpix maps = each map
;                   contains the variance of the noise at that wavelet
;                   scale.
;     cutoffreq   : IDL fltarr(Nscale+1) = values of the multipoles l
;                   that bound the wavelet filters.
;
; EXAMPLE:
;      Estimate the compressed power spectrum from a CMB map with
;      stationary non-stationary noise and no mask:
;      ang_ps  =  compressed_powspec(Map, sigmaNoise=sigmaNoise)
;
;      Estimate the compressed power spectrum from a CMB map with
;      stationary non-stationary noise and mask:
;      ang_ps  =  compressed_powspec(Map, sigmaNoise=sigmaNoise,
;      Mask=Mask)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------


;;========================== compressed_powspec =============================

function compressed_powspec, Map, sigmaNoise=sigmaNoise, filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtwidth=filtwidth, scaledMask=scaledMask, Bmat=Bmat, rescaleMask=rescaleMask, EdgeMask=maskEdge, niter=niter, var_noise=var_noise, cutoffreq=cutoffreq, needt=needt

ang_ps = -1
if N_PARAMS() LT 1  then begin 
   print, 'CALLING SEQUENCE: Res =  compressed_powspec(Map, sigmaNoise=sigmaNoise, filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, scaledMask=scaledMask, Bmat=Bmat, rescaleMask=rescaleMask, EdgeMask=maskEdge, niter=niter, var_noise=var_noise, cutoffreq=cutoffreq)'
   goto, DONE
endif

;; ---------------- Needlet transform  
if not keyword_set(needt) then begin
   print, 'Computing Needlet decomposition over', Nscale, ' scales...'
   ntousi_needlet, Map, needt, NScale=Nscale, Firstl=Firstl, lmax=lmax, filters=bkj, filtwidth=filtwidth
   print, 'Needlet decomposition done'
endif

;; ---------------- Parameters of the needlet transform 
npix = needt.Npix
Nscale = needt.NScale
nside = needt.Nside
lmax = needt.lmax ; lmax = 3*nside
ell = findgen(lmax+1)*2+1
bkj = needt.tabpsi
cutoffreq = needt.cutoffreq

ang_ps = dblarr(Nscale)

;; ---------------- Default parameters  
if not keyword_set(sigmaNoise) then sigmaNoise = 7
if not keyword_set(niter) then niter = 5 ; number of iterations to refine the angular spectrum
;; if (keyword_set(Mask) && not (arg_present(rescaleMask)) ) then rescaleMask = 1

;; ------ pre-compute rescaled masks
if (keyword_set(Mask) && keyword_set(rescaleMask) && not (keyword_set(scaledMask)) ) then begin
   print, 'Mask detected'
   maxaddpix = 15
   filename = 'scaledmask_'+strcompress(Nside,/remove_all)+'_'+strcompress(maxaddpix,/remove_all)+'.save'
   filepresent = file_search(filename,count=c)
   if c eq 0 then begin
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
      save, scaledMask, filename=filename
   endif else begin
      restore, filename
   endelse
   nmask = total(scaledMask,1)
endif

;; ;------ I a map of the variances at each scale is not given then it
;;         is computed from the map sigmaNoise
if not keyword_set(var_noise) then var_noise = estime_var_noise(sigmaNoise, bkj, Nscale=Nscale, npix=npix, lmax=lmax)

;; ---------------- Estimate angular Power Spectrum 

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

;;       addpix = round(1000/needt.tabwidth)
;;       addpix = round(1000/(filtwidth+fltarr(Nscale)))
;;       ind = where (cutoffreq/lmax gt 0.8, c)
;;       if c gt 0 then addpix(ind)=0
;;       ind = where(addpix gt maxaddpix, c)
;;       if c gt 0 then addpix(ind) = maxaddpix
      
      ;; ;------ Compute the integrated power spectrum
      for i=0, Nscale-1 do begin
         whichMask = addpix[i]
         ang_ps[i] = total(scaledMask[*,whichMask]*(needt.coef[*,i]^2-var_noise[*,i]))/nMask[whichMask]
      endfor

;;          ;; ; Refining of the angular Power Spectrum
;;          for iter=1,niter do begin
;;             ;; ; compute optimal weights
;;             for i=0,Nscale-1 do begin
;;                whichMask = ind2element[i]
;;                num = scaledMask[*,whichMask]*(ang_ps[i] + var_noise[*,i])^(-2)
;;                wkj[*,i] =  num/total(num)
;;                ang_ps[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
;;             endfor
;;          endfor
   
   endif else begin
      
      ;; ;------ If the rescaling of the mask option is not present
      nMask = total(total(Mask))
      
      for i=0, Nscale-1 do begin
         ang_ps[i] = total(Mask*(needt.coef[*,i]^2-var_noise[*,i]))/nMask
      endfor
      
;;          ;; ; Refining of the angular Power Spectrum
;;          for iter=1,niter do begin
;;             ;; ; compute optimal weights
;;             for i=0,Nscale-1 do begin
;;                num = Mask*(ang_ps[i] + var_noise[*,i])^(-2)
;;                wkj[*,i] =  num/total(num)
;;                ang_ps[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
;;             endfor
;;          endfor
      
      
   endelse
   
endif else begin
   
;; ;------ If no mask is present
   
   
   for i=0, Nscale-1 do begin
      ang_ps[i] = total(needt.coef[*,i]^2-var_noise[*,i])/npix
   endfor

;;       ;; ; Refining of the angular Power Spectrum
;;       for iter=1,niter do begin
;;          ;; ; compute optimal weights
;;          for i=0,Nscale-1 do begin
;;             num = (ang_ps[i] + var_noise[*,i])^(-2)
;;             wkj[*,i] =  num/total(num)
;;             ang_ps[i] = total(wkj[*,i]*(needt.coef[*,i]^2-var_noise[*,i]))
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
return, ang_ps

end

;;  ; example
;; test = [0,0,0,2,2,2,4,4,8]
;; elementsOf, test, elements=elements, ind2element=ind2element, Nelements=Nelements
;;  ; the result is:
;;  ; elements = [0,2,4,8]
;;  ; ind2element = [0,0,0,1,1,1,2,2,3]
;;  ; Nelements = 4
