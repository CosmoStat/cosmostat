;+
; NAME:
;        toucan
;
; PURPOSE:
;  Computes the CMB power spectrum from a CMB map. Takes the map and
;  the noise in input. The noise can be one value, then it is
;  considered to be uniform and stationary. If it is a map, it has the
;  same dimension as the CMB map and the noise is considered to be
;  non-stationary.
;
; CALLING:
;     toucanCl = toucan(Map, sigmaNoise, Mask=Mask, Nscale=Nscale,
;     Firstl=Firstl, lmax=lmax, filters=bkj, var_noise=var_noise,
;     tousil=tousil, scaledMask=scaleMask, WTNscale=WTNscale,
;     Niter=Niter, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos,
;     rescaleMask=rescaleMask, inpainting=inpainting)
;
; INPUTS:
;     Map        -- Healpix map: Input CMB image to estimate the CMB
;                   power spectrum from.
;     SigmaNoise -- Scalar or Healpix map containing the standard
;                   deviation of the noise
;    
; OUTPUTS:
;     toucanCl   -- IDL 1D array  = Estimated Cl array
;
; INPUT KEYWORD: 
;     Mask        : Healpix map = mask applied to the input Map.
;     Nscale      : int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;     Firstl      : int = initial multipole value to be reconstructed
;                   in Cl.
;     lmax        : int = final multipole value to be reconstructed in
;                   Cl.
;     filters     : array of size lmax+1 by Nscale = wavelet filters
;                   in the multipole domain to be used to compute the
;                   compressed measurements.
;     var_noise   : array containing Nscale Healpix maps = each map
;                   contains the variance of the noise at that wavelet
;                   scale.
;     tousil      : int = value of the multipole l below which the
;                   solution is an average between Tousi and Toucan.
;     scaledMask  : array containing Nscale Healpix maps = mask to be
;                   used in the computation of the compressed
;                   measurements at each wavelet scale.
;     WTNscale    : int = number of wavelet scales used to impose
;                   sparsity during reconstruction.
;     Niter       : int = number of iterations for the reconstruction.
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
;     RescaleMask : bool = if set, a rescaling of the mask is
;                   performed at each wavelet scale during the
;                   estimation of the compressed measurements.
;     inpainting  : bool = if set, an impainting of the map is
;                   performed before estimation of the compressed
;                   measurements.
;
; OUTPUT KEYWORD: 
;
; EXAMPLE:
;      Estimate the CMB power spectrum from a CMB map with stationary
;      noise with standard deviation of 5:
;      toucan_Cl =  toucan(Map, 5)
;
;      Estimate the CMB power spectrum from a map with non-stationary
;      instrumental noise:
;      toucan_Cl =  toucan(Map, sigmaNoise)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------------------------------------

;;============================================== Toucan ===================================================

function ntousi, Map, sigmaNoise, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filters=bkj, var_noise=var_noise, tousil=tousil, scaledMask=scaleMask, WTNscale=WTNscale, Niter=Niter, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos, rescaleMask=rescaleMask, inpainting=inpainting

COMMON C_PLANCK

Result = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: Res =  toucan(Map, sigmaNoise, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filters=bkj, var_noise=var_noise, tousil=tousil, scaledMask=scaleMask, WTNscale=WTNscale, Niter=Niter, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos, rescaleMask=rescaleMask, inpainting=inpainting)'
   goto, DONE
endif

npix = N_elements(Map)
nside = npix2nside(npix)

;;---------------------------------------- Set default parameters ----------------------------------------
if not keyword_set(Niter) then Niter = fix(10*exp(float(Nscale)/10)) ; Niter = 500 ;

if not keyword_set(lmax) then begin
   lmax = long( nside )  * 3l
   if lmax GT P_Lmax then  lmax = P_Lmax  
endif
;; if we don't have many l, we don't want to use the DCT, and we consider more wavelet scales
if lmax lt 1000l then OWT = 1

if not keyword_set(WTNscale) then begin
   if not keyword_set(OWT) then  WTNscale = fix(alog(float(lmax)))  - 1 $
   else  WTNscale = fix(alog(float(lmax)))  + 1
endif

if not keyword_set(Firstl) then Firstl=0
if not keyword_set(Nscale) then Nscale = round(4*alog(lmax))
if not keyword_set(tousil) then tousil = min([max([100,2*Firstl]), lmax/2.]) ;=  min([100, lmax]) ******* CHANGE STUFF
if keyword_set(scaledMask) then rescaleMask = 1
if (arg_present(Mask) and not keyword_set(rescaleMask) and not arg_present(inpainting)) then inpainting = 1

if not arg_present(gen2) then gen2 = 1
if not arg_present(SWT) then SWT = 1
if not arg_present(HWT) then HWT = 1

;;------------------------------------------ Inpainting if set ------------------------------------------------------
if keyword_set(inpainting) then begin
   print, 'inpainting the mask...'
   rescaledMask = 0
   mrs_write,'in_mask.fits',mask
   mrs_write,'in_data.fits',Map
   spawn,'mrs_matmask in_mask.fits out_mat.fits'
   spawn,'mrs_alm_inpainting -m3 -M out_mat.fits  in_data.fits in_mask.fits  out_inp_data.fits'
   Map = rims('out_inp_data.fits')
   print, 'inpainting done'
endif

;;------------------ Intergrated Power Spectrum over the scales of a needlet transform -------------------
print, 'computing the compressed power spectrum...'
angSpec = compressed_powspec(Map, sigmaNoise=sigmaNoise, filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, scaledMask=scaledMask, Bmat=Bmat, var_noise=var_noise, rescaleMask=rescaledMask)
print, 'compressed power spectrum computed'

;;------------------------------------------ Print Toucan options ------------------------------------------------------

print,  "SWT = ", SWT
print, "HWT = ", HWT

;;---------------- Normalizing the input --------------------
sigmaSpec = sigma(angSpec)
angSpec = angSpec/sigmaSpec
if keyword_set(Firstl) then Bmat = Bmat[Firstl:*,*]
Phi = Bmat ## transpose(Bmat)
tau = maxeigval(Phi) ; step size for the gradient
sigl = lmax+1 - Firstl ; length of the Cl to recover
;;------------------------------------------------------------

;;------------------ Choose threshold level ------------------
Dirac = dblarr(sigl)
dirac[sigl/2] = 1
w  = uwt1d(dirac, Nscale=WTNscale, Gen2=Gen2, info=info, OptWT=OptWT)
TabPSWTNorm = dblarr(WTNscale)
for j=0, WTNscale-1 do TabPSWTNorm[j] = sqrt( total(w[*,j]^2) )

LevelDetect = fltarr(WTNscale)
Nsigma = total(sigmaNoise)/N_elements(sigmaNoise) ; better way to determine the threhsold at each scale ?
for j = 0, WTNscale-1 do BEGIN
    if j NE 0 then LevelDetect(j) = NSigma * TabPSWTNorm(j)  $
    else LevelDetect(j) = (NSigma+1) * TabPSWTNorm(j)
 endfor
print, 'threshold chosen'
;;------------------------------------------------------------

;;------------ Set the masks for the HWT ----------------
WTMask = fltarr(sigl, WTNscale) 
WTMask[*, WTNscale-1] = 1.
print, 'Compute estimated Power Spectrum with Master to estimate HWT mask...'
if not keyword_set(mask) then Mcl = mrs_powspec(map, lmax=lmax) else Mcl = mrs_master_powspec(map,mask,lmax=lmax)
Signal =double(Mcl[Firstl:*])
SigmaSignal = sigma(Signal)
Signal =  Signal  / SigmaSignal
powspecInit = mrs_tousi(Mcl[0:2*tousil]) 

l = findgen(sigl) + firstl & lfactor = l * (l+1) / 2. / !dpi
if firstl EQ 0 then l[0]=1
if firstl EQ 0 then lfactor[0]=1
LastThreshold=0
FirstThreshold=  max( sqrt(Signal) * lfactor)
DeltaThreshold = FirstThreshold - LastThreshold
DCTDeltaThreshold =  FirstThreshold

indl = 2. * l + 1.
if firstl EQ 0 then indl[0]=1
NormL = (1.+1/indl)
Fa =  -0.25839012       
Fb = 1.1820104
x = findgen(501) / 500.
yy = Fa*x + Fb
ind = where(yy lt 1, c)
if c GT 0 then yy[ind] = 1.
NormL[0:min([sigl-1,500])] = NormL[0:Min([sigl-1,500])] * yy

if keyword_set (HWT)  then begin
  MR_Trans = uwt1d(Signal, gen2=gen2, Nscale=WTNscale, info=info, OptWT=OptWT) ;***
  for j = 0, WTNscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
     Scale     = MR_Trans[*,j]   
     
     ScaleMask = WTMask[*,j]
;;      Level = NormL * LevelDetect[j]
     Level = max([LevelDetect,0.1])
     
     ind       = where(abs(Scale) GT Level , count) ; Sig should be 1 in the GuassData space 
     if count GT 0 then ScaleMask[ind] = 1  
     WTMask[*,j] = ScaleMask  
  END 
  print, 'HWT mask computed'
endif

;;---------------------------------------------------------------

;;------------------ Initialization --------------------------
Resi =  (transpose(Bmat) ## angSpec)/(sqrt(tau)) ; Residual for the gradient step
Resi = Resi[*]
SigmaResi_old = sigma(Resi)
Result = fltarr(sigl) ; Initial guess for the power spectrum
step = 1
StopDecrease = 0
proj_resi = fltarr(Niter)
;;-------------------------------------------------------------

for iter=1, Niter do begin

   if StopDecrease eq 0 then Lambda = ( 1.-erf(2.5*float(iter) / min([Niter,100])) )

   Resi_old = Resi
   Result_old = Result
   step_old = step
   step = (1. + sqrt(1.+4.*step^2))/2.
   Z = Result + (step_old - 1.)*(Result-Result_old)/step
   Result_old = Result
   Resi_old = Resi
   Resi =  (transpose(Bmat) ## (angSpec - Bmat ## Z) )/tau
   Resi = Resi[*]  
   
   ;;========= WT - HARD Thresholding

   if keyword_set (HWT)  then begin
      MR_Trans = uwt1d(Resi, gen2=gen2, Nscale=WTNscale, info=info, OptWT=OptWT)
      if iter EQ 1 then MaxHWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
         Scale     = MR_Trans[*,j]   
         ScaleMask = WTMask[*,j]
         MR_Trans[*,j] = Scale * ScaleMask
      END 
      Result =  Result  + iuwt1d(MR_Trans, gen2=gen2, info=info)
   endif
   
 
   ;;========= WT - SOFT Thresholding

   if keyword_set(SWT) then begin  
      MR_Trans = star1d(Result*lfactor, gen2=gen2, Nscale=WTNscale)
      if iter EQ 1 then MaxWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do begin 
         SoftThreshold =   MaxWT  *  Lambda *  TabPSWTNorm[j] 
         Scale  = reform(MR_Trans[*,j]) 
         softthreshold, Scale, SoftThreshold
         MR_Trans[*,j] = Scale 
      endfor
      Result = istar1d(MR_Trans, gen2=gen2)/lfactor
   endif
   
   ;;========= Positivity constrain
   if keyword_set(Pos) then begin
      ind = where(Result LT 0, c)
      if c GT 0 then Result[ind] = 0 
   endif
   
   sigmaResi = sigma(Resi)
   proj_resi(iter-1) = sigmaResi
   
   print, Iter, '/', Niter, ':  Lambda = ', Lambda, ' | SigmaResi = ', sigma(Resi), ' | SigmaResult = ', sigma(Result)    
;;    print, Iter, ':  Linf(ResiDiff)= ', max(abs(Resi-Resi_old))

   if sigmaResi_old lt sigmaResi-(1.e-8) then begin
      print, 'Increasing residual -> iteration interrupted ==> try fewer iterations'
;;       interupted = 1
;;       goto, DONE
   endif
   sigmaResi_old = sigmaResi

endfor

DONE:
result = [fltarr(Firstl),result*sigmaSpec]
uresult = result
weight = 1+(exp(5.)-1)*(tousil-findgen(tousil+1))/tousil
weight = alog(weight)
weight2 = (1.-weight/5.)
result[0:tousil] = (weight2*result[0:tousil] + weight*powspecInit[0:tousil])/(weight+weight2) 

return, Result

end
