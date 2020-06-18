;+
; NAME:
;        mrs_iter_wiener
;
; PURPOSE:
;	Iterative Wiener denoising taking into account non stationary noise and missing data.
;       The power spectrum of the noise-free signal is assumed to be known (Cl keyword), 
;       but an estimation is done if it is not given to the routine.
;       Non stationary Gaussian noise can be considered using the keyword RMSNoiseMap,
;       or the RMS map can be estimated if a noise realization map is given using the 
;       the keyword NoiseReaMap.
;       If RMSNoiseMap and NoiseReaMap are not given, a white Gaussian noise is assume.
;       If Cl is not given or if the keyword EstimatCL is set, the power spectrum Cl is
;       estimated from the data
;
; CALLING:
;
;     mrs_iter_wiener, Imag, Filter, RMSNoiseMap=RMSNoiseMap, Mask=Mask,
;                   SigmaNoise=SigmaNoise, Niter=Niter, Cl=Cl, lmax=lmax, 
;                   EstimatCL=EstimatCL, NoiseReaMap=NoiseReaMap
;    
; INPUT:
;     Imag -- IDL array of healpix map: Input image be filtered 
;
; OUTPUT:
;     Filter -- IDL array of healpix map: reconstructed denoising image   
;
; INPUT KEYWORDS:
;      Cl : int array = Number of scales (default is 4)
;      NSigma: float = Level of thresholding (default is 3)
;      SigmaNoise: float = If set, white Gaussian noise is assumed with a   
;                          Noise standard deviation equal to SigmaNoise.
;                          Default is automatically estimated
;      RMSNoiseMap: IDL array of healpix map: noise RMS per pixel.
;                   if set, the noise is assumed to be Gaussian but non stationary.
;       NoiseReaMap: IDL array of healpix map: noise realization map
;      Mask: healpix map: if set, missing data are considered. 
;                         Mask[k] = 0 for missing data, and 1 for good data
;      niter: int: number of iterations used in the filtering. Default is 20.
;      EstimatedCl: scalar  -- if set, Cl power spectrum is estimate from the data
;
; INPUT/OUTPUT KEYWORDS:
;      EstimatedCl --  : Threshold wavelet decomposition of the input image
;
; EXTERNAL CALLS:
;       
;
; EXAMPLE:
;        
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2008
;-
;==================================================================



pro mrs_iter_wiener, Imag, Filter, RMSNoiseMap=RMSNoiseMap, Mask=Mask, $ 
                   SigmaNoise=SigmaNoise, Niter=Niter, Cl=Cl, lmax=lmax, $
                   EstimatCL=EstimatCL, NoiseReaMap=NoiseReaMap, WienerNoiseRea=WienerNoiseRea


COMMON C_PLANCK
 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_iter_wiener, Imag, Filter, RMSNoiseMap=RMSNoiseMap, Mask=Mask,  SigmaNoise=SigmaNoise, Niter=Niter, Cl=Cl, lmax=lmax,  EstimatCL=EstimatCL, NoiseReaMap=NoiseReaMap '
        goto, DONE
        end

if not keyword_set(Niter) then Niter=20
 
if keyword_set (mask) then OutMask = OutMask*mask
npix = (size(imag))[1]
nside = npix2nside(npix)


if not keyword_set(lmax) then begin
   lmax = long( nside )  * 3l
   if lmax GT P_Lmax then  lmax =  P_Lmax
end

; White Gaussian noise
if not keyword_set(Mask) then begin
  Mask = fltarr(npix)
  Mask[*] = 1
end
IndMaskOK = where( Mask EQ 1)
Filter = fltarr(npix)

; RMS noise map creation
if not keyword_set(RMSNoiseMap) then  BEGIN 
   ; Automatic white Gaussian estimation (this has to be improved!)
   if keyword_set(NoiseReaMap) then begin
     print, "NEED to estimated the RMSNoiseMap from NoiseReaMap "
   end else begin  ; Here it is white Gaussian noise
      if not keyword_set(SigmaNoise) then SigmaNoise = sigma(Imag)
      RMSNoiseMap = fltarr(npix)
     RMSNoiseMap[*] = SigmaNoise
   end
END
MaxRms = max(RMSNoiseMap[IndMaskOK])
MinRms = max(RMSNoiseMap[IndMaskOK])
Ind = where( RMSNoiseMap NE 0, c)
if c GT 0 then MinNonZeroRMSMap = min( RMSNoiseMap[Ind] ) else MinNonZeroRMSMap = 1.
print, "RMSmap: Min = ", MinRMS, ", Max = ", MaxRms, " MinNZ = ", MinNonZeroRMSMap
; Question: que fait-on si RMSMap[i] = 0 ?

Resi = Imag
Resi = Resi * Mask
 NoiseSpectrum = fltarr(lmax+1)

if keyword_set(NoiseReaMap) then begin
   OldFilter=0
   NoiseFilter =  fltarr(npix)
   NoiseResi =  fltarr(npix)
   WienerNoiseRea =  fltarr(npix)
end

SigmaData = sigma(Imag)
DeltaSig=0.
for i=0, Niter-1 do begin
    ; Normalized Residual Computation
   Resi[*] = 0
   Resi[IndMaskOK] = (Imag[IndMaskOK] - Filter[IndMaskOK]) / RMSNoiseMap[IndMaskOK] * MinNonZeroRMSMap
   NoiseSpectrum[1:*] = MinNonZeroRMSMap^2
   
   if keyword_set(NoiseReaMap) then NoiseResi[IndMaskOK] = (NoiseReaMap[IndMaskOK] - WienerNoiseRea[IndMaskOK]) / RMSNoiseMap[IndMaskOK] * MinNonZeroRMSMap
 
   ; NoiseSpectrum must be defined first: it is either the variance of the noise
   ; or the noise power spectrum
   ; Cl is the power spectrum of the data.
   mrs_wiener, Resi, NoiseSpectrum, Recons, lmax=lmax, WienerFilter=WienerFilter, /Cole, NoiseRea=NoiseResi,  WienerNoiseRea=WienerNoiseResi
   Filter = Filter + Recons
   
   if keyword_set(NoiseReaMap) then begin
     WienerNoiseRea = WienerNoiseRea + WienerNoiseResi
     mrs_almtrans, WienerNoiseRea, AlmNoise, /tab, /complex, /norm,lmax=lmax 
     NoiseSpectrum = mrs_alm2spec(AlmNoise)
     mrs_wiener, Filter, NoiseSpectrum, Filter2, lmax=lmax, WienerFilter=WienerFilter2, /Cole, NoiseRea=WienerNoiseRea,  WienerNoiseRea=WienerNoiseFilter
     Filter = Filter2
     DeltaSig = sigma(OldFilter - Filter)
     OldFilter = Filter
     WienerNoiseRea = WienerNoiseFilter
  ;    tvs, WienerNoiseFilter, tit='JK '+STRC(i+1)
 ;     tvs, Filter, tit='Sol '+STRC(i+1)
   end
   ; Wiener filtering of Resi
      
     print, "Iter ", i+1, ": SigmaResi = ", sigma(Resi), ", SigmaSol = ", sigma(Filter), ", DeltaSig = ", DeltaSig / SigmaData * 100.

  end

; if keyword_set(EstimatCL) then 

DONE:

END



