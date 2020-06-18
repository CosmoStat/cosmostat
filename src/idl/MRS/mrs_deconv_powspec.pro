;+
; NAME:
;        mrs_deconv_powspec
;
; PURPOSE:
;   Computes the power spectrum of a map with missing data, 
;   using the HEALPix representation (nested data
;   representation by default). The MASTER method is used, but using an iterative technique instead
;   of a brute force matrix inversion. If a noise power spectrum is given, it is removed from the solution.
;   The noise power spectrum is assumed to be unbiased, i.e. free of missing area.
;  
; CALLING:
;     P = mrs_deconv_powspec( PowSpecMaskedData,  Mask,  MatMask=MatMask,  Niter=Niter,  NoisePS=NoisePS, lmax=lmax)
;
; INPUTS:
;     PowSpecMaskedData -- IDL array : Power spectrum of the masked data   
;     Mask -- IDL array of healpix map: Input Mask of missing data (Mask[k] =0, if pixel k is missing) 
;     
; OUTPUTS:
;     P -- 1D IDL fltarr:  Power Spectrum corrected from the mask,
;
; OPTIONAL INPUT KEYWORDS:
;     PSCMB_Mask -- IDL array:   power  spectrum of the mask
;     Niter -- int: Number of iteration. Default is 20.
;     NoisePS -- IDL 1D array: Noise power spectrum. Default is none.
;
; INPUT/OUTPUT KEYWORDS:
;      MatMask -- IDL 2D array: Master matrix related to the mask (see mrs_matmask.pro)
;     
; EXTERNAL CALLS:
;       mrs_matmask.pro
;
; EXAMPLE:
;       Deconvolve the power spectrum of an image. 
;               P = mrs_deconv_powspec(PowSpecMaskedData, Mask) 
;        
;       or directly from the image
;               P = mrs_deconv_powspec(HealpixImage, Mask, /map) 
  
; HISTORY:
;	Written: Jean-Luc Starck 
;	October, 2011 File creation
;       June 16, 2014 Call C++ software call
;--------------------------------------------------------------------------------------------------------

function mrs_deconv_powspec,  PData,  Mask,  MatMask=MatMask, PRes, Niter=Niter, TrueCl=TrueCL, NoisePS=NoisePS, lmax=lmax, mu=mu, zeromonodip=zeromonodip, verb=verb, fsky=fsky, inv=inv, ImatMask=IMatmask, IdealBeam=IdealBeam,  IDL=IDL, map=map,noPositivity=noPositivity
COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: p = mrs_deconv_powspec( PowSpecMaskedData,  Mask,  MatMask=MatMask, niter=niter, NoisePS=NoisePS, lmax=lmax, map=map)'
	    Pres =-1
        goto, DONE
        end
        
if keyword_set(HealpixCXX) and keyword_set(map) and not keyword_set(IDL) then begin
        FN1 = gettmpfilename()
        FN2 = gettmpfilename()
        FN3 = gettmpfilename()
        FN4 = gettmpfilename()

        mrs_write, FN1, PData
        cmd = BIN_ISAPCXX + '/mrs_powspec  ' 
        mrs_write, FN2, Mask
        cmd = cmd + ' -m ' + FN2

z = rim(FN1)
z = rims(FN2)

       if keyword_set(lmax) then cmd = cmd + ' -l ' + STRC(lmax)
       if keyword_set(MatMask) then begin
            writefits, FN3, MatMask
            cmd = cmd + ' -M ' + FN3
       end
      if keyword_set(Niter) then cmd = cmd + ' -i ' + STRC(Niter)
       cmd = cmd + ' ' +  FN1 + ' ' + FN4
;        print, cmd
        spawn, cmd
        Pres = readfits(FN4)
       delete, FN1
       delete, FN2
       delete, FN3
       delete, FN4
 end else  begin
 if keyword_set(map) then PowSpecMaskedData = mrs_powspec(PData) $
else PowSpecMaskedData = PData
PDat = PowSpecMaskedData

if keyword_set(lmax) then begin
    vs = size(PowSpecMaskedData)
    Nl = vs[1] +1 
    if lmax lt Nl then PDat = Pdat[0:lmax]
    PowSpecMaskedData = PowSpecMaskedData[0:lmax]
end else begin
 vs = size(PowSpecMaskedData)
 lmax = vs[1] - 1
end

if not keyword_set(NoisePS) then NoisePS=0.

; Brute force Inversion of the matrix
if keyword_set(inv) or keyword_set(IMatmask) then begin
   if not keyword_set(IMatmask) then begin
       if not keyword_set(MatMask) then MatMask = mrs_matmask(Mask, lmax=lmax)
       IMatmask = invert (matmask)
   end
   PRes = IMatmask # PowSpecMaskedData
   if keyword_set(NoisePS) then PRes = PRes - IMatmask # NoisePS
    ; Positivity constraint
   if not keyword_set(noPositivity) then begin
       ind = where(Pres LT 0, c)
       if c GT 0 then Pres[ind] = 0
   endif
end else begin

if not keyword_set(niter) then niter=40
if not keyword_set(MatMask) then MatMask = mrs_matmask(Mask, lmax=lmax)
if not keyword_set(mu) then mu = 1.


if not keyword_set(fsky) and keyword_set(Mask) then Fsky = float( N_ELEMENTS(Mask) / float ( total(mask)) ) 
if not keyword_set(fsky) then Fsky=1

; print, fsky

; Firt initialisation: Fsky correction and noise subtraction
PRes = PDat*Fsky - NoisePS
;Energy = total(Pres)
;print, ", Energ = ", total(Pres)

Tmat = transpose(MatMask)
; mu = 1./ fsky^2
if keyword_set(IdealBeam) then begin
	vs = size(PDat)
	Nl = vs[1]
	LowPassFilter = dblarr(Nl)
	vs = size(IdealBeam)
	Nl1 = vs[1]
	if NL1 gt Nl then LowPassFilter = IdealBeam[0:Nl-1] $
	else LowPassFilter[0:Nl1-1] = IdealBeam
end

for i=0,Niter-1 do begin
   Resi = PDat - MatMask # ( PRes + NoisePS )
   if keyword_set(zeromonodip) then PRes[0:1] = 0
   
   if keyword_set(IdealBeam) then  PRes += mu * (Tmat # resi) * LowPassFilter  $
   else PRes += mu * (Tmat # resi) 
   
   if  keyword_set(verb) then print, i+1, ': Resi = ', sigma(Resi), ", Energ = ", total(Pres)
   ; Positivity constraint
   if not keyword_set(noPositivity) then begin
       ind = where(Pres LT 0, c)
       if c GT 0 then Pres[ind] = 0
   endif
   ; Pres = Pres * Energy / total(Pres)
   if keyword_set(TrueCl) then begin
      plotcl, TrueCl, thick=2
      oplotcl, Pres
   end
end
end 
end
 
DONE:
 
  RETURN, Pres
end


