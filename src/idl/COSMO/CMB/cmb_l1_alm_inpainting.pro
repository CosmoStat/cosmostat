;+
; NAME: 
;       CMB_L1_ALM_INPAINTING
;
; PURPOSE:
;        Apply an inpainting to a spherical map using the l_1 minimization of the spherical harmonics coefficients.
;        If a mask is not provided, all pixels with a zero value are considered as missing pixels.
;       If file names are given for the input image and mask, these two images are not loaded into IDL.
;
; CALLING:
;     InpaintMap = cmb_l1_alm_inpainting(Imag, Mask=Mask, Niter=Niter, OutPowSpec=OutPowSpec, lmax=lmax, gauss=gauss)
;
; INPUTS:
;     Imag -- IDL 1D array: Input Healpix image to be inpainted 
;     
; OUTPUTS:
;     InpaintMap -- IDL 1D array: Output inpainted Healpix image   
;          
; INPUT KEYWORDS:
;      niter: int: number of iterations used in the reconstruction
;      FNin: Filename containing the input image. The input image won't be read.
;      FNMask: Filename containing the mask. The mask won't be read.
;      FNOut: Filename containing the results. By default, nothing is written on the disk.
;      FNPS: Filename containing the mask. The mask won't be read.
;
; OUTPUT KEYWORDS:
;     OutPowSpec: IDL 1D array: Cl of the inpainted map. 
;       
; EXTERNAL CALLS:
;      mrs_mca
;
; EXAMPLE:
;      
; HISTORY:
;       Written : Jean-Luc Starck   2009.
;-
;-----------------------------------------------------------------

function cmb_l1_alm_inpainting, Imag, Mask=Mask,  niter=niter, to_glesp=to_glesp, FNImag=FNImag, FNMask=FNMask, FNOut=FNOut, FirstIterConstr=FirstIterConstr, OutPowSpec=OutPowSpec
COMMON C_PLANCK

DEF_ALM_FAST = 0
DEF_ALM_NITER = 0

if not keyword_set(niter) then niter=40
 
  if keyword_set(FNMask) then Mask = mrs_read(FNMask)
  if keyword_set(FNImag) then Imag = mrs_read(FNImag)

  if not keyword_set(Mask)  then BEGIN
      Mask = Imag
      if  type_code(Imag) EQ 8 then begin
          Mask.t_sky[*] = 0
          ind = where(Imag.t_sky NE 0, c)
         if c GT 0 then Mask.t_sky[ind] = 1
      end else begin
         ind = where(Imag NE 0, c)
         if c GT 0 then Mask[ind] = 1
     end
    END
 
 ; We remove the mean value from the available input data
 if  type_code(Imag) NE 8 then BEGIN
  ind = where(mask NE 1, c)
  if c GT 0 then mask[ind] = 0
  indMask = where( mask NE 0, cMask)
  Map = Imag
  MeanVal = mean(map[indMask])
  if cMask NE 0 then map[indMask] = map[indMask] -  MeanVal
  map = map * mask
END

if keyword_set(to_glesp) then begin
  GImag = healpix2glesp(Map, /alm)
  GMask = healpix2glesp(mask, optnx=GImag.nx, optnp=GImag.np)
  mrs_mca, GImag, ImagOut, mask=GMask,  niter=niter , residual=residual, selecttrans=[4],/expo,/cstsigma;/fit, , FirstIterConstr=FirstIterConstr    ;,/cstsigma;,/disp;/expo;,/disp;,/cstsigma,/disp
end else mrs_mca, Map, ImagOut, mask=Mask,  niter=niter, residual=residual, selecttrans=[4],/expo ,/cstsigma , /nomean, FirstIterConstr=FirstIterConstr

if  type_code(Imag) NE 8 then BEGIN
if cMask NE 0 then ImagOut[indMask] = ImagOut[indMask] +  MeanVal
END

if keyword_set(FNOut) then mrs_write, FNOut, ImagOut

OutPowSpec = mrs_powspec(ImagOut, /nono)

return, ImagOut
end


