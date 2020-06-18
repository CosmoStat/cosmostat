;+
; NAME:
;        mrs_whitening
;
; PURPOSE:
;     Apply a whitening to a in image, i.e. it divides its spherical harmonics coefficients 
;     by the sqrt root of its power spectrum:
;          Alm = Alm / sqrt( Cl)
;     If the power spectrum is not given, it is computed on the map.
;  
;     If the keyword INV is set and the power spectrum Cl is given, then the routine performs
;     the inverse whitening transform:
;            Alm = Alm * sqrt( Cl)
;
; CALLING:
;     W_Ima = mrs_whitening(Imag,  Cl=Cl, inv=inv, ALM=ALM, lmax=lmax)
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be whitened 
;    
; OUTPUTS:
;     W_Ima -- IDL array of healpix map: Output whitened image
;
; INPUT KEYWORDS:
;     Cl -- 1D IDL fltarr: Power Spectrum. P[k] = Mean(  POWSPECTRUM[*,l]  )
;     ALM -- ALM structure (see mrs_almtrans.pro): if set, the spherical harmonic are not computed, 
;                           and the input variable Imag is not used.
;     lmax -- int: maximum multipole used in the spherical harmonic decomposition
;     inv -- scalar: if set, apply the inverse whitening.
;
; OUTPUT KEYWORDS:
;
; EXTERNAL CALLS:
;       mrs_almtrans.pro,  mrs_alm2spec.pro
;
; EXAMPLE:
;       Apply the whitening to an image. 
;               P = mrs_whitening(Imag) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2009
;	April, 2009 File creation
;-
;--------------------------------------------------------------------------------------------------------

 
function mrs_whitening, Ima, Cl=Cl, inv=inv, ALM=ALM, lmax=lmax, HC=HC

if N_PARAMS() LT 1 and not keyword_set(ALM)  then begin 
        print, 'CALLING SEQUENCE:  mrs_whitening, Ima, Cl=Cl, inv=inv, ALM=ALM, lmax=lmax'
	Rec =-1
        goto, DONE
        end

if not keyword_set(Cl) and keyword_set(inv)  then begin
        print, 'Error: CL keyword must be set when the INV keyword is used ....'
	Rec =-1
        goto, DONE
        end

if not keyword_set(ALM) then begin
  if keyword_set(inv) then mrs_almtrans, IMA, Alm, /tab, lmax=lmax, /norm $
  else mrs_almtrans, IMA, Alm, /tab, lmax=lmax
end
 ; info, Alm.alm

if not keyword_set(Cl) then Cl = mrs_alm2spec(Alm)

Coef = sqrt(Cl[2:*])
ind = where(Coef EQ 0, c)
if c GT 0 then Coef[ind] = 1.
if not keyword_set(inv) then Coef = 1. / Coef
for l=2, ALM.LMAX do ALM.alm[l,*,*] = Alm.alm[l,*,*] * Coef[l-2]
ALM.alm[0:1,*,*] = 0
if keyword_set(inv) then ALM.norm =0 $
else  ALM.norm =  1

mrs_almrec, ALM, Rec
ALM.norm = 0

if keyword_set(HC) then  begin
   TabStat =  get_stat(Rec, HCIma=HCIma)
   Rec = HCIma
end

DONE:
return, Rec
END
