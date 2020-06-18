;+
; NAME:
;        mrs_crosspowspec
;
; PURPOSE:
;   Computes the cross power spectrum of two maps,  using the HEALPix representation (NESTED data
;   representation by default). 
;
; CALLING:
;     P = mrs_crosspowspec(Map1, Map2, a1=a1, a2=a2, StdPS=StdPS, lmax=lmax)
;
; INPUTS:
;     Map1 -- IDL array of healpix map: First Input image  
;     Map2 -- IDL array of healpix map: Second Input image 
;    
; OUTPUTS:
;     P -- 1D IDL fltarr: Cross Power Spectrum. P[k] = Mean(  Alm[Map1] Alm*[Map2]  )
;
; INPUT/OUTPUT KEYWORDS:
;		a1  -- IDL Alm structure (see mrs_almtrans):  Alm of the image Map1. If this keyword is set, then Map1 is not used.
;		a2  -- IDL Alm structure(see mrs_almtrans):  Alm of the image Map2. If this keyword is set, then Map2 is not used.
;
; OUTPUT KEYWORDS:
;		StdPS -- 1D IDL float array: contains the estimated standard deviation on the compute cross power spectrum.
;
; EXTERNAL CALLS:
;       mrs_almtrans.pro
;
; EXAMPLE:
;       Compute the cross power spectrum of two images. 
;               P = mrs_crosspowspec( Imag ) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2010
;	June, 2010 File creation
;--------------------------------------------------------------------------------------------------------
 
;=====================================================================
 
function mrs_crosspowspec, I1, I2, a1=a1, a2=a2, StdPS=StdPS, lmax=lmax

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: p = mrs_crosspowspec(Map1, Map2, a1=a1, a2=a2, StdPS=StdPS, lmax=lmax)'
	Ret=-1
        goto, DONE
        end

if not keyword_set(a1) then mrs_almtrans, i1, a1, /tab, /complex, lmax=lmax
if not keyword_set(a2) then mrs_almtrans, i2, a2, /tab, /complex, lmax=lmax

NbrAlm = a1.LMAX+1
TabM = a1.TABNBRM
C1 = a1.ALM[*,*]
C2 =  a2.ALM[*,*]
StdPS = fltarr(NbrAlm)
Ret = fltarr(NbrAlm)

CrossPSIma = float(C1) *  float(C2) + imaginary(C1) *  imaginary(C2)
for i=0, NbrAlm-1 do  begin
   allalm = 2*TabM[i]-1
   if allalm GT 1 then begin
      P = fltarr(allalm)
      P[0:TabM[i]-2] = reverse( CrossPSIma[i, 1:TabM[i]-1])
      P[TabM[i]-1:*] = CrossPSIma[i, 0:TabM[i]-1]
      Ret[i] = mean( P)  
     StdPS[i] = sigma(P)
   end else begin
       Ret[i] = CrossPSIma[i, 0]
       StdPS[i] = 0.
   end
end

DONE:
return, Ret
end


