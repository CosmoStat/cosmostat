;+
; NAME:
;        mrs_cl_to_intcl
;
; PURPOSE:
;  Computes the integration of the Power Spectrum Cl over several
;  wavelet filters.
;
; CALLING:
;     intCl = mrs_cl_to_intcl(cl, bkj)
;
; INPUTS:
;     Cl         -- fltarr(lmax+1) = true power Spectrum Cl.
;                   power spectrum from.
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;                   Cl.
;    
; OUTPUTS:
;     intCl      -- fltarr(Nscale) = Power Spectrum Cl integrated over the wavelet
;                   filters bkj.
;
; EXAMPLE:
;      Computes intcl over filters bkj from Cl:
;      intCl = mrs_cl_to_intcl(cl, bkj)
;
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------------------------------------

;;============================================== mrs_cl_to_intcl ===================================================

function mrs_cl_to_intcl, cl, bkj

Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1
ell = findgen(lmax+1)*2.+1.
ellcl = ell*cl[0:lmax]

intcl = fltarr(Nscale)
for i=0,Nscale-1 do begin
   intcl(i) = total(bkj[*,i]^2*ellcl)/(4*!dpi)
endfor
