;+
; NAME:
;        mrs_compute_bmat
;
; PURPOSE:
;  Computes the sensing matrix that returns the integrated Cl when
;  multiplying the Cl.
;
; CALLING:
;     mrs_compute_Bmat(bkj)
;
; INPUTS:
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain to used to compute the integrated
;                   Cl.
;    
; OUTPUTS:
;     Bmat       -- dblarr(lmax+1,Nscale) = matrix formed from the
;                   filters. When multiplying Cl by Bmat, one obtains
;                   the integrated Cl.
;
; EXAMPLE:
;      Compute the sensing matrix:
;      Bmat = mrs_compute_Bmat(bkj)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;--------------------------------------------------------------------------------------------------------

function mrs_compute_Bmat, bkj

Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1

print, '--- Compute sensing matrice Bmat...'

Bmat = dblarr(lmax+1,Nscale)
lcl = findgen(lmax+1)*2+1
for i=0, Nscale-1 do begin
   Bmat(*,i) = bkj[*,i]^2*lcl
endfor
Bmat = Bmat / (4. * !dpi)

print, '--- Bmat computed.'
   
return, Bmat

end
