;+
; NAME:
;        mrs_mask_edge
;
; PURPOSE:
;     Estimates the edge of the input mask.
;
; CALLING:
;     edge = mrs_mask_edge(mask)
;
; INPUTS:
;     mask       -- Healpix map = mask composed of zeros and ones.
;    
; OUTPUTS:
;     edge       -- Healpix map = map containing ones where the edge
;                   of the input mask is.
; EXAMPLE:
;      edge = mrs_mask_edge(mask)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------

function mrs_mask_edge, mask

edge = -1
if N_PARAMS() LT 1  then begin 
   print, 'CALLING SEQUENCE: res =  mrs_mask_edge(mask)'
   goto, DONE
endif

npix  = n_elements(mask)
nside = npix2nside(npix)

edge = dblarr(npix)

ind = where(mask eq 0, c)

if c gt 0 then begin
   for i=0L, c-1 do begin
    ; Get neighboring pixels
      neighbours_nest, nside, ind[i], lpix
      
      indk = where( mask(lpix) ne 0, ck)
      if ck gt 0 then edge(ind[i]) = 1
   endfor
endif

DONE:
return, edge

end
