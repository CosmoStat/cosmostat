;+
; NAME:
;        mrs_resize_mask_vec
;
; PURPOSE:
;     Compute extended versions of the input mask.
;
; CALLING:
;     rmask = mrs_resize_mask_vec(mask, daddpix, maskEdge=maskEdge)
;
; INPUTS:
;     mask       -- Healpix map = mask composed of zeros and ones.
;     daddpix    -- IDL intarr = number of pixels the mask must be
;                   extended of.
;    
; OUTPUTS:
;     rmask      -- IDL intarr(npix,max(daddpix)) = contains the
;                   resized masks for each integer in daddpix, stored
;                   at the corresponding index.
;
; INPUT KEYWORDS:
;     maskEdge    : Healpix map = map containing ones where the edge
;                   of the input mask is (for faster resizing only).
;
; EXAMPLE:
;      rmask = mrs_mask_edge(mask, [3,6,10])
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;-------------------------------------------------------------------------

function mrs_resize_mask_vec, mask, daddpix, maskEdge=maskEdge

rmask = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: res =  mrs_resize_mask_vec(mask, daddpix, maskEdge=maskEdge)'
   goto, DONE
endif

mask = mask[*]
Nbr = max(daddpix)
npix  = n_elements(mask)
nside = npix2nside(npix)
rmask = repmat(mask,1,Nbr+1) ;dblarr(npix,Nbr)
fnside = double(nside)

radius = pixel_size(nside)*(daddpix)/60*!Dtor 

if not keyword_set(maskEdge) then maskind = where(mask eq 0, maskc) else maskind = where(maskEdge eq 1, maskc)

ind = where(radius ne 0,c)
if c eq 0 then goto, DONE

if maskc gt 0 then begin
   for i=0l, maskc-1 do begin
      ;; Get the xyz coordinate
      pix2vec_nest, nside, maskind[i], vec
            
      for k=0, c-1 do begin
         indk = ind[k]
         ;; Get the disk of radius halfsize
         query_disc, nside, vec, radius[indk], lpix, /Nest, /inclusive
         
         ;; Get neighboring pixels
         ;;neighbours_nest, nside, ind[i], lpix, nneigh
         
         rmask[lpix,daddpix[indk]] = 0
      endfor
   endfor
endif

DONE:
return, rmask

end
