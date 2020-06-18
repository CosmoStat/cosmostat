;+
; NAME:
;        mrs_median
;
; PURPOSE:
;  Apply a median transform on a healpix spherical image using a disk as structured element 
;
; CALLING:
;     MedIma = mrs_median(Imag, HalfSize=HalfSize)
;
; INPUTS:
;     Imag -- IDL array of a healpix map: Input image to be transformed 
;    
; OUTPUTS:
;     MedIma -- IDL array of a healpix map.
;
; KEYWORD: 
;     WindowSize: float = window size for median operator, default is 5.
;     HalfSize: float = halfsize of the disk in arcmin, by default HalfSize=pixel_size(nside) * WindowSize / 2.
;
; EXAMPLE:
;       Apply an erosion to an image.  The result is stored in Output
;               Output = mrs_median(imag) 
;         
; HISTORY:
;       Written:Francois Xavier Dupe & Jean-Luc Starck, 2011
;--------------------------------------------------------------------------------------------------------

function mrs_median, imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize
rmap=-1
if N_PARAMS() NE 1  then begin 
        print, 'CALLING SEQUENCE: ErdodIma = mrs_median(imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize)'
        goto, DONE
        end

nested=1
npix  = n_elements(imag)
nside = npix2nside(npix)
imag2 = imag

imag2 = reorder(imag,out='ring',in='nest')
rmap = imag2

if not keyword_set( WindowSize) then WindowSize = 5.

if not keyword_set( halfsize) then halfsize=pixel_size(nside) * WindowSize / 2.

; Apply erosion operator
for i = 0ul, npix-1 do begin
     ; Get the xyz coordinate
     pix2vec_ring, nside, i, vec

     ; Get the disk of radius halfsize
     query_disc, nside, vec, float(halfsize)/60.0, lpix, /deg

     rmap[i]  = median( imag2[lpix] )
end     

rmap = reorder(rmap,in='ring',out='nest')

DONE:

return,rmap
end
