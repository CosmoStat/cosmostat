;+
; NAME:
;        mrs_dilation
;
; PURPOSE:
;  Apply a mathematical morphological dilation on a healpix spherical image using a disk as structured element 
;
; CALLING:
;     DilatIma = mrs_dilation(Imag, HalfSize=HalfSize)
;
; INPUTS:
;     Imag -- IDL array of a healpix map: Input image to be dilated 
;    
; OUTPUTS:
;     DilatIma -- IDL array of a healpix map.
;
; KEYWORD: 
;     WindowSize: float = window size for dilation operator, default is 5.
;     HalfSize: float = halfsize of the disk in arcmin, by default HalfSize=pixel_size(nside) * WindowSize / 2.
;
; EXAMPLE:
;       Apply a dilation to an image.  The result is stored in Output
;               Output = mrs_dilation(imag) 
;         
; HISTORY:
;       Written:Francois Xavier Dupe & Jean-Luc Starck, 2011
;--------------------------------------------------------------------------------------------------------

function mrs_dilation, imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize
rmap=-1
if N_PARAMS() NE 1  then begin 
        print, 'CALLING SEQUENCE: ErdodIma = mrs_dilation(imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize)'
        goto, DONE
        end

nested=1
npix  = n_elements(imag)
nside = npix2nside(npix)
imag2 = imag

imag2 = reorder(imag,out='ring',in='nest')
rmap = imag2

if not keyword_set( WindowSize) then WindowSize = 5.

if not keyword_set(halfsizeArcMin) then halfsizeArcMin=pixel_size(nside) * WindowSize / 2.
halfsize = halfsizeArcMin
print, ' halfsize = ', halfsize
; Apply erosion operator
for i = 0ul, npix-1 do begin
     ; Get the xyz coordinate
     pix2vec_ring, nside, i, vec

     ; Get the disk of radius halfsize
     query_disc, nside, vec, float(halfsize)/60.0, lpix, /deg

     rmap[i]  = max( imag2[lpix] )
end     

rmap = reorder(rmap,in='ring',out='nest')

DONE:

return,rmap
end
