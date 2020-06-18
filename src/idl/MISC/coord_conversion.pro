function coord_conversion, map_in, in_coordinate=in_coordinate, out_coordinate=out_coordinate
;+
; NAME:
;       COORD_CONVERSION
;
; PURPOSE:
;       Tool to transform HEALPix maps in Galactic/Ecliptic/Equatorial  coordinates,
;       derived from SISSA/ISAS software (routine coordinate_change).
;
; CALLING SEQUENCE:
;        MapOut = coord_conversion(MapIn, coordinate_in=coordinate_in,coordinate_out=coordinate_out)
;
; INPUTS:
;       map_in = input HEALPix map (nested format)
;
; INPUT KEYWORDS:
;       coordinate_in = coordinates in the output, either 'E' or 'G' or 'Q'
;       coordinate_out = coordinates in the output, either 'E' or 'G' or 'Q' different from coordinate_in
;  
; HISTORY:
;    File create date: April 2013,  Jean-Luc Starck
;-

map_out=-1
coordinate_in='G'
coordinate_out='E'
if keyword_set(in_coordinate) then coordinate_in = in_coordinate
if keyword_set(out_coordinate) then coordinate_out = out_coordinate


if (n_params() lt 1) then begin
    PRINT, 'Error: Wrong number of arguments in coord_conversion'
    print,'Syntax:  MapOut = coord_conversion(MapIn, coordinate_in=coordinate_in,coordinate_out=coordinate_out)'
    goto, DONE
endif

if (coordinate_out eq  coordinate_in) then begin
   print,'Error: Input and output coordinates must be different'
   goto, DONE
endif 


preprocessed = 0
polmap = 0
ordering_in='nested'
nmaps=1
nside = gnside(map_in)
npix=nside2npix(nside)

   
; Coordinate vectors
vec_coordinate_out=fltarr(npix,3L)
ii=lindgen(nside2npix(nside))
pix2vec_nest,nside,ii,vec_coordinate_out
   
  ; Conversion 
;  print,'Conversion'
vec_coordinate_in=rotate_coord(vec_coordinate_out,inco=coordinate_out,outco=coordinate_in)
vec_coordinate_out=0.
vec2pix_nest,nside,vec_coordinate_in,i_coordinate_in
 
map_out=map_in

ii=lindgen(nside2npix(nside))
map_out(ii)=map_in(i_coordinate_in)

DONE:

return, map_out

end

