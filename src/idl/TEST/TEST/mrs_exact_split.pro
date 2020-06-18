;+
; NAME:
;       mrs_exact_split
;
; PURPOSE:
;       Decompose an healpix map (nested format) into patches.
;       The pixels in patches are non regularly distruibuted and
; 			The position and pixel value are stored in a catalog
;
; CALLING:
;
;       Cube = mrs_exact_split(Imag, frac=frac, nx=nx, whanning=whanning)
;
; INPUTS:
;       map -- healpix map (nested format) to be decomposed into patches 
;    
; OUTPUTS:
;       Cube -- IDL structures with the following fields:  
;         Nside : healpix resolution
;         map_hd : the header of one of the small maps that will cover the sphere
;         NMaps : int = number of patchs 
;         Frac: float = overlapping factor between patches
;         hmap: string = header of a small patch
;         Nx,Ny : int = size of each patch
;         Map : fltarr[3,Nx*Ny,NMaps] =  patches (catalog)
;         Npix_map : nb of pixels per patch 
;         Lon : flarr[NMaps] = Longitude of each patch
;         Lat : flarr[NMaps] = Latitude of each patch
;         PixelSize: float = Pixel size in arc minute
;         MapSize: float = Map size in arcmin
;
; KEYWORDS:
;      Frac: overlapping factor  between patches. Default is 0.05
;      Nx: Nb of pixels per patch side. Default is automatically estimated.
;      wind : a window function is applied
;      proj : define the projection from sphere to patch
;
; EXTERNAL CALLS:
;     EXACT_DIVIDE
;
; EXAMPLE:
;     Decompose an image into patches with default options
 ;             Cube = mrs_exact_split(map)
;         
; HISTORY:
;	Written by Sandrine Pires, Oct. 2009
;------------------------------------------------------------------------------
function mrs_exact_split,map,frac=frac,nx=nx,wind=wind, proj=proj
a = systime(1)
if not keyword_set(frac) then frac = 0.05

vs = size(Map)
npix = vs[1]
nside = npix2nside(npix)
; dpix = output pixel resolution in arc minute
dpix = sqrt(4*!pi/!dtor^2/npix)*60d0
dmap = 10.*60d0

if not keyword_set(nx) then nx = fix(dmap/dpix+1) 
print, 'nx=', nx
dmap = dpix*nx

map_hd=generate_header(nx=nx, ny=nx,lon_c=0,lat_c=0,pixdeg=dpix/60d0,angle=0.,coord='G')

exact_divide, map, map_hd=map_hd, frac=frac, dmap=dmap, dpix=dpix, out=out, wind=wind, proj=proj
print, 'time =', systime(1)-a
return, out
end
