;+
; NAME:
;        mrs_split
;
; PURPOSE:
;       Decompose an healpix map (NESTED format) into a cube of small patches.
;
; CALLING:
;
;     PatchTrans = mrs_split( Imag, frac=frac, nx=nx, SizePatchDegrees=SizePatchDegrees, exrec=exrec, PixelSizeParam=PixelSizeParam )
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be decomposed into patches 
;    
; OUTPUTS:
;     PatchTrans -- IDL structures with the following fields:  
;                     NMaps : long = number of patches 
;                     map_hd: string = header of a small patch
;                     Nx,Ny : long = size of each patch
;                     Lon : flarr[NMaps] = Longitude of each patches if the position (Nx/2.,Ny/2.) is the center of the patches.
;                     Lat : flarr[NMaps] = Latitude of each patches if the position (Nx/2.,Ny/2.) is the center of the patches.
;                     Map : fltarr[Nx,Nx,NMaps] =  patches cube
;                     PixelSize: float = Pixel size in arc minute
;                     MapSize: float = Map size in arcmin
;                     Frac: float = overlapping factor between patches
;					  Nside: long = nside parameter of the input imag
;
; KEYWORDS:
;		Frac : overlapping factor  between patches. Default is 0.05
;		Nx : Number of pixels per patch along both x and y axis. Default is automatically estimated.
;		SizePatchDegrees: Size (in degrees) of each patch. Default is 10 degrees.
;		exrec: if set, the pixel size is smaller in order to have an exact reconstruction using mrs_invsplit
;		PixelSizeParam: if exrec, the pixel size is multiplied by a factor (< 1) in order to smaller pixel for exact reconstruction.
;
; EXTERNAL CALLS:
;     SPH_DIVIDE
;
; EXAMPLE:
;       Decompose an image into patches with default options
 ;             Patch = mrs_split(map)
;         
; HISTORY:
;	Written:  Jean-Luc Starck, June 2007
;--------------------------------------------------------------------------------------------------------


function mrs_getpatchheader, PatchTrans, z=z
if not keyword_set(z) then z=0
Nx = long(PatchTrans.nx)
Ny = long(PatchTrans.ny)
map_hd = GENERATE_HEADER(nx=nx, ny=ny, lon_c=PatchTrans.Lon[z], lat_c=PatchTrans.Lat[z], pixdeg=PatchTrans.PIXELSIZE/60d0, angle=0., coord='G')
return, map_hd
end

;=======================================

function mrs_xyz2lb, PatchTrans, x, y, Frame, map_hd_Frame=map_hd_Frame
nx = N_ELEMENTS(x)
ny = N_ELEMENTS(y)
nz = N_ELEMENTS(Frame)
LB=-1
if N_PARAMS() NE 4 then begin
  print, "ERROR in mrs_xyz2lb: number of paramters must be equal to 4 ..."
  goto, DONE
end
if (Nx NE Nz) and (Nz NE 1)  then begin
  print, "ERROR in mrs_xyz2lb: Nz must be equal to Nx or equal to 1 ..."
  goto, DONE
end

if not keyword_set(map_hd_Frame) then map_hd_Frame=mrs_getpatchheader(PatchTrans, z=Frame)

if (Nx NE Nz) and (Nz EQ 1) then begin
   TabZ = x
   TabZ[*] = Frame
end else TabZ = Frame
   
SPH_LB2XY, l, b, TabZ, x, y, map_hd=map_hd_Frame, frac=PatchTrans.frac, /reverse
np = N_ELEMENTS(l)
LB = fltarr(np, 2)
LB[*,0] = l
LB[*,1] = b  
DONE:

return, LB
end

;=======================================

function mrs_lb2xyz, PatchTrans, l, b, map_hd_Frame=map_hd_Frame

if keyword_set(map_hd_Frame) then map_hd=map_hd_Frame $
else map_hd = mrs_getpatchheader(PatchTrans)

SPH_LB2XY, l, b, Frame, x, y, map_hd=map_hd, frac=PatchTrans.frac 
np = N_ELEMENTS(x)
LB = fltarr(np, 3)
LB[*,0] = x
LB[*,1] = y  
LB[*,2] = Frame  
return, LB
end

;=======================================

function mrs_getpatch, Map, l, b, nx=nx, SizePatchDegrees=SizePatchDegrees

vs = size(Map)
npix = vs[1]
nside = npix2nside(npix)
dpix = SQRT(4*!pi/!dtor^2/npix)*60d0

if not keyword_set(SizePatchDegrees) then  SizePatchDegrees = 10. 

SPH_LOADDATA_HEALPIX, $
 Map, select_in, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
 /online,/NESTED,UNITS=units,COORD=coord,FLIP=flip, $
ROT=rot,QUADCUBE=quadcube,LOG=log,ERROR=error, $
POLARIZATION=polarization, FACTOR=factor, OFFSET=offset

dmap = SizePatchDegrees*60d0 ; 10 degres de cote a la resolution full sky
;print, dmap/dpix
if not keyword_set(nx) then nx = fix(dmap/dpix+1) 
dmap = dpix*nx

  rotab = [l,b]
  rot_ang = rotab
  if defined(rot_ang) then rot_ang = ([rot_ang,0.,0.])(0:2) else rot_ang = [0., 0., 0.]
  eul_mat = euler_matrix_new(rot_ang(0), -rot_ang(1), rot_ang(2), /Deg, /ZYX)
  do_rot = (TOTAL(ABS(rot_ang)) GT 1.e-5)

  SPH_GNOMVIEW, Map, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, rot=rot_ang, reso_arcmin=dpix, pxsize=nx, pysize=nx, /online, /nested, outline=outline, grid=grid

return, grid

end

;=======================================

function mrs_split, map, frac=frac, nx=nx, SizePatchDegrees=SizePatchDegrees, exrec=exrec, PixelSizeParam=PixelSizeParam

; dpix = output pixel resolution in arc minute
; dmap = output pixel size in degrees
; Fraction of overlapping (5 percent)
if not keyword_set(frac) then frac = 0.05
if not keyword_set(PixelSizeParam) then PixelSizeParam = 0.7

vs = size(Map)
npix = vs[1]
nside = npix2nside(npix)
dpix = SQRT(4*!pi/!dtor^2/npix)*60d0
if keyword_set(exrec) then dpix = dpix*PixelSizeParam

if not keyword_set(SizePatchDegrees) then  SizePatchDegrees = 10. 

dmap = SizePatchDegrees*60d0 ; 10 degres de cote a la resolution full sky
;print, dmap/dpix
if not keyword_set(nx) then nx = fix(dmap/dpix+1) 
dmap = dpix*nx

map_hd = GENERATE_HEADER(nx=nx, ny=nx, lon_c=0, lat_c=0, pixdeg=dpix/60d0, angle=0., coord='G')
view=0
verbose=0
SPH_DIVIDE, map, map_hd=map_hd, frac=frac, view=view, verbose=verbose, /online, /nested, out=out, dmap=dmap, dpix=dpix

return, out
end

;=======================================

