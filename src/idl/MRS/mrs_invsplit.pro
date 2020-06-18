function mrs_invsplit, PatchTrans
;+
; NAME:
;        mrs_split
;
; PURPOSE:
;       Reconstruct a healpix map (NESTED format) from a cube of small patches obtained by the routine mrs_split.
;
; CALLING:
;
;     Map = mrs_invsplit(PatchTrans)
;
; INPUT:
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
; OUTPUT:
;     Imag -- IDL array of healpix map: output image reconstructed from the patches 
;    
;
; KEYWORDS:
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Reconstruct an image from patches  
 ;             map = mrs_invsplit(Patch)
;         
; HISTORY:
;	Written:  Jean-Luc Starck, August 2007
;--------------------------------------------------------------------------------------------------------

;Weight = PatchTrans.map
;Weight[*]=0
Nx = long(PatchTrans.nx)
Ny = long(PatchTrans.ny)
Nz = long(PatchTrans.nmaps)
x = lindgen(Nx)
y = lindgen(Ny)
ImaX = fltarr(Nx,Ny)
Imay = fltarr(Nx,Ny)
Wima = fltarr(Nx,Ny)
for i=0,Ny-1 do ImaX[*,i] = x
for i=0,Nx-1 do ImaY[i,*] = y

D = sqrt( float(((Nx)/2.)^2+((Ny)/2.)^2.))
Dx = sqrt( float(((Nx)/2.)^2))
Dy = sqrt( float(((Ny)/2.)^2))

for i=0,Ny-1 do $
for j=0,Nx-1 do Wima [j,i] = cos( float(ABS(i+0.5-Ny/2.))/Dy*!PI/2.) * cos( float(ABS(j+0.5 - Nx/2.)) / Dx *!PI/2.) 
 
ImaX = reform(ImaX, Nx*Ny)
ImaY = reform(ImaY, Nx*Ny)
Wima = reform(Wima , Nx*Ny)

Map = fltarr(12L*PatchTrans.nside*PatchTrans.nside)
Weight = Map
; map_hd = mrs_getpatchheader(PatchTrans)
ok=1
if OK EQ 1 then begin
for iz = 0L, Nz-1 do begin
  z=iz
  PatchIma = PatchTrans.map[*,*,iz] * Wima
  map_hd = mrs_getpatchheader(PatchTrans, z=z)
  lb =  mrs_xyz2lb(PatchTrans, ImaX, ImaY, iz, map_hd=map_hd)
  ; if z EQ 20 then print, max( lb[*,0]), max( lb[*,1])
  lb2ang, lb[*,0] ,lb[*,1], theta,phi
  ; if z EQ 20 then print, max(theta), max(phi)
  ang2pix_nest, PatchTrans.nside, theta, phi, ipix
  Map[ipix] = Map[ipix] + PatchIma   ;  PatchTrans.map[*,*,iz] ; *Wima[ipix]
  Weight[ipix] = Weight[ipix] + Wima   ; 1. ; Wima[ipix]
end
 ind = where( Weight NE 0, c)
 if c GT 0 then Map[ind] = Map[ind] / Weight[ind]
 ind = where( Weight EQ 0, c)
 if c GT 0 then BEGIN
   print, "Unseen pixels = ", c
   pix2lb,  PatchTrans.nside, ind, l, b
   map_hd = mrs_getpatchheader(PatchTrans)
   XYZ = mrs_lb2xyz(PatchTrans, l, b, map_hd_Frame=map_hd_Frame)
   XYZ[*,0] = long(XYZ[*,0] + 0.5) < (nx-1)
   XYZ[*,1] = long(XYZ[*,1] + 0.5) < (ny-1)
   XYZ[*,2] = long(XYZ[*,2] + 0.5) < (nz-1)
   for p=0l,N_ELEMENTS(Ind) -1 do  Map[Ind[p]] = PatchTrans.map[XYZ[p,0], XYZ[p,1], XYZ[p,2]]
 END
end else begin
  ipix = lindgen(N_ELEMENTS(Map))
  pix2lb,  PatchTrans.nside, ipix, l, b
  map_hd = mrs_getpatchheader(PatchTrans)
  XYZ = mrs_lb2xyz(PatchTrans, l, b, map_hd_Frame=map_hd_Frame)
   
  XYZ[*,0] = long(XYZ[*,0] + 0.5) < (nx-1)
  XYZ[*,1] = long(XYZ[*,1] + 0.5) < (ny-1)
  XYZ[*,2] = long(XYZ[*,2] + 0.5) < (nz-1)
  for p=0l,N_ELEMENTS(Map) -1 do  Map[p] = PatchTrans.map[XYZ[p,0], XYZ[p,1], XYZ[p,2]]
end


DONE:

return, Map
end

pro tt
;c1 = fltarr(2048L^2*12)
;c1[*]=1.
c1 = rims('Spheremap_128.fits')
tc1 = mrs_split(c1, /exrec, SizePatchDegrees=50,PixelSizeParam=0.68)     

tc1 = mrs_split(d, sizepatchdegrees=30, frac=0.5, /ex, PixelSizeParam=0.6)
rc1 =  mrs_invsplit(tc1)
tvs, rc1
info, rc1
info, c1 - rc1
tvs, c1 - rc1
; mrs_find_obj, d,  t, NbrObj=NbrObj, level=level,  imaEll=ie , /plot 
end

