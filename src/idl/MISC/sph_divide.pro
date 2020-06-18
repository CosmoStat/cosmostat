PRO SPH_DIVIDE, map_hd=map_hd, frac=frac, dir_in=dir_in, name_in, dir_out=dir_out, name_out=name_out, view=view, verbose=verbose, online=online, nested=nested, out=out, lontab=lontab, latab=latab, dmap=dmap, dpix=dpix
;+
; NAME:
;     SPH_DIVIDE
; ONE LINE HELP:
;     Divide the sphere into small maps
; PURPOSE:
;     Divide the sphere into small maps
; CALLING SEQUENCE:
;     SPH_DIVIDE, map_hd=map_hd, frac=frac, dir_in=dir_in, name_in=name_in, dir_out=dir_out, name_out=name_out, view=view, verbose=verbose
; INPUT:
;     The Healpix map whose name is dir_in/name_in
; OUTPUT:
;     The small maps written as dir_out/name_out
; KEYWORDS:
;     map_hd : the header of one of the small maps that will cover the sphere
;     frac : the fraction of overlapping for the small maps
;       NOTE : map_hd and frac are sufficient to define the sphere division
;     dir_in : the directory where the full sky Healpix map is located
;     name_in : the name of the full sky Healpix map
;     dir_out : the directory where you want to write the small maps
;     name_out : the name of the small maps (the number of the small map will be added automatically)
;                the format of the output is : dir_out+name_out+'_*.fits'
;     view : display the sphere division
;     verbose : gives comments
; EXAMPLE :
;     
; RESTRICTIONS:
;     Not fully tested yet. Report any bugs to Jean-Baptiste Melin
;     jean-baptiste.melin@cea.fr
;
; PROCEDURES CALLED:
;     SPH_CENTERS, GNOMVIEW
;
; REVISION HISTORY
;     Written, Jean-Baptiste Melin, October 11 2006
;     Position of the GNOMVIEW window fixed to the origin, March 16 2007
;-

SPH_CENTERS, lontab, latab, map_hd=map_hd, frac=frac
 
nx = SXPAR(map_hd, 'NAXIS1')
dpix = ABS(SXPAR(map_hd, 'CDELT1'))*60d0

nmaps = N_ELEMENTS(lontab) 

TabMap = fltarr(nx,nx,nmaps)

IF KEYWORD_SET(view) THEN SIMPLEVISU, BYTARR(nx,nx)

; SPH_LOADDATA_HEALPIX extracted from SPH_GNOMVIEW
  if not keyword_set(online) then begin
SPH_LOADDATA_HEALPIX, $
dir_in+name_in, select_in,$
data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
SAVE=save,ONLINE=online,NESTED=nested,UNITS=units,COORD=coord,FLIP=flip, $
ROT=rot,QUADCUBE=quadcube,LOG=log,ERROR=error, $
POLARIZATION=polarization, FACTOR=factor, OFFSET=offset
 end else begin
 SPH_LOADDATA_HEALPIX, $
 name_in, select_in,$
data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
SAVE=save,ONLINE=online,NESTED=nested,UNITS=units,COORD=coord,FLIP=flip, $
ROT=rot,QUADCUBE=quadcube,LOG=log,ERROR=error, $
POLARIZATION=polarization, FACTOR=factor, OFFSET=offset
 end
if error NE 0 then return
nside=0

IF KEYWORD_SET(verbose) THEN  print, ' Number of patches = ', nmaps
FOR i=0L, nmaps-1 DO BEGIN
    IF KEYWORD_SET(verbose) THEN BEGIN
    IF i/10 EQ i/10. THEN BEGIN
      PRINT, 'Cutting map ', STRTRIM(STRING(i),1), ' over ', STRTRIM(STRING(nmaps-1),1), '...'
    ENDIF
  ENDIF
  rotab = [lontab[i],latab[i]]
  
  rot_ang = rotab
  ; From SPH_LOADDATA_HEALPIX
  ; ---- extra rotation ----
  if defined(rot_ang) then rot_ang = ([rot_ang,0.,0.])(0:2) else rot_ang = [0., 0., 0.]
  ; eul_mat = euler_matrix(-rot_ang(0), -rot_ang(1), -rot_ang(2), /Deg, /ZYX)
  eul_mat = euler_matrix_new(rot_ang(0), -rot_ang(1), rot_ang(2), /Deg, /ZYX)
  do_rot = (TOTAL(ABS(rot_ang)) GT 1.e-5)

  if not keyword_set(online) then begin
     ; SPH_GNOMVIEW, dir_in+name_in, rot=rotab, reso_arcmin=dpix, pxsize=nx, pysize=nx, fits=dir_out+name_out+'_'+STRTRIM(STRING(i),1)+'.fits', online=online, nested=nested, outline=outline $
     SPH_GNOMVIEW, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, rot=rot_ang, reso_arcmin=dpix, pxsize=nx, pysize=nx, fits=dir_out+name_out+'_'+STRTRIM(STRING(i),1)+'.fits'
  end else begin
     vs = size(name_in)
     npix = vs[1]
     nside = npix2nside(npix)
      SPH_GNOMVIEW, name_in, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, rot=rot_ang, reso_arcmin=dpix, pxsize=nx, pysize=nx, online=online, nested=nested, outline=outline, grid=grid
      TabMap[*,*,i] = grid
 end
   IF KEYWORD_SET(view) THEN BEGIN
    map = READFITS(dir_out+name_out+'_'+STRTRIM(STRING(i),1)+'.fits', mhd, /SILENT)
    SIMPLEVISU, map, /nonew
  ENDIF
ENDFOR

if not keyword_set(dpix) then dpix = 0
if not keyword_set(dmap) then dmap = 0

if keyword_set(online) then  out= { nside: nside, map_hd: map_hd, nmaps:nmaps, frac:frac, Nx:Nx, Ny:Nx, Map: TabMap, lon: lontab, lat: latab, pixelsize:dpix, mapsize:dmap}


IF KEYWORD_SET(view) THEN WDELETE

END
