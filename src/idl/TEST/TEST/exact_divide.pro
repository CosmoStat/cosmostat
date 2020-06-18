;+
; NAME:
;     exact_divide
;
; PURPOSE:
;     Divide the sphere into small maps
;
; CALLING SEQUENCE:
;     exact_divide, name_in,,map_hd=map_hd, frac=frac, dmax=dmax, dpix=dpix, out=out
;
; INPUT:
;     The Healpix map name_in
;     
; OUTPUT:
;     
; KEYWORDS:
;     map_hd : the header of one of the small maps that will cover the sphere
;     frac : the fraction of overlapping for the small maps
;       NOTE : map_hd and frac are sufficient to define the sphere division
;     out : result structure
;     dpix : output pixel resolution in arcmin
;     dmap : output map size in arcmin
;
; PROCEDURES CALLED:
;     sph_centers
;
; REVISION HISTORY
;     Written by Sandrine Pires Oct. 2009
;-
pro exact_divide, name_in, map_hd=map_hd, frac=frac, dmap=dmap, dpix=dpix, out=out, wind=wind, proj=proj

sph_centers, map_hd=map_hd, frac=frac, lontab, latab

nx_i = SXPAR(map_hd, 'NAXIS1')
nx=long(floor(nx_i/0.7))
vs = size(name_in)
npix = vs[1]
nside = long(npix2nside(npix))
dpix_i = ABS(SXPAR(map_hd, 'CDELT1'))*60d0
dpix = dpix_i*0.7


nmaps=n_elements(lontab)
print, 'nmaps=', nmaps
npix_map=dblarr(nmaps)
TabMap = fltarr(3, nx_i*nx_i, nmaps)
test = fltarr(npix)

for i=0L, nmaps-1 DO BEGIN
  print, 'nmaps_i = ', i
  rotab = [lontab[i],latab[i]]
  rot_ang = rotab
  if defined(rot_ang) then rot_ang = ([rot_ang,0.,0.])(0:2) $
  else rot_ang = [0., 0., 0.]
  eul_mat = euler_matrix_new(rot_ang(0), -rot_ang(1), rot_ang(2), /Deg, /ZYX)
  do_rot = (TOTAL(ABS(rot_ang)) GT 1.e-5)
  exact_data2gnom, name_in, i, test, nside, do_rot, eul_mat, $
    rot=rot_ang, reso_arcmin=dpix, pxsize=nx, pysize=nx, grid=grid, count, wind=wind, proj=proj
  TabMap[*,0:nx_i*nx_i-1,i] = grid[*, 0:nx_i*nx_i-1]
	npix_map[i] = count
endfor
if not keyword_set(dpix) then dpix = 0
if not keyword_set(dmap) then dmap = 0

out= {nside: nside, map_hd: map_hd, nmaps:nmaps, frac:frac, Nx:nx_i, Ny:nx_i, $
Map: TabMap, npix_map : npix_map, lon: lontab, lat: latab, pixelsize:dpix_i, mapsize:dmap}

end
