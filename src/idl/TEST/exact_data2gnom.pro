;+
; NAME:
;     EXACT_DATA2GNOM
;
; PURPOSE:
;     This code is a data2gnom version which does not display any comment
;    turns a Healpix into a Gnomonic rectangular map
;
;     exact_data2gnom, data, nside, do_rot, eul_mat, pxsize=, pysize=, rot=,$
;            reso_arcmin=, grid
; INPUT :
;      data --- healpix data
;      nside --- healpix resolution
;      do_rot --- do rotation boolean 
;      eul_mat --- euler matrix to do the rotation
; OUTPUT :
;      grid --- catalog with position and pixel value 
;
; KEYWORDS
;      pxsize --- pixel size on x
;      pysize --- pixel size on y
;      rot --- rotation angle
;      reso_arcmin --- resolution in arcmin
;
;  REVISION HISTORY
;  written by Sandrine Pires, Oct. 2009
;----------------------------------------------------------------------------------------------
pro exact_data2gnom, data, nmap, test, nside, do_rot, eul_mat, PXSIZE=pxsize,PYSIZE=pysize,$
 ROT=rot, RESO_ARCMIN=reso_arcmin, grid=grid, count, wind=wind, proj=proj

if not keyword_set(proj) then proj = 'gnome'
npix = n_elements(data)
print, 'nside=', nside
Nmax = reso_arcmin/60.*!DtoR*pxsize/2.
print, 'Nmax=', Nmax
; -------------------------------------------------------------
; create the rectangular window
; -------------------------------------------------------------
if defined(pxsize) then xsize = pxsize*1L else xsize = 500L
if defined(pysize) then ysize = pysize*1L else ysize = xsize
if defined(reso_arcmin) then resgrid = reso_arcmin/60. else resgrid = 1.5/60.
dx      = resgrid * !DtoR

grid = FLTARR(3, ceil(xsize*0.7)*ceil(ysize*0.7))
; -------------------------------------------------------------
; makes the projection around the chosen contact point
; -------------------------------------------------------------
; position on the planar grid  (1,u,v)
x0 = +1.
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur)  ; & deltax = (xur - xc)
yc = 0.5*(yll+yur)  ; & deltay = (yur - yc)

nband = ysize*1L
npb = xsize*1L * nband
u = -(findgen(xsize) - xc)# replicate(dx,nband); minus sign = astro convention
v = replicate(dx,xsize) # (findgen(nband) - yc)
x = replicate(x0, npb)

vector = [[x],[reform(u,npb)],[reform(v,npb)]] ; non normalised vector
;should use gnome_inv x,y --> theta, phi
if (do_rot) then vector = vector # eul_mat 
vec2pix_nest, nside, vector, id_pix
itest = dblarr(12.*nside*nside)
itest[id_pix] = 1
ind = where(itest eq 1, count)
print, 'npix per patch=', count
if count gt ceil(xsize*0.7)*ceil(ysize*0.7) then count = ceil(xsize*0.7)*ceil(ysize*0.7)
rot = rot*!dtor
for i = 0l, count-1 do begin
  pix2ang_nest, nside, ind[i], theta, phi
	if proj eq 'gnome' then exact_gnome, theta, phi, -rot[1]+!pi/2., rot[0], xx, yy
  if proj eq 'stereo' then exact_stereo, theta, phi, -rot[1]+!pi/2., rot[0], xx, yy
  grid[0, i] = -yy
  grid[1, i] = xx
  if keyword_set(wind) then begin
    ii = xx*0.5/Nmax
    jj = yy*0.5/Nmax
		if wind eq 'hamming' then grid[2, i] = data[ind[i]]*(0.54-0.46*cos(2*!pi*(ii+0.5)))*(0.54-0.46*cos(2*!pi*(jj+0.5)))
	  if wind eq 'hanning' then grid[2, i] = data[ind[i]]*cos(!pi*ii)^2*cos(!pi*jj)^2
	  if wind eq 'cosine' then grid[2, i] = data[ind[i]]*cos(!pi*ii)*cos(!pi*jj)
		if wind eq 'bartlett' then grid[2, i] = data[ind[i]]*(1.-ii-jj)
		if wind eq 'rect' then grid[2, i] = data[ind[i]]*1.
  endif else begin
	  grid[2, i] = data[ind[i]]
  endelse
endfor
test[ind] = test[ind] + 1

end

