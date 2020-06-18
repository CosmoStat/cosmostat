
function WG7refmask, NSIDE = nside, PHYSICAL = physical, region=region, name=name

; This routine to generate a mask of interesting areas for Galactic
; science. Drawn from the Planck Early Papers and to be used as
; reference targets for CompSep work.

; ---------------------------------------------------------------------

; regions
; 1  AG: (l,b)=(164.9,65.5), 5.1x5.1 degree patch
; 2  BOOTES: (l,b)=(57.8,68.6), 12.3x4.0 degree patch
; 3  DRACO: (l,b)=(92.2,38.4), 5.1x5.1 degree patch
; 4  G86: (l,b)=(87.9,59.1), 5.1x5.1 degree patch    
; 5  MC: (l,b)=(57.0,-81.5), 6.0x5.1 degree patch
; 6  N1: (l,b)=(85.4,44.3), 5.1x5.1 degree patch
; 7  NEP: (l,b)=(96.4,30.0), 12.2x12.1 degree patch
; 8  POL: (l,b)=(125.0,27.4), 6.0x10.0 degree patch     
; 9  POLNOR: (l,b)=(125.0,37.4), 6.0x10.0 degree patch    
; 10 SP: (l,b)=(132.4,47.5), 5.1x5.1 degree patch
; 11 SPC: (l,b)=(135.7,29.3), 12.0x8.6 degree patch
; 12 SPIDER: (l,b)=(135.0,39.9), 10.3x10.2 degree patch      
; 13 UMA: (l,b)=(144.2,38.5), 9.0x9.0 degree patch
; 14 UMAEAST: (l,b)=(155.8,37.0), 10.2x6.1 degree patch
; 15 M42: (l,b)=(209.01,-19.38), 3x3 degree patch
; 16 NGC7293/Helix nebula: (l,b)=(36.16,-57.12)
; 17 Perseus molecular cloud: (l,b)=(160.2,-18.5), 5x5 degree patch
; 18 Rho Ophiuchus molecular cloud: (l,b)=(353.0,+17.0), 5x5 degree patch
; 19 Galactic plane region: (l,b)=(37.5,0.0), 20x20 degree patch
; 20 LMC: (l,b)=(279.16,-33.84), 5x5 degree patch
; 21 SMC: (l,b)=(301.61,-44.24), 2x2 degree patch
; 22 Gum1: (l,b)=(245.5,-15), 21x30 degree patch
; 23 Gum2: (l,b)=(268.0,-15), 24x30 degree patch
; 24 Gum3: (l,b)=(245.5,15), 21x30 degree patch
; 25 Gum4: (l,b)=(268.0,15), 24x30 degree patch
; 26 Taurus Molecular Cloud: (l,b)=(170,-17.5), 20x15 degree patch


; map defaults
if n_elements(nside) eq 0 then nside = 2048L
npix = nside2npix(nside)
mask = fltarr(npix)  
pix2ang_ring,nside,lindgen(npix),theta,phi

; coordinates of regions 
lc = [164.9,57.8,92.2,87.9,57.0,85.4,96.4,125.0,125.0,132.4,135.7, $
      135.0,144.2,155.8,209.01,36.16,160.2,353.0,37.5,279.16,301.61, $
      245.5,268.0,245.5,268.0,170.0] * !DtoR
bc = (90.0 - [65.5,68.6,38.4,59.1,-81.5,44.3,30.0,27.4,37.4,47.5,29.3,39.9, $
      38.5,37.0,-19.38,-57.12,-18.5,17.0,0.0,-33.84,-44.24, $
      -15.0,-15.0,15.0,15.0,-17.5]) * !DtoR
dx = [5.1,12.3,5.1,5.1,6.0,5.1,12.2,6.0,6.0,5.1,12.0,10.3,9.0,10.2, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $ 
      21.0,24.0,21.0,24.0,20.0] * !DtoR
dy = [5.1,4.0,5.1,5.1,5.1,5.1,12.1,10.0,10.0,5.1,8.6,10.2,9.0,6.1, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $
      30.0,30.0,30.0,30.0,15.0] * !DtoR

if not keyword_set(region) then region = 0
IRGN = region
; loop over regions
  bu = bc(irgn) - dy(irgn)/2.
  bd = bc(irgn) + dy(irgn)/2.
  if keyword_set(physical) then begin
    ll = lc(irgn) - dx(irgn)/2/sin(theta)
    lr = lc(irgn) + dx(irgn)/2/sin(theta)
  endif else begin
    ll = lc(irgn) - dx(irgn)/2
    lr = lc(irgn) + dx(irgn)/2
  endelse
  in = where( (phi ge ll) and (phi le lr) and (theta ge bu) and (theta le bd) )
  mask(in) = irgn + 1
 
; write mask
mask = reorder(mask, in='ring', out='nest')
  
; Exit
return, mask
end

;----------------------------------------------------


FUNCTION WG7refrgn, region, NSIDE = nside, NESTED = nested, PHYSICAL = physical

; This routine to return the pixels corresponding to a region of
; interest for WG7 Galactic science. Drawn from the Planck Early
; Papers and to be used as reference targets for CompSep work.

; ---------------------------------------------------------------------

; regions
; 1  AG: (l,b)=(164.9,65.5), 5.1x5.1 degree patch
; 2  BOOTES: (l,b)=(57.8,68.6), 12.3x4.0 degree patch
; 3  DRACO: (l,b)=(92.2,38.4), 5.1x5.1 degree patch
; 4  G86: (l,b)=(87.9,59.1), 5.1x5.1 degree patch    
; 5  MC: (l,b)=(57.0,-81.5), 6.0x5.1 degree patch
; 6  N1: (l,b)=(85.4,44.3), 5.1x5.1 degree patch
; 7  NEP: (l,b)=(96.4,30.0), 12.2x12.1 degree patch
; 8  POL: (l,b)=(125.0,27.4), 6.0x10.0 degree patch     
; 9  POLNOR: (l,b)=(125.0,37.4), 6.0x10.0 degree patch    
; 10 SP: (l,b)=(132.4,47.5), 5.1x5.1 degree patch
; 11 SPC: (l,b)=(135.7,29.3), 12.0x8.6 degree patch
; 12 SPIDER: (l,b)=(135.0,39.9), 10.3x10.2 degree patch      
; 13 UMA: (l,b)=(144.2,38.5), 9.0x9.0 degree patch
; 14 UMAEAST: (l,b)=(155.8,37.0), 10.2x6.1 degree patch
; 15 M42: (l,b)=(209.01,-19.38), 3x3 degree patch
; 16 NGC7293/Helix nebula: (l,b)=(36.16,-57.12)
; 17 Perseus molecular cloud: (l,b)=(160.2,-18.5), 5x5 degree patch
; 18 Rho Ophiuchus molecular cloud: (l,b)=(353.0,+17.0), 5x5 degree patch
; 19 Galactic plane region: (l,b)=(37.5,0.0), 20x20 degree patch
; 20 LMC: (l,b)=(279.16,-33.84), 5x5 degree patch
; 21 SMC: (l,b)=(301.61,-44.24), 2x2 degree patch
; 22 Gum1: (l,b)=(245.5,-15), 21x30 degree patch
; 23 Gum2: (l,b)=(268.0,-15), 24x30 degree patch
; 24 Gum3: (l,b)=(245.5,15), 21x30 degree patch
; 25 Gum4: (l,b)=(268.0,15), 24x30 degree patch
; 26 Taurus Molecular Cloud: (l,b)=(170,-17.5), 20x15 degree patch

; map defaults
if n_elements(nside) eq 0 then nside = 2048L
npix = nside2npix(nside)

if keyword_set(nested) then begin
  pix2ang_nest,nside,lindgen(npix),theta,phi
endif else begin
  pix2ang_ring,nside,lindgen(npix),theta,phi
endelse

; coordinates of regions 
lc = [164.9,57.8,92.2,87.9,57.0,85.4,96.4,125.0,125.0,132.4,135.7, $
      135.0,144.2,155.8,209.01,36.16,160.2,353.0,37.5,279.16,301.61, $
      245.5,268.0,245.5,268.0,170.0] * !DtoR
bc = (90.0 - [65.5,68.6,38.4,59.1,-81.5,44.3,30.0,27.4,37.4,47.5,29.3,39.9, $
      38.5,37.0,-19.38,-57.12,-18.5,17.0,0.0,-33.84,-44.24, $
      -15.0,-15.0,15.0,15.0,-17.5]) * !DtoR
dx = [5.1,12.3,5.1,5.1,6.0,5.1,12.2,6.0,6.0,5.1,12.0,10.3,9.0,10.2, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $ 
      21.0,24.0,21.0,24.0,20.0] * !DtoR
dy = [5.1,4.0,5.1,5.1,5.1,5.1,12.1,10.0,10.0,5.1,8.6,10.2,9.0,6.1, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $
      30.0,30.0,30.0,30.0,15.0] * !DtoR

; pick out region
irgn = region - 1
bu = bc(irgn) - dy(irgn)/2.
bd = bc(irgn) + dy(irgn)/2.
if keyword_set(physical) then begin
  ll = lc(irgn) - dx(irgn)/2/sin(theta)
  lr = lc(irgn) + dx(irgn)/2/sin(theta)
endif else begin
  ll = lc(irgn) - dx(irgn)/2
  lr = lc(irgn) + dx(irgn)/2
endelse
pix = where( (phi ge ll) and (phi le lr) and (theta ge bu) and (theta le bd) )

; Exit
return, pix
end




;----------------------------------------------------

function lb_region, reg, TabRegionName=TabRegionName, dx=dx, dy=dy, lc=lc, bc=bc
; ---------------------------------------------------------------------

; regions
; 1  AG: (l,b)=(164.9,65.5), 5.1x5.1 degree patch
; 2  BOOTES: (l,b)=(57.8,68.6), 12.3x4.0 degree patch
; 3  DRACO: (l,b)=(92.2,38.4), 5.1x5.1 degree patch
; 4  G86: (l,b)=(87.9,59.1), 5.1x5.1 degree patch    
; 5  MC: (l,b)=(57.0,-81.5), 6.0x5.1 degree patch
; 6  N1: (l,b)=(85.4,44.3), 5.1x5.1 degree patch
; 7  NEP: (l,b)=(96.4,30.0), 12.2x12.1 degree patch
; 8  POL: (l,b)=(125.0,27.4), 6.0x10.0 degree patch     
; 9  POLNOR: (l,b)=(125.0,37.4), 6.0x10.0 degree patch    
; 10 SP: (l,b)=(132.4,47.5), 5.1x5.1 degree patch
; 11 SPC: (l,b)=(135.7,29.3), 12.0x8.6 degree patch
; 12 SPIDER: (l,b)=(135.0,39.9), 10.3x10.2 degree patch      
; 13 UMA: (l,b)=(144.2,38.5), 9.0x9.0 degree patch
; 14 UMAEAST: (l,b)=(155.8,37.0), 10.2x6.1 degree patch
; 15 M42: (l,b)=(209.01,-19.38), 3x3 degree patch
; 16 NGC7293/Helix nebula: (l,b)=(36.16,-57.12)
; 17 Perseus molecular cloud: (l,b)=(160.2,-18.5), 5x5 degree patch
; 18 Rho Ophiuchus molecular cloud: (l,b)=(353.0,+17.0), 5x5 degree patch
; 19 Galactic plane region: (l,b)=(37.5,0.0), 20x20 degree patch
; 20 LMC: (l,b)=(279.16,-33.84), 5x5 degree patch
; 21 SMC: (l,b)=(301.61,-44.24), 2x2 degree patch
; 22 Gum1: (l,b)=(245.5,-15), 21x30 degree patch
; 23 Gum2: (l,b)=(268.0,-15), 24x30 degree patch
; 24 Gum3: (l,b)=(245.5,15), 21x30 degree patch
; 25 Gum4: (l,b)=(268.0,15), 24x30 degree patch
; 26 Taurus Molecular Cloud: (l,b)=(170,-17.5), 20x15 degree patch
; 27 SEP: (l,b) = (-83.6, -29.8), 20x15 degree patch
; 28 GAC: (l,b) = (180., 0), 20x15 degree patch
; 29 Cluster A2143: (l,b)=(45, 60), 24x30 degree patch
; 30 Coma: (l,b)=(0,90), 20x15 degree patch
; 31 PsA2163, (l,b) =(45, −60), 20x15 degree patch
; 32 Virgo (l,b) = (135, −45).: , 20x15 degree patch
; 33 Galactic plane region 2: (l,b)=(80,0.0), 20x20 degree patch
; 34 Galactic plane region 3: (l,b)=(340,20), 20x20 degree patch

TabRegionName =['AG', 'BOOTES', 'DRACO', 'G86', 'MC', 'N1', 'NEP', 'POL', 'POLNOR', 'SP', $
         'SPC', 'SPIDER', 'UMA', 'UMAEAST',  'M42', 'NGC7293/Helix nebula', 'Perseus molecular cloud', 'Rho Ophiuchus molecular cloud', 'Galactic plane region', $
	 'LMC', 'SMC', 'Gum1', 'Gum2', 'Gum3', 'Gum4', 'Taurus Molecular Cloud', 'SEP', 'GAC', $
	 'A2143', 'Coma', 'PsA2163', 'Virgo', 'Galactic plane region 2', 'Galactic plane region 3']

; coordinates of regions 
lc = [164.9,57.8,92.2,87.9,57.0,85.4,96.4,125.0,125.0,132.4,135.7, $
      135.0,144.2,155.8,209.01,36.16,160.2,353.0,37.5,279.16,301.61, $
      245.5,268.0,245.5,268.0,170.0, -83.6, 180., 45., 70, 45, 135, 80, 340]  
bc = [65.5,68.6,38.4,59.1,-81.5,44.3,30.0,27.4,37.4,47.5,29.3,39.9, $
      38.5,37.0,-19.38,-57.12,-18.5,17.0,0.0,-33.84,-44.24, $
      -15.0,-15.0,15.0,15.0,-17.5, -29.8, 0, 60., 88., -60., -45., 0, 20] 
dx = [5.1,12.3,5.1,5.1,6.0,5.1,12.2,6.0,6.0,5.1,12.0,10.3,9.0,10.2, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $ 
      21.0,24.0,21.0,24.0,20.0,20.0, 20., 10, 10, 10, 10, 5, 2]  
dy = [5.1,4.0,5.1,5.1,5.1,5.1,12.1,10.0,10.0,5.1,8.6,10.2,9.0,6.1, $
      3.0,3.0,5.0,5.0,20.0,5.0,2.0, $
      30.0,30.0,30.0,30.0,15.0, 15.0, 15., 10, 10, 10, 10, 5, 2]  

InfoReg = [ lc[reg], bc[reg]]
return, InfoReg
end

;===========================================================

pro WG7view, Map, region=region, PHYSICAL = physical, res=res, Mask=Mask, patch=patch, title=title, min=min, max=max, amin=amin, Ima=Ima, png=png, Nl=Nl, Nc=Nc, Filter=Filter
NP =   N_PARAMS() 
 
Nested=1
Nside = gnside(Map)
if not keyword_set(res) then res=3
if not keyword_set(region) then region=0

InfoReg = lb_region(region, TabRegionName=TabRegionName, dx=dx, dy=dy, lc=lc, bc=bc)
 
TabName = TabRegionName
; This routine to return the pixels corresponding to a region of
; interest for WG7 Galactic science. Drawn from the Planck Early
; Papers and to be used as reference targets for CompSep work.

if NP LT 1 then begin
  vs = size(TabName)
  N = vs[1]
  for i=0,N-1 do print, i, ' Reg = ', TabName[i], ",     (l,b) = ", "(", STRC(lc[i]), ",", STRC(bc[i]), ")"
  goto, DONE
end	 
; map defaults
if n_elements(nside) eq 0 then nside = 2048L
npix = nside2npix(nside)

if keyword_set(nested) then begin
  pix2ang_nest,nside,lindgen(npix),theta,phi
endif else begin
  pix2ang_ring,nside,lindgen(npix),theta,phi
endelse

pxsize=long(dx[region]/float(res)*60.)
if pxsize LT 512 then pxsize = 512
if keyword_set(Nc) then pxsize=Nc
pysize=long(dy[region]/float(res)*60.)
Name = TabName[region]
if pysize LT 512 then pysize = 512
if keyword_set(Nl) then pysize=Nl
print, Name, '(l,b)=', lc[region], bc[region] 

Tit = Name 
if keyword_set(title) then tit = tit + ':' + title
if keyword_set(amin) then min=-amin
if keyword_set(amin) then max=amin

print, 'PP ' , pxsize, pysize

if keyword_set(Filter) then begin
   mrs_convol, Map, Filter, MapF
   gnomview, /nested, MapF, rot=[ lc[region],  bc[region]], RESO_ARCMIN=res,  pxsize=pxsize, pysize=pysize, title=tit, min=min, max=max,  map_out=ima, png=png
end else $
gnomview, /nested, Map, rot=[ lc[region],  bc[region]], RESO_ARCMIN=res,  pxsize=pxsize, pysize=pysize, title=tit, min=min, max=max,  map_out=ima, png=png

; Mask = WG7refmask(nside=nside, region=region)
; gnomview, /nested, Mask, rot=[ lc[region],  bc[region]], RESO_ARCMIN=res, pxsize=pxsize, pysize=pysize

; Exit
DONE:

end

;----------------------------------------------------

function mrs_get_area, Map, lb_coord,  Nx=Nx, Ny=Ny, Res=Res, Region=Region,  min=min, max=max, WINDOW=WINDOW, NScale=NScale, Filter=Filter
if not keyword_set(lb_coord) then lb_coord = [135., -45.]  ; VIRGO
if not keyword_set(Nx) then Nx = 512
if not keyword_set(Ny) then Ny = 512
if not keyword_set(Res) then Res = 0.8
if keyword_set(Region) then lb_coord = lb_region(Region, TabRegionName=TabRegionName)
if keyword_set(Region) then tit = TabRegionName [Region]

if keyword_set(Filter) then begin
   mrs_convol, Map, Filter, MapF
   gnomview, /nested, MapF, rot=[ lb_coord[0],  lb_coord[1]], RESO_ARCMIN=res,  pxsize=Nx, pysize=Ny,  min=min, max=max, map_out=map_out, title=tit,WINDOW=WINDOW
end else $
gnomview, /nested, Map, rot=[ lb_coord[0],  lb_coord[1]], RESO_ARCMIN=res,  pxsize=Nx, pysize=Ny,  min=min, max=max, map_out=map_out, title=tit,WINDOW=WINDOW

if keyword_set(NScale) then begin
   W = star2d(map_out, nscale=nscale)
   for i=0, NScale-1 do begin
      window, i, xsize=512, ysize=512, tit=tit +': Scale '+STRC(i+1)
      tvscl, W[*,*,i] 
   end
end

return, map_out
end

;----------------------------------------------------
; tvreg, res=1.5, reg=33, i, filter=f, amin=85

pro tvreg, Map, region=region, PHYSICAL = physical, res=res, Mask=Mask, patch=patch, title=title, min=min, max=max, amin=amin, Ima=Ima, png=png, Nl=Nl, Nc=Nc, Filter=Filter, F10=F10
if  N_PARAMS() EQ 0 then begin
    WG7view
    print, "Calling Syntax: tvreg, Map, region=region, PHYSICAL = physical, res=res, Mask=Mask, patch=patch, title=title, min=min, max=max, amin=amin "
end else begin
   if keyword_set(F10) then begin
      b5 = getbeam(fwhm=5)
      b10 = getbeam(fwhm=10)
      Filter = b5 - b10
      Res=1.5
   end
  WG7view, Map, region=region, PHYSICAL = physical, res=res, Mask=Mask, patch=patch,min=min, max=max, amin=amin , title=title, Ima=Ima, png=png, Nl=Nl, Nc=Nc, Filter=Filter
end
end

;----------------------------------------------------
  
 
