PRO SPH_CENTERS, lontab, latab, map_hd=map_hd, frac=frac
;+
; NAME:
;     SPH_CENTERS
; ONE LINE HELP:
;     Give the center coordinates (lontab,latab) of the small maps covering the sphere
; PURPOSE:
;     Give the center coordinates (lontab,latab) of the small maps covering the sphere
; CALLING SEQUENCE:
;     SPH_CENTERS, lontab, latab, map_hd=map_hd, frac=frac
; INPUT:
;
; OUTPUT:
;     lontab : longitude array of the centers of the small maps 
;     latab : latitude array of the centers of the small maps
; KEYWORDS:
;     map_hd : the header of one of the small maps covering the sphere
;     frac : the fraction of overlapping for the small maps
;       NOTE : map_hd and frac are sufficient to define a sphere division
; EXAMPLE :
;     
; RESTRICTIONS:
;     Not fully tested yet. Report any bugs to Jean-Baptiste Melin
;     jean-baptiste.melin@cea.fr
;
; PROCEDURES CALLED:
;
; REVISION HISTORY
;     Written, Jean-Baptiste Melin, October 11 2006
;-

IF DEFINED(map_hd) EQ 0 THEN BEGIN
  PRINT, 'No small-maps header map_hd. Please define it! Returning...'
  GOTO, closing
ENDIF

IF DEFINED(frac) EQ 0 THEN BEGIN
  PRINT, 'No fraction of overlapping frac. Please define it! Returning...'
  GOTO, closing
ENDIF

nx = SXPAR(map_hd, 'NAXIS1')
d = ABS(SXPAR(map_hd, 'CDELT1'))
theta = FLOAT(nx*d) ; deg
delta = theta-theta*FLOAT(frac)

bmin = 0
bmax = 90
btab = (FINDGEN(CEIL((bmax-bmin)/delta))*delta)+bmin
nb = N_ELEMENTS(btab)

cosb = COS((btab-theta/2)*!dtor)
cosb[0] = 1
lmin = 0
lmax = 360

;tot = 0B ; DEBUG line

FOR i=0, nb-1 DO BEGIN
  deltal = delta/cosb[i]
  ltab = (FINDGEN(CEIL((lmax-lmin)/deltal))*deltal)+lmin
  nl = N_ELEMENTS(ltab)
  ;print, i, nl ; DEBUG line
  ;print, (360-MAX(ltab))/deltal ; DEBUG line
  IF (360-MAX(ltab))/deltal GT 1 THEN PRINT, 'Warning : there is a bug !'
  ;tot = tot+N_ELEMENTS(ltab) ; DEBUG line
  IF i EQ 0 THEN BEGIN
    lontab = ltab
    latab = REPLICATE(btab[0],nl)
  ENDIF ELSE BEGIN
    lontab = [lontab,ltab,ltab]
    latab = [latab,REPLICATE(btab[i],nl),REPLICATE(-btab[i],nl)]
  ENDELSE
ENDFOR

lontab = [lontab,0,0]
latab = [latab,90,-90]

closing:
END