PRO SPH_LB2XY, l, b, num, x, y, map_hd=map_hd, frac=frac, reverse=reverse
;+
; NAME:
;     SPH_LB2XY
; ONE LINE HELP:
;     Convert spherical coordinates (l,b) to small maps coordinates (num,x,y)
; PURPOSE:
;     Convert spherical coordinates (l,b) to small maps coordinates (num,x,y) and the opposite if reverse is set
; CALLING SEQUENCE:
;     SPH_LB2XY, l, b, num, x, y, map_hd=map_hd, frac=frac, reverse=reverse
; INPUT:
;     l, b : longitude and latitude in the spherical coordinates
;     OR
;     num : number of the small map
;     x, y : pixel of the small map
; OUTPUT:
;     num : number of the small map
;     x, y : pixel of the small map
;     OR
;     l, b : longitude and latitude in the spherical coordinates
; KEYWORDS:
;     map_hd : the header of one of the small maps covering the sphere
;     frac : the fraction of overlapping for the small maps
;       NOTE : map_hd and frac are sufficient to define a sphere division
;     reverse : set this keyword to compute (l,b) from (num,x,y)
; EXAMPLE :
;     
; RESTRICTIONS:
;     Not fully tested yet. Report any bugs to Jean-Baptiste Melin
;     jean-baptiste.melin@cea.fr
;
; PROCEDURES CALLED:
;     SPH_CENTERS, WCSXY2SPH
;
; REVISION HISTORY
;     Written, Jean-Baptiste Melin, October 11 2006
;-

IF KEYWORD_SET(reverse) THEN BEGIN
   nxx = N_ELEMENTS(x)
   nyy = N_ELEMENTS(y)
   nnum = N_ELEMENTS(num)
   IF (nxx NE nyy) OR (nxx NE nnum) OR (nyy NE nnum) THEN BEGIN
     PRINT, 'x, y and num must have the same dimension. Returning...'
     GOTO, closing
   ENDIF
   SPH_CENTERS, lontab, latab, map_hd=map_hd, frac=frac
   l = FLTARR(nxx)
   b = FLTARR(nxx)
   reso = FLOAT(ABS(SXPAR(map_hd,'CDELT1')))
   nx = SXPAR(map_hd,'NAXIS1')
   xx = (x-(nx-1)/2.)*reso
   yy = (y-(nx-1)/2.)*reso
   yy = -yy
   FOR i=0L, nxx-1 DO BEGIN
     center = [lontab[num[i]],latab[num[i]]]
     longpole= 0
     ctype12 = ['GLON-TAN','GLAT-TAN']
     WCSXY2SPH, xx[i], yy[i], lo, la, ctype=ctype12, crval=center, longpole=longpole
     l[i] = lo
     b[i] = la
   ENDFOR
ENDIF ELSE BEGIN
   nl = N_ELEMENTS(l)
   nb = N_ELEMENTS(b)
   IF nl NE nb THEN BEGIN
     PRINT, 'l and b must have the same dimension. Returning...'
     GOTO, closing 
   ENDIF
   SPH_CENTERS, lontab, latab, map_hd=map_hd, frac=frac
   num = INTARR(nl)
   x = FLTARR(nl)
   y = FLTARR(nl)
   lontab_rad = lontab*!dtor
   l_rad = l*!dtor
   latab_rad = latab*!dtor
   b_rad = b*!dtor
   nmaps = 3 ; maximum number of maps investigated
   FOR i=0L, nl-1 DO BEGIN
     ; distance on the sphere
     theta=2*ASIN(SQRT((SIN((latab_rad-b_rad[i])/2))^2+COS(latab_rad)*COS(b_rad[i])*(SIN((lontab_rad-l_rad[i])/2))^2))
     so = SORT(theta)
     is_in = 0
     index = 0
     WHILE (is_in EQ 0) AND (index LE nmaps-1) DO BEGIN
       ;plot, lontab, latab, psym=4 ; DEBUG line
       ;oplot, [l[i]], [b[i]], psym=5 ; DEBUG line
       whmin = so[index]
       ;oplot, lontab[so[0:nmaps-1]], latab[so[0:nmaps-1]], psym=2 ; DEBUG line
       ;oplot, [lontab[whmin]], [latab[whmin]], psym=1, color=0 ; DEBUG line
       center = [lontab[whmin],latab[whmin]]
       longpole = 0
       ctype12 = ['GLON-TAN','GLAT-TAN']  
       WCSSPH2XY, l[i], b[i], xx, yy, ctype=ctype12, crval=center, longpole=longpole
       reso = FLOAT(ABS(SXPAR(map_hd,'CDELT1')))
       nx = SXPAR(map_hd,'NAXIS1')
       yy = -yy
       xx = FLOAT(xx/reso+(nx-1)/2.)
       yy = FLOAT(yy/reso+(nx-1)/2.)
       ;print, whmin, xx, yy ; DEBUG line
       IF (xx GE -0.5) AND (xx LE nx-0.5) AND (yy GE -0.5) AND (yy LE nx-0.5) THEN BEGIN
         is_in=1
         num[i] = whmin
         x[i] = xx
         y[i] = yy
       ENDIF ELSE BEGIN
         IF index EQ 0 THEN BEGIN
           num[i] = whmin
           x[i] = xx
           y[i] = yy
         ENDIF
         index=index+1
       ENDELSE
     ENDWHILE
     ;IF index GE nmaps THEN BEGIN
     ;  print, i, ' is out of small maps'
     ;ENDIF
   ENDFOR
ENDELSE

closing:
END