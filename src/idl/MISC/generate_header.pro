FUNCTION GENERATE_HEADER, nx=nx, ny=ny, lon_c=lon_c, lat_c=lat_c, pixdeg=pixdeg, angle=angle, $
	coord=coord
;+
; NAME:
;     GENERATE_HEADER
; ONE LINE HELP:
;     Generate a header for a small rectangular map, tangent projection
; PURPOSE:
;     Generate a header for a small rectangular map, tangent projection, for use in
;     local map representation (e.g. IRAS type)
; CALLING SEQUENCE:
;     h = GENERATE_HEADER(nx=nx, ny=ny, lon_c=lon_c, lat_c=lat_c, pixdeg=pixdeg, angle=angle, coord=coord)
; INPUT:
;     
; OUTPUT:
;     a header for standard fits rectangular maps, tangential projection
; KEYWORDS:
;     nx=nx, ny=ny, lon_c=lon_c, lat_c=lat_c, pixdeg=pixdeg, angle=angle, coord=coord
;     keywords to replace interactive filling of header fields
; EXAMPLES:
;     
; ALGORITHM:
;     read in values interactively
; SIDE EFFECTS:
;     
; RESTRICTIONS:
;     Report any bugs to Jacques Delabrouille
;     j.delabrouille@cdf.in2p3.fr
; PROCEDURES CALLED:
;     SXPAR
; REVISION HISTORY
;     Written, Jacques Delabrouille July 2000
;     Bug in text for parameters inputs fixed, JD, 24-07-2000
;     Now possible to enter formula for pixel size, JD, 24-07-2000
;     Commented field for units, JD, 26/07/2000
;     Keywords added, JD, april 2001
;-

inputstring = ''
inputlong = 0L
inputdouble = 0D

IF NOT KEYWORD_SET(nx) THEN BEGIN
	PRINT, 'ENTER NUMBER OF PIXELS IN X DIRECTION'
	READ, 'NUMBER OF PIXELS IN X DIRECTION : ',inputlong
	mx=inputlong
ENDIF ELSE mx=nx
IF NOT KEYWORD_SET(ny) THEN BEGIN
	PRINT, 'ENTER NUMBER OF PIXELS IN Y DIRECTION'
	READ, 'NUMBER OF PIXELS IN Y DIRECTION : ',inputlong
	my=inputlong
ENDIF ELSE my=ny
MKHDR, head, 4, [mx,my]
SXADDPAR, head,'CRPIX1',mx/2.
SXADDPAR, head,'CRPIX2',my/2.

IF DEFINED(lon_c) EQ 0 THEN BEGIN
	PRINT, 'ENTER LONGITUDE OF MAP CENTER (DEGREES)'
	READ, 'LONGITUDE OF MAP CENTER : ',inputdouble
	lon_c = inputdouble
ENDIF 
SXADDPAR, head,'CRVAL1',lon_c

IF DEFINED(lat_c) EQ 0  THEN BEGIN
	PRINT, 'ENTER LATITUDE OF MAP CENTER (DEGREES)'
	READ, 'LATITUDE OF MAP CENTER : ',inputdouble
	lat_c = inputdouble
ENDIF 
SXADDPAR, head,'CRVAL2',lat_c

IF NOT KEYWORD_SET(pixdeg) THEN BEGIN
	res = 0
	WHILE res NE 1 DO BEGIN
	   PRINT, 'ENTER PIXEL SIZE (DEGREES)'
	   READ, 'PIXEL SIZE : ',inputstring
	   res = EXECUTE('pixdeg = '+inputstring)
	ENDWHILE
ENDIF
SXADDPAR, head,'CDELT1',-pixdeg
SXADDPAR, head,'CDELT2',pixdeg

IF DEFINED(angle) EQ 0 THEN BEGIN
PRINT, 'ENTER NORTH POLE POLAR ANGLE ON MAP (DEGREES)'
READ, 'NORTH POLE POLAR ANGLE : ',inputdouble
angle=inputdouble
ENDIF
SXADDPAR, head,'LONGPOLE',-90.-angle

IF KEYWORD_SET(coord) EQ 0 THEN BEGIN
	PRINT, 'ENTER COORDINATES (AND EPOCH)'
	PRINT, 'G FOR GALACTIC, E FOR ECLIPTIC, Q1987.5 FOR EQUATORIAL AT EPOCH 1987.5'
	READ, 'COORDINATES : ', inputstring
ENDIF ELSE inputstring = STRLOWCASE(coord)
CASE STRLOWCASE(STRMID(inputstring,0,1)) OF
'g': BEGIN
	cootype1='GLON-TAN'
	cootype2='GLAT-TAN'
END
'e': BEGIN
	cootype1='ELON-TAN'
	cootype2='ELAT-TAN'
END
'q': BEGIN
	cootype1='RA---TAN'
	cootype2='DEC--TAN'
	IF STRLEN(inputstring) GT 1 THEN epoch = FLOAT(STRMID(inputstring,1)) ELSE epoch = 2000.
	SXADDPAR, head,'EPOCH',epoch
END
ENDCASE
SXADDPAR, head,'CTYPE1',cootype1
SXADDPAR, head,'CTYPE2',cootype2
head = head(WHERE(STRLEN(head) NE 0))
RETURN, head

END
