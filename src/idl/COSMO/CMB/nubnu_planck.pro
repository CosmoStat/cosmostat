PRO NUBNU_PLANCK, LAMBDA, TEMP, SPECTRUM, YSZ= YSZ, DT_T= DT_T
;+
;=-----------------------------------------------------------------------------
; NAME:   
;     NUBNU_PLANCK
; PURPOSE:
;   Computes Planck function as NuBnu in Wm-2sr-1
;	or SZ or dT/T distortion
; CALLING SEQUENCE:
;   NUBNU_PLANCK, LAMBDA, TEMP, SPECTRUM, [YSZ= YSZ], [DT_T= DT_T]
; INPUTS:
;   LAMBDA	: wavelength vector IN Meters
;   TEMP	: Temperature vector in Kelvin
; OPTIONAL INPUT:
; KEYWORD INPUT : 
;   YSZ		: if present compute the Sunyaev-Zeldovich distortion
;   DT_T	: if present compute the deltaT/T type of distortion
; OUTPUTS:
;   SPECTRUM	: Brightness vector per log frequency, NuBnu in Wm-2sr-1
; SIDE EFFECT:
; RESTRICIONS:
;   Either LAMBDA or TEMP can be a vector, not both
; PROCEDURE CALLS:
; HISTORY:
;	08-jun-1993	version 1 FXD IAS
;	16-jan-1994	V1.1	  FXD IAS  Correct the case where T is a vector
;=-----------------------------------------------------------------------------
;-
IF N_PARAMS() EQ 0 THEN BEGIN
	PRINT, 'Call is '
	PRINT, $
  'NUBNU_PLANCK, LAMBDA (meters), TEMP (Kelvin), SPECTRUM (Wm-2sr-1), '
	PRINT, $
  '   YSZ= YSZ, DT_T= DT_T'
	GOTO, CLOSING
ENDIF

HC_K= 0.0143879	;in meter*K
PLANCK_PREFACTOR= 1.191E-16	; 2*h*c2= 1.191E-16   in wm-2sr-1 m4
STUPID = WHERE( TEMP* LAMBDA LE 0., COUNT_STUPID)
X= TEMP* LAMBDA > 1.E-33
X= HC_K/ X
NUBNU_CN= X* 0. ; Define output vector
LAMAUX= X* 0.	; Define intermediate vector of wavelength
LAMAUX(*)= LAMBDA  
IF COUNT_STUPID NE 0 THEN $
	PRINT, 'W - Some nonsense Lambda or temperatures'
X= X< 70.
EXX= EXP( X)
EX1= EXX- 1.
U= WHERE( X LT 1E-5, COUNT)
IF COUNT GT 0 THEN EX1( U)= X( U)
OK_VALUE = WHERE( TEMP* LAMBDA GT 0., COUNT_OK)
IF COUNT_OK GT 0 THEN $
  NUBNU_CN( OK_VALUE) = PLANCK_PREFACTOR/ LAMAUX( OK_VALUE)^4/ EX1( OK_VALUE)

IF KEYWORD_SET( YSZ) THEN $
	SPECTRUM= YSZ* (X * EXX)/ EX1* (X* (EXX+ 1)/EX1 - 4.)* NUBNU_CN ELSE $
IF KEYWORD_SET( DT_T) THEN $
	SPECTRUM= DT_T* X* EXX/ EX1* NUBNU_CN ELSE $
	SPECTRUM= NUBNU_CN	

CLOSING:
RETURN
END
