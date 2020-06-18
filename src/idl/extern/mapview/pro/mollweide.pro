;+
; NAME:
;     mollweide
;
; PURPOSE:
;     Computes the Mollweide angle theta for a given array of latitudes delta.

; CALLING SEQUENCE:
;     theta = mollweide(delta)
;
; INPUTS:
;     delta - latitudes, in radians, scalar or vector, float or double
;
; OUTPUTS:  
;     theta = mollweide angle, in radians, same type and number of 
;             elements as delta
;
; KEYWORDS:
;     None.
;
; COMMENTS:
;     Useful primarily as lower level routine called by mollweide_xy.
;     
; MODIFICATION HISTORY:
;     initial version, G. Hinshaw, 29 Apr 1997
;     Use double precision constants   W. Landsman    December 2002
;
;-
;======================================================================

FUNCTION MOLLWEIDE, Delta

; Computes the Mollweide angle theta for a given array of latitudes delta.
; Requires an iteration with a trial value of delta.

Theta = Delta
N = N_ELEMENTS( Delta )

         FOR I = 0,5 DO BEGIN
             FOR J = 0L,(N-1L) DO BEGIN
                 IF( Delta[J] LT -1.57 OR Delta[J] GT 1.57 )THEN BEGIN
                   Theta[J] = Delta[J] 
                 ENDIF ELSE BEGIN 
		  Theta2 = 2.0d*theta[j]
                   Theta[J] = Theta[J] $
                   + (!DPi*SIN(Delta[J]) - Theta2 - SIN(Theta2))$
                     / ( 2.0d + 2*COS(Theta2) )
                 ENDELSE
             ENDFOR
         ENDFOR

if size(delta,/n_dimen) EQ 0 then return, theta[0] else return, theta

END
