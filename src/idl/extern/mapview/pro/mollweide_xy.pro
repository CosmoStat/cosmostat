;+
; NAME:
;     MOLLWEIDE_XY
;
; PURPOSE:
;     Convert from (lon,lat) to Mollweide projection (x,y) positions. 
;     
; CALLING SEQUENCE:
;     mollweide_xy, lon, lat, x, y
;
; INPUTS:
;     lon  - longitudes, in degrees, numeric scalar or vector
;
;     lat - latitudes in degrees, numeric scalar or vector
;
; OUTPUTS:  
;     x   - x position on the Mollweide projection, scaled such that
;                    -2 < x < 2, always at least floating pt.
;
;     y   - y position on the Mollweide projection, scaled such that
;                    -1 < y < 1, always at least floating pt.
;
; KEYWORDS:
;     None.
;
; EXAMPLE:
;     Find the X,Y mapping of the Galactic poles and the Galactic Center.
;     IDL> mollweide_xy,[0,0,0],[-90,0,90],x,y & print,x,y
;          ===>     0.00000     -0.00000      0.00000
;                  -1.00000      0.00000      1.00000
;
;    and we see that the Galacti North and South poles map into [0,1] and 
;   [0,-1], respectively, and the Galactic center maps into [0,0].
; ROUTINES CALLED:
;     mollweide.pro
;
; COMMENTS:
;     Used by higher level map overlay routines.
;     
; MODIFICATION HISTORY:
;     initial version, G. Hinshaw, 29 Apr 1997
;     Update for integer and double precision values of Lon, Lat W.L Dec 2002
;
;-
;======================================================================

PRO Mollweide_XY, Lon, Lat, X, Y

 if N_params() LT 3 then begin
      print,'Syntax - Mollweide_XY, Lon, Lat, X, Y'
      print,'     Input Lon and Lat in *degrees*'
      return
 endif

 dtor = !dpi/180.0d
 
; This routine converts an input set of longitude and latitude values,
; and converts them to Mollweide projected image values in the range:

;    -2 < x < 2
;    -1 < y < 1

; The inputs Lon and Lat may be arrays (of equal size).


; Reset the longitude scale to [-180,+180]
 bad = where( lon GT 180., Nbad)
 IF Nbad GT 0 THEN Lon[bad] = Lon[bad] - 360.0d


; Transform the latitude to the Mollweide angle
Psi = MOLLWEIDE( DtoR*Lat )

; Transform to image coordinates
X = -(Lon/90.0d)*COS(Psi)    ; X in [-2,+2]
Y = SIN(Psi)                 ; Y in [-1,+1]

; Return single precision?
 double = (size(lon,/tname) EQ 'DOUBLE') or (size(lat,/tname) EQ 'DOUBLE')
 if not double then begin
      x = float(x) & y = float(y)
 endif
RETURN
END
