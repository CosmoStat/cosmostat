;+
; NAME:
;     ZEA_XY
;
; PURPOSE:
;     Converts from (lon,lat) to Zenithal equal-area projection (x,y) positions. 
;     
; CALLING SEQUENCE:
;     zea_xy, lon, lat, x, y
;
; INPUTS:
;     lon - longitude(s), in degrees,  numeric scalar or vector
;
;     lat - latitude(s), in degrees, numeric scalar or vector
;
; OUTPUTS:  
;     x   - x position on the Zenithal Equal Area projection, scaled such
;                     that -2 < x < 2
;
;     y   - y position on the Zenithal equal area projection, scaled such
;                     that -1 < y < 1
;
; KEYWORDS:
;     None.
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     None.
;
; COMMENTS:
;     Used by higher level map overlay routines.
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 22 Apr 1999
;     Internal computations in Double precision W. Landsman December 2002
;
;-
;======================================================================
PRO zea_xy, Lon, Lat, X, Y

; This routine converts an input set of longitude and latitude values,
; and converts them to Zenithal equal area projected image values in the range:

;    -2 < x < 2    (N pole is at x = -1, S pole is at x = +1)
;    -1 < y < 1
;                 
; projection is such that lon=90 is following -y axis, and lon=0 is at
; projection center (lon=180 at outer edges).
;
; eqn source:  Lambert Azimuthal Equal Area projection
;              Snyder and Voxland, "An album of map projections", p. 226
;              GA110.S575 1989

; The inputs Lon and Lat may be arrays (of equal size).

 sqrt2 = sqrt(2.0d)

; Transform to image coordinates
 dtr = !dpi/180.0d
 north = where(lat ge 0, Nnorth)
 south = where(lat lt 0, Nsouth)

x = dblarr(n_elements(lat)) & y = x

 if (Nnorth GT 0) then begin
   X[north] = sqrt2*sin((45.0d - lat[north]/2.)*dtr)*cos(lon[north]*dtr) 
   X[north] = X[north] -1.
   Y[north] = -1.*sqrt2 * sin((45.0d - lat[north]/2.)*dtr)*sin(lon[north]*dtr)
 endif
 if (Nsouth GT 0) then begin
   X[south] = sqrt2*sin((45.0d + lat[south]/2.)*dtr)*(-1.)*cos(lon[south]*dtr) 
   X[south] = X[south] +1.
   Y[south] = -1.*sqrt2*sin((45.0d + lat[south]/2.)*dtr)*sin(lon[south]*dtr)
 endif

; Return single precision?
 double = (size(lon,/tname) EQ 'DOUBLE') or (size(lat,/tname) EQ 'DOUBLE')
 if not double then begin
      x = float(x) & y = float(y)
 endif

return
END



