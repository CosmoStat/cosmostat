;+
; NAME:
;     projxy2coord
;
; PURPOSE:
;     Converts the normalized (x,y) position on an image projection into
;     (lon,lat) coordinates.
;     
; CALLING SEQUENCE:
;     projxy2coord,x0,y0,lon0,lat0 [,proj={'M','Z'}]
;
; INPUTS:
;     x0 - float, range 0-1 - Normalized x position on the image.
;     y0 - float, range 0-1 - Normalized y position on the image.
;
; OUTPUTS:  
;     lon0 - float - longitude of position, in degrees.
;     lat0 - float - latitude of position, in degrees.
;
; KEYWORDS:
;     proj - char string - Single character specifying the projection type
;                          of the image.  Either 'M' (Mollweide) or
;                          'Z' (Zenithal Equal area) are allowed.  Case 
;                          insensitive.   Defaults to 'M'.
; 
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     None.
;
; EXAMPLE:
;     See e.g., use within CURSORCOORD.PRO.
;     
; COMMENTS:
;     Longitude and latitude are computed for the native projection of the
;     image (i.e., it assumes (lon,lat) = 0,0 at x0,y0 = .5,.5.
;
;     The routine may bomb if the user attempts to input arrays of x and y,
;     rather than single scalar values for x and y.
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 30 Apr 1999
;
;-
;======================================================================
pro projxy2coord,x0,y0,lon0,lat0,proj=proj
;
; x and y are NORMALIZED coordinated in the range [0,1].
; lon and lat are returned in degrees, assuming a native coord system
; where (0,0) is at the center of the projection.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
; Projection Type
;
if keyword_set(proj) then begin
   proj = strupcase(strtrim(proj,2))
   if ((proj ne 'M') and (proj ne 'Z')) then begin
     print, 'Only Mollweide (M) or ZEA (Z) projections allowed.'
     return
   endif
endif else begin
   proj = 'M'                      ;default = Mollweide
endelse

;
; Convert (x0,y0) to (lon0,lat0) for specified projection
;

case proj of

'M': begin
     ; 
     ; conversion eqn as per make_heal_mollweide_lut.f90
     ; eqn assumes range of x & y is -1 to 1
     ;
     xo = (x0-.5)/.5
     yo = (y0-.5)/.5
     if ( (xo^2 + yo^2) le 1 ) then begin        ;on sky
        psi  = asin(yo)
        lon0 = -(!pi*xo) / cos(psi)
        lat0 = asin ((2.*psi + sin(2.*psi))/!pi)
        lon0 = lon0 * 180./!pi
        lat0 = lat0 * 180./!pi
     endif else begin
        print,'hey, that point is not on the sky projection!'
        lon0 = -999
        lat0 = -999
        return
     endelse
     end

'Z': begin
     ; 
     ; conversion eqn as per inversion of zea_xy.pro
     ; eqn assumes range of x & y is -1 to 1 in EACH hemisphere
     ;
     if (x0 le .5) then begin                     ;northern hemi
        xo = (x0-.25)/.25
        yo = (y0-.5)/.5
        if ( (xo^2 + yo^2) le 1 ) then begin        ;on sky
           rho =  sqrt(xo^2 + yo^2)/sqrt(2.)
           lon0 = atan(-yo,xo)
           lat0 = 2.*(!pi/4. - asin(rho))
           lon0 = lon0 * 180./!pi
           lat0 = lat0 * 180./!pi
        endif else begin
           print,'hey, that point is not on the sky projection!'
           lon0 = -999
           lat0 = -999
           return
        endelse
     endif else begin                             ;southern hemi
        xo = (x0-.75)/.25
        yo = (y0-.5)/.5
        if ( (xo^2 + yo^2) le 1 ) then begin        ;on sky
           rho  = sqrt(xo^2 + yo^2)/sqrt(2.)
           lon0 = atan(-yo,-xo)
           lat0 = 2.*asin(rho) - (!pi/2.)
           lon0 = lon0 * 180./!pi
           lat0 = lat0 * 180./!pi
        endif else begin
           print,'hey, that point is not on the sky projection!'
           lon0 = -999
           lat0 = -999
           return
        endelse
     endelse
     end

endcase

;
;
;
end
