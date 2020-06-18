pro starast,ra,dec,x,y,cd, righthanded=right,hdr=hdr     
;+
; NAME:
;       STARAST 
; PURPOSE:
;       Compute astrometric solution using positions of 2 or 3 reference stars
; EXPLANATION:
;       Computes an exact astrometric solution using the positions and 
;       coordinates from 2 or 3 reference stars.   If 2 stars are used, then
;       the X and Y plate scales are assumed to be identical, and the
;       axis are assumed to be orthogonal.   Use of three stars will
;       allow a unique determination of each element of the CD matrix.
;
; CALLING SEQUENCE:
;       starast, ra, dec, x, y, cd, [/Righthanded, HDR = h]
;
; INPUTS:
;       RA - 2 or 3 element vector containing the Right Ascension in DEGREES
;       DEC- 2 or 3 element vector containing the Declination in DEGREES
;       X -  2 or 3 element vector giving the X position of reference stars
;       Y -  2 or 3 element vector giving the Y position of reference stars
; OUTPUTS:
;       CD - CD (Coordinate Description) matrix (DEGREES/PIXEL) determined 
;               from stellar positions and coordinates.
; OPTIONAL INPUT KEYWORD:
;       /RightHanded - If only 2 stars are supplied, then there is an ambiguity
;               in the orientation of the coordinate system.   By default,
;               STARAST assume the astronomical standard left-handed system
;               (R.A. increase to the left).   If /Right is set then a 
;               righthanded coordinate is assumed.  This keyword has no effect
;               if 3 star positions are supplied
; OPTIONAL INPUT-OUTPUT KEYWORD:
;        HDR - If a FITS header string array is supplied, then an astrometry 
;              solution is added to the header using the CD matrix and star 0
;              as the reference pixel (see example).   Equinox 2000 is assumed.
; EXAMPLE:
;        To use STARAST to add astrometry to a FITS header H;
;
;        IDL> starast,ra,dec,x,y,cd       ;Determine CD matrix
;        IDL> crval = [ra[0],dec[0]]      ;Use Star 0 as reference star
;        IDL> crpix = [x[0],y[0]] +1      ;FITS is offset 1 pixel from IDL
;        IDL> putast,H,cd,crpix,crval     ;Add parameters to header
;
;        This is equivalent to the following command:
;        IDL> starast,ra,dec,x,y,hdr=h      
;  
; METHOD:
;       The CD parameters are determined by solving the linear set of equations
;       relating position to local coordinates (l,m)
; REVISION HISTORY:
;       Written, W. Landsman             January 1988
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Added /RightHanded and HDR keywords   W. Landsman   September 2000
;-
 if N_params() LT 4 then begin
        print,'Syntax - STARAST, ra, dec, x, y, cd. [/Right, HDR =h]'
        return                         
 endif

 cdr = !DPI/180.0D 

 nstar = min( [N_elements(ra), N_elements(dec), N_elements(x), N_elements(y)])
 if (nstar NE 2) and (nstar NE 3) then $
        message,'Either 2 or 3 star positions required'

; Convert RA, Dec to Eta, Xi

 ra_rad = ra*cdr       
 dec_rad = dec*cdr
 delx1 = x[1] - x[0] 
 dely1 = y[1] - y[0]     
 delra = ra_rad[1] - ra_rad[0]
 cosdec = cos(dec_rad)  
 sindec = sin(dec_rad) 
 cosra =  cos(delra) 
 sinra = sin(delra)
 denom = sindec[1]*sindec[0] + cosdec[1]*cosdec[0]*cosra
 l1 = cosdec[1]*sinra/denom
 m1 = (sindec[1]*cosdec[0] - cosdec[1]*sindec[0]*cosra)/denom

if nstar EQ 3 then begin

        delx2 = x[2] - x[0] & dely2 = y[2] - y[0]
        delra = ra_rad[2] - ra_rad[0]
        cosra = cos(delra) 
        sinra = sin(delra)
        denom = sindec[2]*sindec[0] + cosdec[2]*cosdec[0]*cosra
        l2 = cosdec[2]*sinra/denom
        m2 = (sindec[2]*cosdec[0] - cosdec[2]*sindec[0]*cosra)/denom
        b = double([l1,m1,l2,m2])
        a = double( [ [delx1, dely1, 0,    0    ], $
                      [0    , 0,     delx1,dely1], $
                      [delx2, dely2, 0,    0    ], $
                      [0    , 0    , delx2,dely2] ] )
        a = transpose(a)        ;Make IDL conform to matrix notation
endif else begin

        b = double( [l1,m1] )
        if keyword_set(right) then  $
              a = double( [ [delx1,dely1], [-dely1,delx1] ] ) else $
              a = double( [ [delx1,-dely1], [dely1,delx1] ] )

endelse

 cd = invert(a)#b        ;Solve linear equations

 if nstar EQ 2 then begin
           if keyword_set(right) then $ 
               cd = [ [cd[0],cd[1]],[-cd[1],cd[0]] ] else $
               cd = [ [cd[0],cd[1]],[cd[1],-cd[0]] ] 
 endif else cd = dblarr(2,2) + cd

 cd = cd/cdr
 
 if N_elements(hdr) GT 0 then begin
        crval = [ra[0],dec[0]]         ;Use Star 0 as reference star
        crpix = [x(0),y(0)] + 1.       ;FITS is offset 1 pixel from IDL
        putast,hdr,cd,crpix,crval,equi=2000.0        ;Add parameters to header
 endif
     
 return
 end
