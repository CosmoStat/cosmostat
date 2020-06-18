pro getrot, hdr, rot, cdelt, DEBUG = debug      ;Get rotation and scale factor from header
;+
; NAME:
;    GETROT
; PURPOSE:
;     Return the rotation and plate scale of an image from its FITS header
; EXPLANATION:
;     Derive the counterclockwise rotation angle, and the X and Y scale
;     factors of an image, from a FITS image header.   Input parameter 
;     may be either a FITS image header or an astrometry structure (as 
;     obtained by EXTAST.PRO)
;
; CALLING SEQUENCE:
;     GETROT, Hdr, [ Rot, CDelt, /DEBUG  ]   
;             or 
;     GETROT, Astr, Rot, CDelt, /DEBUG ]       
;
; INPUT PARAMETERS:
;     HDR - FITS Image header (string array).  Program will extract the 
;             astrometry structure
;              or
;     ASTR -  ASTROMETRY structure, of the type returned by EXTAST.
;             See the documentation for EXTAST.PRO for details.
;
; OPTIONAL OUTPUT PARAMETERS:
;       ROT - Scalar giving the counterclockwise rotation of NORTH in DEGREES 
;               from the +Y axis of the image.
;       CDELT- 2 element vector giving the scale factors in DEGREES/PIXEL in 
;               the X and Y directions.  Values correspond to the FITS 
;               parameters CDELT1 and CDELT2 
;
;       If no output variables are supplied (or /DEBUG is set), then GETROT 
;       will display the rotation and plate scale at the terminal.
;
; OPTIONAL INPUT KEYWORD
;       /DEBUG - if DEBUG is set, GETROT will print the rotation for both the 
;       X and Y axis when these values are unequal.  If DEBUG is set to 2, 
;       then the output parameter ROT will contain both X and Y rotations.
;
; PROCEDURE:
;       If the FITS header already contains CDELT (and CD or CROTA) keyword,
;       (as suggested by the proposed Greisen & Calabretta FITS standard) 
;       then this is used for the scale factor.
;       
;       If the header contains CD keywords but no CDELT keywords (as in IRAF
;       headers) then the scale factor is derived from the CD matrix. 
;       
; REVISION HISTORY:
;       Written W. Landsman STX January 1987 
;       Convert to IDL V2. M. Greason, STX, May 1990
;       Option to return both rotations added.  J. D. Offenberg, STX, Aug 1991
;       Use new astrometry structure   W. Landsman  Mar 1994
;       Recognize a GSSS header        W. Landsman  June 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Correct rotation determination with unequal CDELT values WL October 1998
;       Consistent conversion between CROTA and CD matrix  WL  October 2000
;-
 On_error,2

 if N_params() EQ 0 then begin
        print,'Syntax: GETROT, Hdr, [ Rot, CDelt ]'
        print,'    OR: GETROT, Astr, [ Rot, CDelt ]'
        return
 endif

 if not keyword_set(DEBUG) then debug = 0
 radeg = 180.0/!DPI
 sz = size(hdr)                ;Did user supply a FITS header or a CD matrix?

 if ((sz[0] eq 1) and (sz[2] eq 7)) then begin                ;FITS header?

        extast,hdr,astr             ;Extract astrometry from header,
        if strmid(astr.ctype[0],5,3) EQ 'GSS' then begin
                hdr1 = hdr
                gsss_stdast, hdr1
                extast, hdr1, astr
        endif
        cd = astr.cd/RADEG      ;then extract CD matrix from astrometry.
        if N_elements(cd) NE 4 then $
            message,'ERROR - Header is missing astrometry keywords CD or CDELT'

endif else $
 if ((sz[sz[0]+1] eq 8) AND (sz[sz[0]+2] eq 1)) then begin   ;ASTROMETRY structure
      astr = hdr
      cd = astr.cd/RADEG
 endif else message, $
        'ERROR - First parameter must be an image header or astrometry structure'

 if astr.cdelt[0] NE 1.0 then begin
        cdelt = astr.cdelt
        cd[0,0] = cd[0,0]*cdelt[0] & cd[0,1] =   cd[0,1]*cdelt[0]
        cd[1,1] = cd[1,1]*cdelt[1] & cd[1,0] =   cd[1,0]*cdelt[1]
 endif else  cd = astr.cd/RADEG                                        

 det = cd[0,0]*cd[1,1] - cd[0,1]*cd[1,0]
 if det LT 0 then sgn = -1 else sgn = 1
 if det GT 0 then $
   message,'WARNING - Astrometry is for a right-handed coordinate system',/INF
 cdelt = fltarr(2)
 if (cd[1,0] eq 0) or (cd[0,1] eq 0) then begin ;Unrotated coordinates?
   rot = 0.
   rot2 = 0.
   cdelt[0] = cd[0,0]
   cdelt[1] = cd[1,1]
 endif else begin
   rot  = atan(  sgn*cd[0,1],  sgn*cd[0,0] ) 
   rot2 = atan( -cd[1,0],  cd[1,1] )
  
 endelse

 if (abs(rot) NE  abs(rot2)) then begin
      if keyword_set(debug) then $
        print,'X axis rotation:',rot*!RADEG, ' Y axis rotation:',rot2*!RADEG
        if debug eq 2 then rot = [rot,rot2] else $
      if (rot - rot2)*!RADEG LT 2 then rot = (rot + rot2)/2.
 endif

  cdelt[0] =   cd[0,0]/cos(rot)
  cdelt[1] =   cd[1,1]/cos(rot)

 rot = float(rot*RADEG)
 cdelt = float(cdelt*RADEG)

DONE:
 if N_params() EQ 1 or keyword_set(DEBUG) then begin
        if debug LT 2 then print,'Rotation (counterclockwise)',rot,' degrees'
        print,'Sampling interval X axis',cdelt[0]*3600.,' arc seconds/pixel'
        print,'                  Y axis',cdelt[1]*3600.,' arc seconds/pixel'
 endif

 return
 end
