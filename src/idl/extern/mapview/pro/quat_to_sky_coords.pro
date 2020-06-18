PRO Quat_to_Sky_Coords, Q, DA_str, Res, Side = Side, $
    ECL=ecl, GAL=gal, CEL=cel, PIXEL=pixel
;+
;  NAME:
;     QUAT_TO_SKY_COORDS
;  PURPOSE:
;     Extract a time series of astronomical coordinates or of A and B
;     side pixel numbers from an array of input quaternions, for a 
;     given DA.
; CALLING SEQUENCE:
;     QUAT_TO_SKY_COORDS, Q, DA_str, [Res, ECL= , GAL=, CEL=, PIXEL = ]
; INPUTS:
;     Q            An array of quaternions of size 4xN.
; OPTIONAL INPUTS:
;     DA_str       Case-insensitive string containing the DA to process:
;               'K1',Ka1',Q1','Q2','V1','V2','W1','W2','W3','W4'.   If not
;               supplied or set to '', then the mean optical axis is used.
;     Res          Integer(1-13) giving the resolution of the returned nested 
;               Healpix  pixel number.    Only required if the output PIXEL 
;               keyword is supplied.
; OPTIONAL INPUT KEYWORD:
;    Side   - String -- either 'A' or 'B' or 'AB' specifying which side of the 
;           spacecraft to compute coordinates for.    Default is 'AB' to compute
;           for both sides.
;                          
; OUTPUT KEYWORDS:
;     If only 1 side is specified, then the following output arrays will be
;     N X 2 rather than N X 4
;
;  ECL          N X 4 array containing ecliptic long & lat for side A,
;               ecliptic long & lat for side B.
;  GAL          N X 4 array containing galactic long & lat for side A,
;               galactic long & lat for side B.
;  CEL          N X 4 array containing RA & Dec in degrees for side A,
;               RA & Dec in degrees for side B.
;  PIXEL        N X 2 array containing pixel number for side A,
;               pixel number for side B.  The Side keyword is ignored.
; EXAMPLE:
;     Find the Side A pointing of the K band in Galactic coordinates for 
;     the first science frame in a TOD file
;
;     IDL> fits_read_tod, 'MAP_tod_20022162357_20022172357.fits', tod
;     IDL> q  = tod[0].quaternions       ;Extract first 4 x 33 quaternion array
;     IDL> quat_to_sky_coords, q[*,1:30], 'K1',gal = gal
;
;     The subscripts [1:30] are used since the first quaternion is for the 
;     preceding science frame, and the last two are for the subsequent science
;     frame
; PROCEDURES USED:
;     COORTRANS, Q2M                     WMAP Library
;     VEC2PIX_NEST                       HealPix Library
;     FTAB_EXT                           IDLAstro Library
; REVISION HISTORY:
;     Converted from G. Hinshaw's Quat_to_Pixel_Series.  R. S. Hill, SSAI,
;                                  21 Jan 2003
;     Make DA_STR parameter optional           W. Landsman SSAI 06 Mar 2003
;-
 if N_params() LT 1 then begin
     print,'Syntax - Quat_to_Sky_Coords, Q, DA_str, [ Res ]'
     print,'         Output Keywords - ECL=, GAL=, CEL=, PIXEL = '
     return
  endif

 if N_elements(side) EQ 1 then begin
      doA = strpos(strupcase(side),'A') GE 0
      doB = strpos(strupcase(side),'B') GE 0 
 endif else begin
      doA = 1 & doB = 1
 endelse     
 doAB = doA and doB

 if N_elements(DA_STR) EQ 0 then da_str = '' 
; Get the direction cosines of the line of sight
 if da_str EQ '' then begin 
    dir_A_LOS = [0.d0, 0.94264149d0, -0.33380686d0]   ; Optical axis
    dir_B_LOS = [0.d0,-0.94264149d0, -0.33380686d0]    
 endif  else begin
    da_str= strtrim(da_str)
    fstring = da_str + 'A,' + da_str + 'B'
    FTAB_EXT, '$MAP_REF/map_los.fits', fstring, Dir_A_LOS, DIR_B_LOS
 endelse
  
; Extract the rotation matrices
Npts = N_ELEMENTS(Q)/4L
Q2M, Q, M
M = TRANSPOSE(M,[2,0,1])

; Get equatorial coordinates in rectangular form
if doA then Dir_A_LOS_cel = TOTAL(M*REBIN(REFORM(Dir_A_LOS,1,1,3),Npts,3,3),3)
if doB then Dir_B_LOS_cel = TOTAL(M*REBIN(REFORM(Dir_B_LOS,1,1,3),Npts,3,3),3)

; Convert to longlat form
IF( ARG_PRESENT(Cel) )THEN BEGIN
  if doA then if doB then CoorTrans, Dir_A_LOS_cel, Cll_a, 'u2ll' $
                     else CoorTrans, Dir_A_LOS_cel, Cel, 'u2ll'
  if doB then if doA then CoorTrans, Dir_B_LOS_cel, Cll_b, 'u2ll' $
                     else CoorTrans, Dir_B_LOS_cel, Cel, 'u2ll'
  if doAB then Cel = [[temporary(Cll_a)],[temporary(Cll_b)]]
ENDIF

; Convert to galactic
IF( ARG_PRESENT(Gal) OR ARG_PRESENT(Pixel) )THEN BEGIN
  if doA then begin 
        CoorTrans, Dir_A_LOS_cel, Dir_A_LOS_gal, 'c2g'
        if doB then CoorTrans, Dir_A_LOS_gal, Cll_A, 'u2ll' $ 
                   else CoorTrans, Dir_A_LOS_gal, Gal, 'u2ll'
  endif
  if doB then begin
        CoorTrans, Dir_B_LOS_cel, Dir_B_LOS_gal, 'c2g'
        if doA then CoorTrans, Dir_B_LOS_gal, Cll_B, 'u2ll' $ 
                   else CoorTrans, Dir_B_LOS_gal, Gal, 'u2ll'
  endif	
  if doAB then Gal = [[temporary(Cll_a)],[temporary(Cll_b)]]
ENDIF

; Convert to ecliptic
IF ARG_PRESENT(Ecl) THEN BEGIN
  if doA then if doB then CoorTrans, Dir_A_LOS_cel, ecl_A, 'c2e', /lonlat else $
                          CoorTrans, Dir_A_LOS_cel, ecl, 'c2e', /lonlat
  if doB then if doA then CoorTrans, Dir_B_LOS_cel, ecl_B, 'c2e', /lonlat else $
                          CoorTrans, Dir_B_LOS_cel, ecl, 'c2e', /lonlat
  if doAB then Ecl = [[temporary(ecl_A)],[temporary(ecl_B) ]  ]
ENDIF

; Extract the pixel numbers
IF ARG_PRESENT(Pixel)THEN BEGIN
  if N_elements(Res) EQ 0 then message, $
      'ERROR - a resolution must be supplied with the PIXEL keyword'
  Nside = 2L^Res
  if doA then if doB then Vec2Pix_Nest, Nside, Dir_A_LOS_gal, Pixel_A $
                     else Vec2Pix_Nest, Nside, Dir_A_LOS_gal, Pixel
  if doB then if doA then Vec2Pix_Nest, Nside, Dir_B_LOS_gal, Pixel_B $
                     else Vec2Pix_Nest, Nside, Dir_B_LOS_gal, Pixel
  if doAB then Pixel = [[temporary(Pixel_a)],[temporary(Pixel_b)]]
 
ENDIF

RETURN
END
