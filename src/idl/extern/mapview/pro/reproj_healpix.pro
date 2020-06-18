PRO REPROJ_HEALPIX, T, Proj, Subs, MASK=mask, _EXTRA=extra
;+
;   NAME:           
;      REPROJ_HEALPIX
;   PURPOSE:
;       Convert a HealPix image to a flatmap projection 
;   EXPLANATION:  
;      This procedure converts a HealPix image to a flat map projection without
;      doing any intensity scaling, i.e., the pixel values are in their 
;      original units.
;
;   CALLING SEQUENCE:
;       REPROJ_HEALPIX, T, Proj, Subs, [COORD=[1-9], PROJECT=[1,2],SIZE=[1-5],$
;                                        MASK= ] 
;   INPUT ARGUMENT:
;       T   - Vector giving input healpix pixel list, full sky.
;
;   OUTPUT ARGUMENT:
;       Proj   - Output reprojected map, 2-d array
;       Subs   - Subscripts that generate the projection from the healpix 
;                vector.
;
;   INPUT KEYWORDS, 
;       The following keywords are passed to GET_HEAL_LUT via _EXTRA:
;       User will be prompted for the following values if they are not supplied
;       as keywords:
;
;       COORD - Scalar Integer (1-9) giving the coordinate system transformation:
;                                              Display
;                                  Galactic    Celestial    Ecliptic
;                       Native:
;                       Galactic      1            2            3
;                       Celestial     4            5            6
;                       Ecliptic      7            8            9
;
;       PROJECT - Scalar integer (1 or 2) giving the type of projection:
;                       1 - Mollweide
;                       2 - Zenithal Equal Area
;
;       SIZE - Scalar integer giving the size of the output image described by 
;              the table:
;                       1 -- Small    (512 x 256)
;                       2 -- Medium   (1024 x 512)
;                       3 -- Large    (2048 x 1024)
;                       4 -- X large  (4096 x 2048)
;                       5 -- XX large (8192 x 4096, mollweide, native 
;                            coordinates only)
;
;   OUTPUT KEYWORD:
;       Mask  -   Mask outlining the data/nondata ellipse.   Same size as the
;                 the output proj parameter
;   EXAMPLE:
;       Display the K band all-sky map in Galactic coordinates in a medium 
;          (1024 x 512) size Mollweide projection

;       IDL> fits_read_map,'map_k_imap_yr1_v1.fits',t,n
;       IDL> reproj_healpix,t,proj,coord=1,project=1,size=3
;       IDL> tv,bytscl(proj,0,1)
;
;   PROCEDURES USED:
;       GET_HEAL_LUT, GET_HEAL_RES()
;   MODIFICATION HISTORY:
;       Initial delivery: 26 Sept 2001 JP
;       Subs argument added.  RSH
;       Doc header refurbished.  RSH, SSAI, 21 June 2002.
;-
 if N_Params() LT 1 then begin
    print,'Syntax - REPROJ_HEALPIX, T, Proj, [ Subs, '
    print,'              COORD=[1-9], PROJECT=[1,2],SIZE=[1-5], MASK = ]'
    return
  endif
; Obtain the look-up table
 GET_HEAL_LUT, Table, _extra=extra

; Obtain the pixel reduction factor from the size of the input array
 N_Pixels = N_ELEMENTS(T)
 Resolution = GET_HEAL_RES( N_Pixels )
 IF( Resolution GT 10 OR Resolution LT 1 )THEN BEGIN
  PRINT,' N_Pixels, Resolution:', N_Pixels, Resolution
  message,' Map resolution too high or low for the current look-up tables...'
 ENDIF ELSE Factor = 4L^(10-Resolution)

 mask = table gt -1
 subs = (Table > (-1))/Factor
 proj = T[ subs ] * mask
 dd = where(table LT 0, ndd)
 IF ndd GT 0 THEN subs[dd] = -999

return
 end 
