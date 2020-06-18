pro fits_reproj_healpix, infile, outfile,coord=coord,projection=projct,size=size 
;+
;   NAME:           
;      FITS_REPROJ_HEALPIX
;   PURPOSE:
;       Convert a HealPix image to a flatmap projection FITS file with WCS info 
;   EXPLANATION:  
;      Writes a FITS file containing a flatmap projection in either a Mollweide
;      or Zenithal Equal Area projection, and write the appropriate World
;      Coordinate System (WCS, see Calabretta & Greisen, 2002, A&A, 395, 1077)
;      information into the FITS header.    This allows coordinates to be 
;      displayed by standard astronomical display software (e.g. DS9, SAOImage).
;
;      Note that the conversion to a flat map projection does not preserve the
;      noise properties of the CMB, and should not be used for a CMB analysis.
;
; CALLING SEQUENCE:
;     FITS_REPROJ_HEALPIX, INFILE, OUTFILE, [COORD={'G','C','E'}, 
;                 PROJECT=['M','Z'],SIZE=[1-5],$
; INPUT ARGUMENT:
;       INFILE   - Scalar string giving name of the FITS file containing the
;                  image in HealPix projection.
;                 'map_k_imap_yr1_v1.fits'
; OUTPUT ARGUMENT:
;       OUTFILE - Scalar string giving the name of the FITS file containing
;                the flat map projection with astrometry.
;
; OPTIONAL INPUT KEYWORDS:
;       User will be prompted for the following values if they are not supplied
;       as keywords:
;
;       PROJECT - Scalar string specifying either a "Mollweide" or a "Zenithal"
;               Equal Area projection.   Only the first letter is needed.
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
;       COORD - One of the three words "Galactic" (default), "Celestial", or
;               "Ecliptic' giving the output projection.   Only the first
;               letter is needed.
;
;   EXAMPLE:
;       Convert the K band all-sky HealPix map to a medium (1024 x 512) size 
;       Mollweide projection stored in a FITS file moll_k_yr1_v1.fits

;       IDL> infile = 'map_k_imap_yr1_v1.fits'
;       IDL> fits_reproj_healpix, infile, 'moll_k_yr1_v1.fits',size=2,proj='Mol'
;
;   NOTES:
;       A Mollweide projection is written in a single FITS file extension.  A
;       Zenithal Equal Area projection is written as two image in two extensions
;        -- the first centered on the North Pole, and the second on the South 
;       Pole.
;   PROCEDURES USED:
;       FITS_READ_MAP GET_HEAL_LUT, GET_HEAL_RES()
;   MODIFICATION HISTORY:
;-
 if N_Params() LT 1 then begin
    print,'Syntax - fits_reproj_healpix, in, outfile'
    print,"              COORD={'G','E','C'}, PROJECT={'M','Z'},SIZE=[1-5] ]"
    return
  endif
; Obtain the look-up table

 FITS_READ_MAP,infile,t,pheader = ph
 if N_elements(coord) EQ 0 then coord = 'G'
    case strupcase(strmid(coord,0,1)) of 
    'G': lcoord = 1
    'C': lcoord = 2
    'E': lcoord = 3
    else: message,'Unrecognized projection type of ' + strtrim(coord,2)
 endcase

;Make sure original HealPix indexing is in Galactic coordinates

 skycoord = sxpar(ph,'SKYCOORD',Count = N_skycoord)
 galactic =  N_skycoord GE 1
 if not galactic then $
    message,'WARNING - SKYCOORD keyword not found in header',/con
 if galactic then galactic = strupcase(strtrim(skycoord,2)) EQ 'GALACTIC'
 if not galactic then $
     message,/CON,'WARNING - Native Galactic coordinates assumed'

 if size(projct,/TNAME) EQ 'STRING' then begin
      case strupcase(strmid(projct,0,1)) of
      'M': project = 1 
      'Z': project = 2
      else: message,'ERROR - Unrecognized projection type'
      endcase
  endif else if N_elements(projct) GT 0 then project  = projct

; Obtain the appropriate look-up table
 GET_HEAL_LUT, Table, coord=lcoord,project=project,size=size

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
 SXADDHIST, ['FITS_REPROJ_HEALPIX:' + systime() + $
    ' Converted to Flat Map Projection', $
     'This FITS file contains a reprojection of WMAP HealPix data into', $
      'a flat map and should not be used for a CMB analysis'],ph
 

      dimen = size(proj,/dimen)
      case lcoord of 
     1: ctype = ['GLON','GLAT']
     2: ctype = ['RA--','DEC-']
     3: ctype = ['ELON','ELAT']
     endcase
     if  project EQ 2 then begin     ;ZEA projection
         imsize = dimen[1]
         im_north = proj[0:imsize-1, *]
         im_south  = proj[imsize:*, *]
	 mkhdr,hn,im_north,/image
         hs  = hn     
 
         imsize2 = imsize/2.
         crpix = replicate(imsize2 + 0.5,2)
	 cdelt = 180*sqrt(2)/(imsize2*!pi)*[-1,1]
	 ctype = ctype + '-ZEA'	 
         crval = [90,90.0]
         MAKE_ASTR,astr,ctype=ctype,delt=cdelt,crpix=crpix,crval=crval, $
	           longpole=0.
         PUTAST,hn, astr
	 crval = [270, -90]
	 cdelt = -cdelt
         MAKE_ASTR,astr,ctype=ctype,delt=cdelt,crpix=crpix,crval=crval, $
	           longpole=180.
         putast,hs,astr,cd_type=1
 	 
	 writefits,outfile,'',ph	   
	 writefits,outfile,im_north,hn,/append
	 writefits,outfile,im_south,hs,/append
 endif else begin                    ;Mollweide projection
         mkhdr,h,proj,/IMAGE
         crpix = dimen/2. + 0.5
	 cdelt = 360.0d*sqrt(2.0d)/!dpi*[-2./dimen[0],1./dimen[1]]
	 ctype = ctype + '-MOL'	 
         crval = [0.,0.]
         make_astr,astr,ctype=ctype,delt=cdelt,crpix=crpix,crval=crval, $
	           longpole=0.
	putast,h,astr, cd_type = 1
	writefits, outfile, '',ph
	writefits,outfile,proj,h,/append
 endelse
  
return
 end 
