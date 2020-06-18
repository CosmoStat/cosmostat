pro imdbase,hdr,catalogue,list,XPOS=xpos,YPOS=ypos, $
                XRANGE=xrange,YRANGE=yrange, SUBLIST = sublist
;+
; NAME:
;     IMDBASE
; PURPOSE:
;     Find the sources in an IDL database that are located on a given image.
;
; CALLING SEQUENCE:
;    imdbase, hdr, catalogue, [list, XPOS= ,YPOS=, XRANGE= ,YRANGE= , 
;                       SUBLIST = ]  
;
; INPUTS:
;    hdr - FITS image header containing astrometry, and the NAXIS1,
;               NAXIS2 keywords giving the image size
;    catalogue - string giving name of catalogue in database
;              Database must contain the (preferably indexed)
;              fields RA (in hours) and DEC.   Type DBHELP for a 
;              list of the names of available catalogues.
;
; OPTIONAL OUTPUT PARAMETER:
;    LIST - A long vector containing the entry numbers of sources found
;           within the image.   This vector can then be used with other
;           database procedures, e.g. to print specified fields (DBPRINT)
;           or subselect with further criteria (DBFIND)
;
; OPTIONAL OUTPUT KEYWORD PARAMETER:
;     XPOS - REAL*4 vector giving X positions of catalogue sources found 
;            within the image
;     YPOS - REAL*4 vector giving Y positions of catalogue sources found 
;            within the image
;
; OPTIONAL INPUT KEYWORD PARAMETERS
;     XRANGE - 2 element vector giving the X range of the image to consider.
;              The default is to search for catalogue sources within the entire
;             image
;     YRANGE - 2 element vector giving the Y range of the image to consider.
;     SUBLIST - vector giving entries in the database to consider in the
;               search.  If not supplied, or set equal to -1, then all entries
;               are considered.
;
; NOTES:
;     If an output list vector is not supplied, then the found objects are
;     diplayed at the terminal.
;
; EXAMPLE:
;      Find all existing IUE SWP observations within the field of the STSDAS
;      disk file FUV0435FC.  Subselect those taken with the SWP camera
;
;      SXHREAD,'fuv0435fc',H             ;Read header from disk
;      IMDBASE,H,'IUE',list              ;Find IUE obs. within image 
;      LIST2 = DBFIND('CAM_NO=3',LIST)   ;Subselect on SWP images
;
; SIDE EFFECTS:
;      The IDL database is left open upon exiting IMDBASE.
; NOTES:
;      IMDBASE checks the description of the RA item in the database for the
;      string '2000'.    If found, the database RA and Dec are assumed to be 
;      in equinox J2000.   Otherwise they are assumed to be in equinox B1950 
;
; SYSTEM VARIABLES:                                                
;      The non-standard system variable !TEXTOUT is required for use with the
;      database procedures.   
;
; PROCEDURES USED:
;      AD2XY, DBEXT, DB_ITEM, DB_ITEM_INFO(), DBOPEN, DBFIND(), EXTAST, 
;      GET_EQUINOX(), GSSSADXY, GSSSXYAD, HPRECESS, SXPAR(), XY2AD 
; REVISION HISTORY:
;      Written W. Landsman            September, 1988
;      Added SUBLIST keyword          September, 1991
;      Updated to use ASTROMETRY structures    J.D. Offenberg, HSTX, Jan 1993
;      Conversion for precession fixed.   R.Hill, HSTX, 22-Apr-93
;      Check RA description for equinox   W. Landsman  Aug 96
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Call HPRECESS if header equinox does not match DB  W. Landsman Oct. 1998
;-
 On_error,2                                  ;Return to caller

 if N_params() LT 2 then begin               ;Sufficient parameters?
     print,'Syntax - imdbase, hdr, catalogue, [ list, SUBLIST = '
     print,'              XPOS = , YPOS = , XRANGE =, YRANGE =  ]'
     print,'Type DBHELP for available catalogues'
     return
 endif              

; Check if catalogue has preselected output fields

 catname = strupcase(strtrim(catalogue,2)) 

 dbopen,catalogue,unavail=unavail       ;Was database found?
 if unavail EQ 1 then message,'Database ' + catalogue + ' is unavailable'
                
 db_item,'ra',itnum
 descrip = db_item_info('description',itnum[0])
 if strpos(descrip,'2000') GE 0 then cat_year = 2000. else cat_year = 1950.

; Get X and Y of 4 corners of the image

 if N_elements(xrange) NE 2 then begin    
      xmin = 0          & xmax =  sxpar(hdr,'NAXIS1') - 1
 ENDIF ELSE BEGIN
      xmin = xrange[0]  & xmax = xrange[1]
 ENDELSE

 if N_elements(yrange) NE 2 then BEGIN
        ymin=0          & ymax = sxpar(hdr,'NAXIS2') - 1
 ENDIF ELSE BEGIN
     ymin = yrange[0]   & ymax = yrange[1]
 ENDELSE

 x = [xmin,xmax,xmax,xmin]
 y = [ymin,ymin,ymax,ymax]

; Make sure header has astrometry and convert X,Y to Ra, Dec

 extast, hdr, ASTR, noparams
 if noparams LT 0 then message,'Image header does not contain astrometry'

; Compare equinox of image with that of database and precess if necessary

 im_year = GET_EQUINOX(hdr,code)
 if ( code EQ -1 ) then begin 
        message,/inf,'EQUINOX keyword not found in header, assumed to be J2000'
        im_year = 2000.    ;Assume image in 2000 Equinox as default
 endif
 if ( im_year NE cat_year ) then begin      ;Need to precess header?
      hdr1 = hdr
      hprecess,hdr1,cat_year
      extast,hdr1, ASTR, noparams
 endif 

 proj = strmid(astr.ctype[0],5,3)           ;Astrometric projection type

 case proj of
         'GSS': gsssxyad, astr, x, y, ra,dec
         else: xy2ad, x, y, ASTR, ra, dec
 endcase

 ra = ra/15.                            ;Convert from degrees to hours
 ramin = min(ra)  & ramax = max(ra)     ;Get max and min RA values
 decmin = min(dec) & decmax = max(dec)  ;Get max and min Dec values
 if (ramax - ramin) GT 12 then begin    ;Does the RA cross 24 hours?
     newmax = ramin
     ramin = ramax
     ramax = 24.
     redo = 1
endif else redo = 0
if not keyword_set(SUBLIST) then sublist = -1


 search = strtrim(ramin,2)  + ' < ra < ' + strtrim(ramax,2) + ', ' + $
         strtrim(decmin,2) + ' < dec < ' + strtrim(decmax,2)
 print,'IMDBASE: Now searching ',catname,' catalogue - be patient'
 print,search

 list = dbfind(search,sublist,/SILENT, Count = nstar)        ;Search for stars in field
 if redo then begin
     search = '0 < ra < ' + strtrim(newmax,2) + ', ' + $
               strtrim(decmin,2) + '< dec <' + strtrim(decmax,2)
     print,search
stop
     newlist = dbfind(search,sublist,/SILENT, Count = count)
     if count GT 0 then list = [list,newlist]
     nstar = nstar + count
 endif
 print,''

 if nstar GT 0 then begin                     ;Any stars found?
   dbext,list,'ra,dec',ra,dec                ;Extract RA,DEC of stars found
   ra = ra*15. 

   case proj of
           'GSS': gsssadxy, astr,ra,dec,x,y
            else: ad2xy,ra,dec,astr,x,y
   endcase 

   good = where( (x GT xmin) and ( x LT xmax ) $     ;Select stars within field
             and (y GT ymin) and ( y LT ymax), ngood)
   if ngood GT 0 then begin 
       list = list[good]
       xpos = x[good] & ypos = y[good]
      message,strtrim(ngood,2)+' '+ catname +' sources found within image',/INF
       if ( N_params() LT 3 ) then dbprint,list,textout=1   ;List stars found
   endif else GOTO,NO_MATCH
 endif else GOTO,NO_MATCH
return

NO_MATCH: message,'No '+ catname + ' sources found within supplied image',/CON
return

end
