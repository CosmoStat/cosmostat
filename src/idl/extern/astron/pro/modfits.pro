pro MODFITS, filename, data, header, EXTEN_NO = exten_no, ERRMSG = errmsg
;+
; NAME:
;      MODFITS
; PURPOSE:
;      Modify a FITS file by updating the header and/or data array.  
; EXPLANATION:
;      The updated header or array cannot change the size of the FITS file.
;
; CALLING SEQUENCE:
;      MODFITS, Filename, Data, [ Header, EXTEN_NO = ]
;
; INPUTS:
;      FILENAME = Scalar string containing the name of the FITS file  
;                  to be modified.
;
;      DATA - data array to be inserted into the FITS file.   Set DATA = 0
;               to leave the data portion of the FITS file unmodified
;
;      HEADER - FITS header (string array) to be updated in the FITS file.
;
; OPTIONAL INPUT KEYWORDS:
;      EXTEN_NO - scalar integer specifying the FITS extension to modified.  For
;               example, specify EXTEN = 1 or /EXTEN to modify the first 
;               FITS extension. 
; OPTIONAL OUTPUT KEYWORD:
;       ERRMSG - If this keyword is supplied, then any error mesasges will be
;               returned to the user in this parameter rather than depending on
;               on the MESSAGE routine in IDL.   If no errors are encountered
;               then a null string is returned.               
;
   ;
; EXAMPLE:
;     (1) Modify the value of the DATE keyword in the primary header of a 
;             file TEST.FITS.
;
;              IDL> h = headfits('test.fits')      ;Read primary header
;              IDL> sxaddpar,h,'DATE','1994-03-23' ;Modify value of DATE 
;              IDL> modfits,'test.fits',0,h        ;Update header only
;
;       (2) Replace the values of the primary image array in 'test.fits' with 
;               their absolute values
;
;               IDL> im = readfits('test.fits')    ;Read image array
;               IDL> im = abs(im)                  ;Take absolute values
;               IDL> modfits,'test.fits',im        ;Update image array
;
;       (3) Modify the value of the EXTNAME keyword in the first extension
;       
;               IDL> h = headfits('test.fits',/ext)  ;Read first extension hdr
;               IDL> sxaddpar,h,'EXTNAME','newtable' ;Update EXTNAME value
;               IDL> modfits,'test.fits',0,h,/ext    ;Update extension hdr
;
; NOTES:
;       MODFITS performs numerous checks to make sure that the DATA and
;       HEADER are the same size as the data or header currently stored in the 
;       FITS files.    (More precisely, MODFITS makes sure that the FITS file
;       would not be altered by a multiple of 2880 bytes.    Thus, for example,
;       it is possible to add new header lines so long as the total line count 
;       does not exceed the next multiple of 36.)    MODFITS is best
;       used for modifying FITS keyword values or array or table elements.
;       When the size of the data or header have been modified, then a new
;       FITS file should be written with WRITEFITS.
; RESTRICTIONS:
;       (1) Cannot be used to modifiy the data in FITS files with random 
;           groups or variable length binary tables.   (The headers in such
;           files *can* be modified.)
;
; PROCEDURES USED:
;       Functions:   SXPAR(), FXPOSIT(), IS_IEEE_BIG()
;       Procedures:  IEEE_TO_HOST, CHECK_FITS
;
; MODIFICATION HISTORY:
;       Written,    Wayne Landsman          December, 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Fixed possible problem when using WRITEU after READU   October 1997
;       New and old sizes need only be the same within multiple of 2880 bytes
;       Added call to IS_IEEE_BIG()     W. Landsman   May 1999
;       Added ERRMSG output keyword     W. Landsman   May 2000
;       Update tests for incompatible sizes   W. Landsman   December 2000
;-
  On_error,2                    ;Return to user

; Check for filename input

   if N_params() LT 1 then begin                
      print,'Syntax - MODFITS, Filename, Data, [ Header, EXTEN_NO = ]'
      return
   endif

   if not keyword_set( EXTEN_NO ) then exten_no = 0
   if N_params() LT 2 then Header = 0
   nheader = N_elements(Header)
   ndata = N_elements(data)
   printerr =  not arg_present(ERRMSG) 

   if (nheader GT 1) and (ndata GT 1) then begin
        check_fits, data, header, /FITS, ERRMSG = MESSAGE
        if message NE '' then goto, BAD_EXIT
   endif

; Open file and read header information
         
   if exten_no EQ 0 then begin 
         if nheader GT 0 then $
             if strmid( header[0], 0, 8)  NE 'SIMPLE  ' then begin 
                 message = $
                'Input header does not contain required SIMPLE keyword'
                 if printerr then message, 'ERROR ' + message else $
                 errmsg = message
                 return
             endif
         openu, unit, Filename, /GET_LUN, /BLOCK
   endif else begin
         if nheader GT 0 then $
             if strmid( header[0], 0, 8)  NE 'XTENSION' then message, $
             'ERROR - Input header does not contain required XTENSION keyword'
         unit = fxposit( Filename, exten_no) 
         if unit EQ -1 then begin 
                message = 'ERROR opening FITS file'
                goto, BAD_EXIT
         endif
   endelse

   point_lun, -unit, pointlun
   file = fstat(unit)
   nbytesleft = file.size - pointlun 

   hdr = bytarr( 80, 36, /NOZERO )
   if nbytesleft LT 2880 then begin 
           free_lun, unit
           message = 'EOF encountered while reading FITS header'
           if printerr then message,'ERROR - ' + message else errmsg = message
           return
   endif
   readu, unit, hdr
   nbytesleft = nbytesleft - 2880
   oldheader = string( hdr > 32b )

   endline = where( strmid(oldheader,0,8) EQ 'END     ', Nend )
   if Nend GT 0 then oldheader = oldheader[ 0:endline[0] ] 

   while Nend EQ 0 do begin
            if nbytesleft LT 2880 then begin
                free_lun, unit 
                message = 'EOF encountered while reading FITS header'
                if printerr then message,'ERROR - ' + message else $
                           errmsg = message
                return
             endif
   readu, unit, hdr
   nbytesleft = nbytesleft - 2880
   hdr1 = string( hdr > 32b )
   endline = where( strmid(hdr1,0,8) EQ 'END     ', Nend )
   if Nend GT 0 then hdr1 = hdr1[ 0:endline[0] ] 
   oldheader = [ oldheader, hdr1 ]
   endwhile


   if nheader GT 1 then begin
      point_lun, unit, pointlun
      noldheader = N_elements(oldheader)
      if ( (nheader-1)/36) NE ( (Noldheader-1)/36) then begin     ;Updated Dec. 2000
        message = 'FITS header not compatible with existing file '
        message,'Input FITS header contains '+ strtrim(nheader,2) +' lines',/inf
        message,'Current disk FITS header contains ' + strtrim(Noldheader,2) + $
                ' lines',/inf
        goto,BAD_EXIT
        endif

        writeu, unit, byte(header)
   endif 

   if ndata GT 1 then begin
        Naxis = sxpar(oldheader, 'NAXIS')
        bitpix = sxpar( oldheader, 'BITPIX')

        if Naxis GT 0 then begin
            Nax = sxpar( oldheader, 'NAXIS*' )   ;Read NAXES
            nbytes = nax[0]*abs(bitpix/8)
           if naxis GT 1 then for i = 2, naxis do nbytes = nbytes*nax[i-1]
        endif else nbytes = 0

        newbytes = ndata
        case datatype(data) of 
            'INT': newbytes = newbytes*2
            'FLO': newbytes = newbytes*4
            'LON': newbytes = newbytes*4
            'DOU': newbytes = newbytes*8
            'COM': newbytes = newbytes*8
            'UIN': 
            else:
        endcase

        if ((newbytes-1)/2880) NE ( (nbytes-1)/2880) then begin   ;Updated Dec. 2000
        message = 'FITS data not compatible with existing file '
        message,'Input FITS data contains '+ strtrim(newbytes,2) + ' bytes',/inf
        message,'Current disk FITS data contains ' + strtrim(nbytes,2) + $
                ' bytes',/inf
        goto, BAD_EXIT
        endif
        if nheader EQ 0 then begin
                check_fits,data,oldheader,/FITS,ERRMSG = message
                if message NE '' then goto, BAD_EXIT
        endif
        vms = !VERSION.OS EQ "vms"
        Little_endian = not IS_IEEE_BIG()

        junk = fstat(unit)   ;Need this before changing from READU to WRITEU
        if (VMS or Little_endian) then begin
             newdata = data
             host_to_ieee, newdata
             writeu,unit,newdata
        endif else writeu, unit ,data
    endif       

    free_lun, unit
   return 

BAD_EXIT:
    if N_elements(unit) GT 0 then if (unit NE -1) then free_lun, unit
    if printerr then message,'ERROR - ' + message,/CON else errmsg = message
    message,'FITS file ' + filename + ' not modified',/INF
    return
   end 
