pro fitsdir ,directory, TEXTOUT = textout,NoTelescope = NoTelescope, $
   EXTEN = exten
;+
; NAME:
;     FITSDIR 
; PURPOSE:
;     Provide a brief description of the headers of FITS disk files.  
; EXPLANATION:
;     The values of the FITS keywords NAXISi, OBS-DATE (or TDATEOBS or DATE),
;     TELESCOPE (or TELNAME, OBSERVAT, or INSTRUME), OBJECT (or TARGNAME), 
;     EXPTIME (or INTEG or EXPOSURE) are displayed.    All of these are 
;     commonly used FITS keywords and all except EXPTIME are officially 
;     reserved FITS keywords.   Keyword names in parentheses are searched if 
;     the primary keyword is not found.
;
;     In IDL V5.3 or later, FITSDIR will also recognize gzip compressed files
;     (must have a .gz extension).
; CALLING SEQUENCE:
;     FITSDIR , [ directory, TEXTOUT =, /NoTelescope, EXTEN= ] 
;
; OPTIONAL INPUT PARAMETERS:
;     DIRECTORY - Scalar string giving file name, disk or directory to be 
;             searched.   Wildcard file names are allowed.    Examples of 
;             valid VMS or Unix names include '*.fit' or 'tape*'.    An 
;             example of a valid VMS name is  'UIT$USER2:[JONES]*.FIT' while
;                a valid Unix string is 'iraf/*.fits'.
;
;             If not given, FITSDIR searches *.fits files in the default 
;             directory.
;
; OPTIONAL KEYWORD INPUT PARAMETER
;     /NOTELESCOPE - If this keyword is set and non-zero then the value of 
;               the (usually less important) TELESCOPE keyword is not 
;               displayed, and more space is available to display the other 
;               keyword values.    The following table shows the default 
;               output formats with and without the /NOTELESCOPE keyword.
;
;       File name    NAXISi    OBS-DATE    TELESCOPE    OBJECT    EXPTIME
;          A18        A11        A10          A10        A20        F7.1
;          A20        A12        A10                     A29        F7.1
;
;                                                                    
;      TEXTOUT - Controls output device as described in TEXTOPEN procedure
;               textout=1       TERMINAL using /more option
;               textout=2       TERMINAL without /more option
;               textout=3       <program>.prt
;               textout=4       laser.tmp
;               textout=5       user must open file
;               textout=7       Append to existing <program>.prt file
;               textout = filename (default extension of .prt)
;
;      EXTEN - Integer>0 specifying which extension header to first read to 
;              find the FITS keywords.   By default, FITSDIR only looks at the 
;              primary header.   If /Exten is set, then FITSDIR first looks at
;              the first extension (if present) for the keywords, and then
;              looks at the primary header.   The only disadvantage to including
;              the EXTEN keyword is that it slows the processing slightly.
;
; OUTPUT PARAMETERS:
;       None.
;
; RESTRICTIONS:
;       (1) Field values may be truncated if their length exceeds the default
;           format.
;
;       (2) Does not print an information about the columns in FITS tables. 
;           Use FTAB_HELP to get this information.    Use FITS_HELP to get
;           information about every extension in a FITS file.    
;
;       (3) Users may wish to modify the program to display other FITS 
;               keywords of particular interest to them
; EXAMPLES:  
;       IDL> fitsdir          ;Print info on all '*.fits' files in current 
;                               directory.     
;       IDL> fitsdir ,'*.fit'   ;Lists all '*.fit' files in current directory 
;       IDL> fitsdir ,'*.fit.gz',/ext   ;Print info on all *.fit.gz files in 
;                               ;current directory.   Files are asummed to be 
;                               ;compressed.   The first extension header
;                               ;will be read along with the primary header.
;
;       Write info on all *.fits files in the Unix directory /usr2/smith, to a 
;       file 'smith.txt' and don't display the value of the TELESCOPE keyword
;
;       IDL> fitsdir ,'/usr2/smith/*.fits',t='smith.txt', /NoTel 
;
; PROCEDURE:
;       FINDFILE is used to find the specified FITS files.   The header of
;       each file is read, and rejected if the file is not FITS.    Each header 
;       searched for the parameters NAXISi, TELESCOP (or TELNAME or INSTRUME), 
;       OBJECT (or TARGNAME), DATE-OBS (or TDATEOBS or DATE) and 
;       EXPTIME (or TEXPTIME or EXPOSURE or INTEG).  
;
; SYSTEM VARIABLES:
;       The non-standard system variables !TEXTOUT and !TEXTUNIT must be 
;       defined before calling FITS_INFO.   
;
;       DEFSYSV,'!TEXTOUT',1
;       DEFSYSV,'!TEXTUNIT',0
;
;       One way to define these is to call the procedure ASTROLIB.   
;       See TEXTOPEN.PRO for more info
; PROCEDURES USED:
;       FDECOMP, FXMOVE, MRD_HREAD, REMCHAR,  SPEC_DIR(), 
;       TEXTOPEN, TEXTCLOSE, ZPARCHECK
; MODIFICATION HISTORY:
;       Written, W. Landsman,  HSTX    February, 1993
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Search alternate keyword names    W.Landsman    October 1998
;       Avoid integer truncation for NAXISi >32767  W. Landsman  July 2000
;       Don't leave open unit    W. Landsman  July 2000 
;       Added EXTEN keyword, work with compressed files, additional alternate
;       keywords W. Landsman     December 2000
;-
 On_error,2

 if N_params() GT 0 then $
     zparcheck, 'FITSDIR ', directory, 1, 7, 0, 'Directory Name' $
 else directory = '*.fits'
 if N_elements(exten) EQ 0 then exten = 0 

 fdecomp, directory, disk, dir, filename, ext
 if filename EQ '' then begin 
      directory = disk + dir + '*.fits'
      filename = '*'
      ext = 'fits'
 endif else if !VERSION.OS_FAMILY EQ 'unix' then begin
        if (strpos(filename,'*') LT 0) and (ext EQ '') then begin  
        directory = disk + dir + filename + '/*.fits'
        filename = '*'
        ext = 'fits'
        endif
 endif

 if keyword_set(NoTelescope) then begin 
          namelen = 20
          objectlen = 29 
 endif else begin 
          namelen = 18
          objectlen = 20
 endelse

 direct = spec_dir(directory)
 files = findfile( direct, COUNT = n)
 if n EQ 0 then begin                                      ;Any files found?
       message,'No files found on '+ direct, /CON
       return
 endif 

  good = where( strlen(files) GT 0, Ngood)
  if Ngood EQ 0 then message,'No FITS files found on '+ direct $
                 else files = files[good]

; Set output device according to keyword TEXTOUT or system variable !TEXTOUT

 if not keyword_set( TEXTOUT ) then textout= !TEXTOUT

 dir = 'dummy'
 num = 0

 get_lun,unit

 fmt1 = '(a,t20,a,t31,a,t42,a,t52,a,t71,a)'
 fmt2 = '(a,t22,a,t34,a,t45,a,t73,a)'

 for i = 0,n-1 do begin                           ;Loop over each FITS file
   fdecomp, files[i], disk, dir2, fname, ext     ;Decompose into disk+filenamOe
   if !VERSION.RELEASE GT '5.3' then begin 
      if (ext EQ 'gz') then compress = 1 else compress = 0
      openr, unit, files[i], /block, /binary, error = error, compress = compress    
   endif else openr, unit, files[i], /block, /binary, error = error
   if error LT 0 then goto, BADHD
   mrd_hread, unit, h, status, /silent
   if status LT 0 then goto, BADHD

   if exten GT 0 then begin 
         close,unit
         if !VERSION.RELEASE GT '5.3' then $
            openr, unit, files[i], /block, /binary, $
                       error = error, compress = compress else $    
            openr, unit, files[i], /block, /binary, error = error
         stat = fxmove(unit, exten, /silent)
         mrd_hread, unit, h1, status, /silent
         if status EQ 0 then h = [h1,h] 
    endif 

   keyword = strtrim( strmid(h,0,8),2 )       ;First 8 chars is FITS keyword
   value = strtrim( strmid(h,10,20),2 )        ;Chars 10-30 is FITS value

 l= where(keyword EQ 'NAXIS',Nfound)            ;Must have NAXIS keyword
    if Nfound GT 0 then naxis  = fix( value[ l[0] ] ) else goto, BADHD

 if naxis EQ 0 then naxisi = '0' else begin

 l = where( keyword EQ 'NAXIS1', Nfound)         ;Must have NAXIS1 keyword
    if Nfound gt 0 then naxis1  = long( value[l[0] ] ) else goto, BADHD 
    naxisi = strtrim( naxis1,2 )
 endelse

 if NAXIS GE 2 then begin
 l = where(keyword EQ 'NAXIS2', Nfound)          ;Must have NAXIS2 keyword
    if Nfound gt 0 then naxis2  = long(value[l[0]]) else goto, BADHD
    naxisi = naxisi + ' ' + strtrim( naxis2, 2 )
 endif

 if NAXIS GE 3 then begin
 l = where( keyword EQ 'NAXIS3', Nfound )          ;Must have NAXIS3 keyword
    if Nfound GT 0 then naxis3  = fix( value[l[0]] ) else goto, BADHD
    naxisi = naxisi + ' ' + strtrim( naxis3, 2 )
 endif

 if not keyword_set(NoTelescope) then begin
 l= where(keyword EQ 'TELESCOP',Nfound)         ;Search for TELESCOP keyword
 if Nfound EQ 0 then l = where(keyword EQ 'TELNAME',Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'OBSERVAT',Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'INSTRUME',Nfound)
    if Nfound GT 0 then begin 
          telescop = value[l[0]] 
          remchar,telescop,"'"
    endif else  telescop = '   ?      ' 
 endif

 l = where(keyword eq 'EXPTIME', Nfound)           ;Search for EXPTIME keyword
 if Nfound EQ 0 then l = where(keyword EQ 'TEXPTIME',Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'EXPOSURE',Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'INTEG',Nfound)
    if Nfound GT 0 then begin 
       exptim = float(value[l[0]]) 
       exptim = string(exptim, f = '(f7.1)')
    endif else exptim ='    ? '

 l = where(keyword EQ 'OBJECT',Nfound)            ;Search for OBJECT keyword
 if Nfound EQ 0 then l = where(keyword EQ 'TARGNAME',Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'TARGETID',Nfound)
    if Nfound GT 0 then begin 
       object = strtrim(strmid(h[l],10,30),2)
       remchar,object,"'"
   endif else object = '  ?     '

 l = where(keyword EQ 'DATE-OBS', Nfound)         ;Search for DATE-OBS keyword
 if Nfound EQ 0 then l = where(keyword EQ 'TDATEOBS', Nfound)
 if Nfound EQ 0 then l = where(keyword EQ 'DATE', Nfound) 
   if Nfound GT 0 then begin 
       obs = value[l[0]] 
       remchar, obs, "'"
       obs = strmid(strtrim(obs,2),0,10)
   endif else obs = '  ?     '


 num = num + 1
 if num EQ 1 then begin                 ;Print output header

    textopen, 'fitsdir', TEXTOUT=textout,/STDOUT 
    printf,!TEXTUNIT, f = '(a,/)', 'FITS File Directory ' + systime()
    if keyword_set(NoTelescope) then printf, !TEXTUNIT, $              
' NAME                 SIZE      DATE-OBS            OBJECT            EXPTIM' $
    else printf,!TEXTUNIT, $
 ' NAME                SIZE     DATE-OBS   TELESCOP  OBJECT             EXPTIM' 
endif

 if dir2 NE dir then begin                  ;Has directory changed?   
       if disk+dir2 EQ '' then cd,current=dir else dir = dir2
       printf, !TEXTUNIT,format='(/a/)', disk + dir + filename+'.'+ext  
       dir = dir2                                  ;Save new directory
 endif                   

 fname = strmid( fname, 0, namelen ) 
 object = strmid( object, 0, objectlen )


  if keyword_set( NOTELESCOPE) then $ 
  printf,!textunit,f=fmt2, $
      fname, naxisi, obs, object, exptim  else  $
  printf,!textunit,f=fmt1, $
      fname, naxisi, obs, telescop, object, exptim 

 BADHD:  

 close,unit
 endfor
 DONE: 
 free_lun, unit
 if num GT 0 then textclose, TEXTOUT=textout else begin 
           message,'No valid FITS files found on '+ spec_dir(direct),/CON
           return
 endelse 

 return      ;Normal return   
 end
