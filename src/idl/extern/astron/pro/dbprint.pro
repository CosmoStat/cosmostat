pro dbprint,list,items, FORMS=forms, TEXTOUT=textout, NOHeader = noheader
;+
; NAME:
;     DBPRINT
; PURPOSE:
;     Procedure to print specified items from a list of database entries
;
; CALLING SEQUENCE:     
;     dbprint, list, [items, FORMS= , TEXTOUT= , /NoHeader]  
;
; INPUTS:
;     list  - list of entry numbers to be printed, vector or scalar 
;               if list = -1, then all entries will be printed.
;               An error message is returned if any entry number is larger
;               than the number of entries in the database
;
; OPTIONAL INPUT-OUTPUT:
;     items - items to be printed, specified in any of the following ways:
;
;               form 1  scalar string giving item(s) as list of names
;                       separated by commas
;               form 2  string array giving list of item names
;               form 3  string of form '$filename' giving name
;                       of text file containing items (one item per
;                       line)
;               form 4  integer scalar giving single item number or
;                         integer vector list of item numbers
;               form 5  Null string specifying interactive selection.   This
;                       is the default if 'items' is not supplied
;               form 6  '*'     select all items, printout will be in
;                       table format. 
;
;            If items was undefined or a null string on input, then
;            on output it will contain the items interactively selected.
;
; OPTIONAL INPUT KEYWORDS:
;       FORMS - The number of printed lines per page. If forms is not 
;               present, output assumed to be in PORTRAIT form, and 
;               a heading and 47 lines are printed on each page, with
;               a page eject between each page.  For LANDSCAPE form with
;               headings on each page, and a page eject between pages, set 
;               forms = 34.  For a heading only on the first page, and no
;               page eject, set forms = 0.   This is the default for output
;               to the terminal.
;
;       TEXTOUT - Integer (0-7) used to determine output device (see TEXTOPEN
;       for more info).  If not present, the !TEXTOUT system variable is used.
;               textout=0       Nowhere
;               textout=1       if a TTY then TERMINAL using /more option
;                                   otherwise standard (Unit=-1) output
;               textout=2       if a TTY then TERMINAL without /more option
;                                   otherwise standard (Unit=-1) output
;               textout=3       dbprint.prt (file)
;               textout=4       laser.tmp
;               textout=5       user must open file
;               textout=7      same as 3 but text is appended to <program>.prt
;               textout = filename   (default extension of .prt)
;
;       /NOHEADER - If this keyword is set, then the column headers will not
;               be printed
;
; EXAMPLE:
;       The following example shows how a multiple valued item DATAMAX can be 
;       printed as separate columns.   In the WFPC2 target database, DATAMAX
;       is an item with 4 values, one for each of the 4 chips
;
;       IDL> dbopen,'wflog'
;       IDL> dbprint,list,'entry,datamax(0),datamax(1),datamax(2),datamax(3)'
;
; SYSTEM VARIABLES:
;       Output device controlled by non-standard system varaible !TEXTOUT, if 
;       TEXTOUT keyword is not used.    
;
; NOTES:
;       Users may want to adjust the default lines_per_page value given at
;       the beginning of the program for their own particular printer.
;
; HISTORY:
;       version 2  D. Lindler  Nov. 1987 (new db format)
;       Test if user pressed 'Q' in response to /MORE W. Landsman  Sep 1991
;       Apply STRTRIM to free form (table) output    W. Landsman   Dec 1992
;       Test for string value of TEXTOUT         W. Landsman   Feb 1994
;       William Thompson, GSFC, 3 November 1994
;                       Modified to allow ZDBASE to be a path string.
;       W. Landsman, GSFC, July, 1997, Use CATCH to catch errors
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Removed STRTRIM in table format output to handle byte values April 1999
;       Fixed occasional problem when /NOHEADER is supplied   Sep. 1999
;       Only byteswap when necessary for improved performance  Feb. 2000
;       Change loop index for table listing to type LONG  W. Landsman Aug 2000
;-
;
 On_error,2                                ;Return to caller


 if N_params() EQ 0 then begin
       print,'Syntax - dbprint, list, items, [ FORMS = , TEXTOUT =, /NoHeader ]'
       return
 endif

 lines_per_page = 47                 ;Default # of lines per page
 zparcheck, 'DBPRINT', list, 1, [1,2,3,4,5], [0,1], 'Entry List Vector'

 catch, error_status
 if error_status NE 0 then begin 
       print,!ERR_STRING
       return
  endif


; Make list a vector

 nentry = db_info( 'ENTRIES', 0)
 if list[0] EQ -1 then list = lindgen(nentry) + 1 
 dbname = db_info( 'NAME', 0 )
 if !VERSION.OS NE 'vms' then dbname = strlowcase(dbname)

 if max(list) GT nentry then message, dbname + $
     ' entry numbers must be between 1 and ' + strtrim( nentry, 2 )
  nv = N_elements(list)                 ;number of entries requested

; No need for byteswapping if data is not external or it is a big endian machine

  convert= db_info('EXTERNAL')
  noconvert = 1 - convert[0]
  if convert[0] then noconvert = is_ieee_big() 

; Determine items to print

 if N_params() EQ 1 then begin

      file = find_with_def(dbname +'.items', 'ZDBASE')
      if file NE '' then items = '$' + file else items = '' 

 endif
 
 db_item, items, it, ivalnum, dtype, sbyte, numvals, nbytes
 numvals = numvals<1                    ;can't print vectors
 nvalues = db_item_info( 'NVALUES', it )        ;number of values in item
 qnumit = db_info( 'ITEMS' )                    ;number of items
 nitems = N_elements( it )                      ;number of items requested
 qnames = db_item_info( 'NAME', it )
 qtitle = db_info( 'TITLE', 0 )         ;data base title

; Open output text file

 if not keyword_set(TEXTOUT) then textout = !textout  ;use default output dev.

 textopen, dbname, TEXTOUT = textout
 if datatype(TEXTOUT) EQ 'STR' then text_out = 5 else text_out = textout
 if (nitems EQ qnumit)  then begin

; Create table listing of each item specified. -------------------------

 for i = 0L, nv-1 do begin
      dbrd, list[i], entry, noconvert = noconvert   ; read an entry.
      printf, !TEXTUNIT, ' '                        ; print  blank line.

; display name and value for each entry 

      for k = 0, qnumit-1  do begin
         ;
         ; only print entries of reasonable size... < 5 values in item.
         ;
         if ( nvalues[k] LT 5 ) then begin
            somvar =  dbxval(entry,dtype[k],nvalues[k],sbyte[k],nbytes[k]) 
            if dtype[k] EQ 1 then somvar=fix(somvar)
            printf,!textunit,k,') ',qnames[k], strtrim(somvar,2)
                                                        ;display name,value
         endif                                               
       endfor   ; k

    endfor      ; i

 printf,!textunit,' '                         ;Added 11/90
 
 end else begin

; get info on items

   formats = db_item_info( 'FORMAT', it )
   flen = db_item_info( 'FLEN', it )            ;field lengths
   nvals = db_item_info( 'NVALUES', it )        ;larger than one for vector items

; Set up format array

   form = '(' + strtrim(formats,2)      + ')'   ;remove blanks, and add paren

   linelength = total(flen) + nitems            ;length of output lines
   dash = byte('-') & dash = dash[0]
   dashes = ' '+string( replicate( dash, linelength ) )
;
   if not keyword_set( NoHeader) then begin

      title = string( replicate(byte(32), linelength>42) )
      strput, title, qtitle, (linelength-40)/2>1           ;center title

; Extract headers

    headers = db_item_info( 'HEADERS', it )
    c1 = strmid( headers,0,15 )
    c2 = strmid( headers,15,15 )
    c3 = strmid( headers,30,15 )

; Place value numbers for multiple valued items in h3
    for i = 0,nitems-1 do begin
          if nvals[i] GT 1 then $       ;multiple values?
             c3[i] = '(' + strtrim(string(ivalnum[i]),2) + ')'
    endfor        ;i

    h1 = dbtitle( c1,flen )
    h2 = dbtitle( c2,flen )
    h3 = dbtitle( c3,flen )

 endif

; Loop on entries

 if ( N_elements(forms) GT 0 ) then begin
        if ( forms GT 0 ) then pcount = forms $ ;lines per page
        else pcount = N_elements(list)          ;no page breaks
 endif else if text_out LE 2 then pcount = N_elements(list) $
      else pcount = lines_per_page                ;Portrait form default
 limit = pcount - 1

  for j = 0L, N_elements(list)-1 do begin

   if not keyword_set( NoHeader) then begin

        if pcount GT limit then begin           ;new page?
                pcount = 0
                if text_out GT 2 then $
                            printf,!textunit,string(byte(12))   $;eject
                       else printf,!textunit,' '
                printf,!textunit,title                  ;print title
                printf,!textunit,dashes                 ;print headings
                printf,!textunit,h1
                printf,!textunit,h2
                printf,!textunit,h3
                printf,!textunit,dashes
        endif

    endif
        dbrd, list[j], entry, noconvert = noconvert        ;read entry
        ;
        ; loop on items
        ;
        st = ''                                 ;output string
        for i = 0,nitems-1 do  begin

                val = dbxval(entry,dtype[i],numvals[i],sbyte[i],nbytes[i])
                if dtype[i] EQ 1 then val = fix(val)
                if dtype[i] EQ 7 then begin
                   b = byte(val)
                   bad = where(b EQ 0, nbad)
                   if nbad GT 0 then begin
                       b[bad] = 32b
                       val = string(b)
                   endif
                endif
                st = st+' ' + string(val,form[i])

        endfor

        printf, !TEXTUNIT, st                   ;print line
        if text_out EQ 1 then  $       ;Did user press 'Q' in /MORE ?
                if ( !ERR EQ 1 ) then return
        pcount = pcount+1            ;increment line counter
    end                              ; loop on entries

 endelse                             ; N_params > 1

; Clean up

 textclose, TEXTOUT = textout                   ;close text file

 return
 end
