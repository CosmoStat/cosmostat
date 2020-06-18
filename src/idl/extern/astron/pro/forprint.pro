pro forprint, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
      v15,v16,v17,v18,TEXTOUT = textout, FORMAT = format, SILENT = SILENT, $ 
      STARTLINE = startline, NUMLINE = numline
;+
; NAME:
;	FORPRINT
; PURPOSE
;	Print a set of vectors by looping over each index value.
; EXPLANATION:
;	If W and F are equal length vectors, then the statement
;		IDL> forprint, w, f   
;	is equivalent to 
;		IDL> for i = 0L, N_elements(w)-1 do print,w[i],f[i]    
;
; CALLING SEQUENCE:
;	forprint, v1,[ v2, v3, v4,....v18, FORMAT = , TEXTOUT = ,STARTLINE =,
;					   NUMLINE =, /SILENT ] 
;
; INPUTS:
;	V1,V2,...V18 - Arbitary IDL vectors.  If the vectors are not of
;		equal length then the number of rows printed will be equal
;		to the length of the smallest vector.   Up to 18 vectors
;		can be supplied.
;
; OPTIONAL KEYWORD INPUTS:
;
;	TEXTOUT - Controls print output device, defaults to !TEXTOUT
;
;		textout=1	TERMINAL using /more option
;		textout=2	TERMINAL without /more option
;		textout=3	<program>.prt
;		textout=4	laser.tmp
;		textout=5      user must open file
;		textout = filename (default extension of .prt)
;		textout=7	Append to <program>.prt file if it exists
;
;	FORMAT - Scalar format string as in the PRINT procedure.  The use
;		of outer parenthesis is optional.   Ex. - format="(F10.3,I7)"
;		This program will automatically remove a leading "$" from
;		incoming format statments. Ex. - "$(I4)" would become "(I4)".
;	STARTLINE - Integer scalar specifying the first line in the arrays
;		to print.   Default is STARTLINE = 1, i.e. start at the
;		beginning of the arrays.
;	SILENT - Normally, with a hardcopy output (TEXTOUT > 2), FORPRINT will
;		add a time stamp to the output file.    If the SILENT keyword
;		is set and non-zero, then this time stamp is suppressed.
; OUTPUTS:
;	None
; SYSTEM VARIABLES:
;	If keyword TEXTOUT is not used, the default is the nonstandard 
;	keyword !TEXTOUT.    If you want to use FORPRINT to write more than 
;	once to the same file, or use a different file name then set 
;	TEXTOUT=5, and open and close then file yourself (see documentation 
;	of TEXTOPEN for more info).
;	
;	One way to add the non-standard system variables !TEXTOUT and !TEXTUNIT
;	is to use the procedure ASTROLIB
; EXAMPLE:
;	Suppose W,F, and E are the wavelength, flux, and epsilon vectors for
;	an IUE spectrum.   Print these values to a file 'output.dat' in a nice 
;	format.
;
;	IDL> fmt = '(F10.3,1PE12.2,I7)'
;	IDL> forprint, F = fmt, w, f, e, TEXT = 'output.dat'
;
; PROCEDURES CALLED:
;	DATATYPE(), TEXTOPEN, TEXTCLOSE
; REVISION HISTORY:
;	Written    W. Landsman             April, 1989
;	Keywords textout and format added, J. Isensee, July, 1990
;	Made use of parenthesis in FORMAT optional  W. Landsman  May 1992
;	Added STARTLINE keyword W. Landsman    November 1992
;	Set up so can handle 18 input vectors. J. Isensee, HSTX Corp. July 1993
;	Handle string value of TEXTOUT   W. Landsman, HSTX September 1993
;	Added NUMLINE keyword            W. Landsman, HSTX February 1996
;	Added SILENT keyword             W. Landsman, RSTX, April 1998
;	Converted to IDL V5.0            W. Landsman, RSTX, April, 1998
;-            
  On_error,2                               ;Return to caller

  npar = N_params()
  if npar EQ 0 then begin
      print,'Syntax - FORPRINT, v1, [ v2, v3,...v18, FORMAT =, /SILENT, '
      print,'                         STARTLINE = , NUMLINE =, TEXTOUT =]'
      return
  endif

  if not keyword_set( STARTLINE ) then startline = 1l else $
         startline = startline > 1l 

  fmt="F"                 ;format flag
  istart = 2
  npts = N_elements(v1)

  if ( npts EQ 0 ) then message,'ERROR - Parameter 1 is not defined'

;  Remove "$" sign from format string and append parentheses if not 
;  already present

  if N_elements( format ) EQ 1 then begin

     fmt = "T"                                 ;format present
     frmt = format            
     if strmid(frmt,0,1) eq '$' then $
          frmt = strmid(frmt,1,strlen(frmt)-1) ;rem. '$' from format if present

     if strmid(frmt,0,1) NE '(' then frmt = '(' + frmt
     if strmid( frmt,strlen(frmt)-1,1) NE ')' then frmt = frmt + ')'

  endif

  if npar GT 1 then begin         ;Get number of elements in smallest array

      for i = istart, npar do begin 
          tst = execute('npts = min([ npts, N_elements(v'+strtrim(i,2)+')])')
          if npts EQ 0 then $
              message,'ERROR - Parameter ' + strtrim(i,2) + ' is not defined'
      endfor

  endif

  if keyword_set(NUMLINE) then npts = (startline + numline) < npts

  str = 'v' + strtrim( istart-1,2) + '[i]'
  if npar GT 1 then $
       for i = istart, npar do str = str + ',v' + strtrim(i,2) + '[i]'

; Use default output dev.
   if not keyword_set( TEXTOUT ) then textout = !TEXTOUT 
   if datatype( textout) EQ 'STR' then text_out = 6  $      ;make numeric
                                  else text_out = textout

   textopen,'FORPRINT',TEXTOUT=textout
   if ( text_out GT 2 ) and (not keyword_set(SILENT)) then $
	printf,!TEXTUNIT,'FORPRINT: ',systime()

   if fmt EQ "F" then begin            ;Use default formats

      for i = startline-1, npts-1 do begin 

          test = execute('printf,!TEXTUNIT,' + str) 
          if ( text_out EQ 1 ) then $                  ;Did user press 'Q' key?
               if !ERR EQ 1 then goto,DONE

      endfor

   endif else begin                    ;User specified format
 
      for i = startline-1, npts-1 do begin 

         test = execute( 'printf, !TEXTUNIT,  FORMAT=frmt,' + str ) 
         if  ( text_out EQ 1 ) then  $
               if !ERR EQ 1 then goto,DONE

      endfor

  endelse

DONE: 
  textclose, TEXTOUT = textout          ;Close unit opened by TEXTOPEN

  return
  end
