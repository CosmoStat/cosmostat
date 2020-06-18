Pro sxaddpar, Header, Name, Value, Comment, Location, before=before, $
                                       after=after , format=format, pdu = pdu
;+
; NAME:
;	SXADDPAR
; PURPOSE:
;	Add or modify a parameter in a FITS header array.
;
; CALLING SEQUENCE:
;	sxaddpar, Header, Name, Value, [ Comment,  Location,
;				BEFORE =, AFTER = , FORMAT= , /PDU]
;
; INPUTS:
;	Header = String array containing FITS or STSDAS header.    The
;		length of each element must be 80 characters.    If not 
;		defined, then SXADDPAR will create an empty FITS header array.
;
;	Name = Name of parameter. If Name is already in the header the value 
;		and possibly comment fields are modified.  Otherwise a new 
;		record is added to the header.  If name = 'HISTORY' then the 
;		value will be added to the record without replacement.  In 
;		this case the comment parameter is ignored.
;
;	Value = Value for parameter.  The value expression must be of the 
;		correct type, e.g. integer, floating or string.  String values
;		 of 'T' or 'F' are considered logical values.
;
; OPTIONAL INPUT PARAMETERS:
;	Comment = String field.  The '/' is added by this routine.  Added 
;		starting in position 31.    If not supplied, or set equal to 
;		'', then the previous comment field is retained (when found) 
;
;	Location = Keyword string name.  The parameter will be placed before the
;		location of this keyword.    This parameter is identical to
;		the BEFORE keyword and is kept only for consistency with
;		earlier versions of SXADDPAR.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;	BEFORE	= Keyword string name.  The parameter will be placed before the
;		location of this keyword.  For example, if BEFORE='HISTORY'
;		then the parameter will be placed before the first history
;		location.  This applies only when adding a new keyword;
;		keywords already in the header are kept in the same position.
;
;	AFTER	= Same as BEFORE, but the parameter will be placed after the
;		location of this keyword.  This keyword takes precedence over
;		BEFORE.
;
;	FORMAT	= Specifies FORTRAN-like format for parameter, e.g. "F7.3".  A
;		scalar string should be used.  For complex numbers the format
;		should be defined so that it can be applied separately to the
;		real and imaginary parts.  If not supplied then the default is
;               'G19.12' for double precision, and 'G14.7' for floating point.
;
;	/PDU    = specifies keyword is to be added to the primary data unit
;		header. If it already exists, it's current value is updated in
;		the current position and it is not moved.
; OUTPUTS:
;	Header = updated FITS header array.
;
; EXAMPLE:
;       Add a keyword 'TELESCOP' with the value 'KPNO-4m' and comment 'Name
;       of Telescope' to an existing FITS header h.
;
;       IDL> sxaddpar, h, 'TELESCOPE','KPNO-4m','Name of Telescope'
; NOTES:
;       The functions SXADDPAR() and FXADDPAR() are nearly identical, with the
;       major difference being that FXADDPAR forces required FITS keywords
;       BITPIX, NAXISi, EXTEND, PCOUNT, GCOUNT to appear in the required order
;       in the header.    There is no particular reason for having two nearly 
;       identical procedures, but both are too widely used to drop either one.
;
; RESTRICTIONS:
;	Warning -- Parameters and names are not checked
;		against valid FITS parameter names, values and types.
;
; MODIFICATION HISTORY:
;	DMS, RSI, July, 1983.
;	D. Lindler Oct. 86  Added longer string value capability
;	Converted to NEWIDL  D. Lindler April 90
;	Added Format keyword, J. Isensee, July, 1990
;	Added keywords BEFORE and AFTER. K. Venkatakrishna, May '92
;	Pad string values to at least 8 characters   W. Landsman  April 94
;	Aug 95: added /PDU option and changed routine to update last occurence
;		of an existing keyword (the one SXPAR reads) instead of the
;		first occurence.
;	Comment for string data can start after column 32 W. Landsman June 97
;	Make sure closing quote supplied with string value  W. Landsman  June 98
;	Converted to IDL V5.0    W. Landsman   June 98
;	Increase precision of default formatting of double precision floating
;		point values.	C. Gehman, JPL	September 1998
;	Mar 2000, D. Lindler, Modified to use capital E instead of lower case
;		e for exponential formats.
;       Apr 2000, Make user-supplied format upper-case  W. Landsman April 2000
;	
;-
 if N_params() LT 3 then begin             ;Need at least 3 parameters
      print,'Syntax - Sxaddpar, Header, Name,  Value, [Comment, Postion'
      print,'                      BEFORE = ,AFTER = , FORMAT = ]'
      return
 endif

; Define a blank line and the END line

 ENDLINE = 'END' +string(replicate(32b,77))	;END line
 BLANK = string(replicate(32b,80))	       ;BLANK line
;
;  If Location parameter not defined, set it equal to 'END     '
;
 if ( N_params() GT 4 ) then loc = strupcase(location) else $
 if keyword_set( BEFORE) then loc = strupcase(before) else $
 if keyword_set( AFTER)  then loc = strupcase(after) else $
 if keyword_set( PDU) then loc = 'BEGIN EX' else $
                             loc = 'END'

 while strlen(loc) lt 8 do loc = loc + ' '

 if N_params() lt 4 then comment = ''      ;Is comment field specified?

 n = N_elements(header)	                 ;# of lines in FITS header
 if (n EQ 0) then begin	                 ;header defined?
	  header=strarr(10)              ;no, make it.
	  header[0]=ENDLINE
	  n=10
 endif else begin
	  s = size(header)               ;check for string type
	      if (s[0] ne 1) or (s[2] ne 7) then $
		  message,'FITS Header (first parameter) must be a string array'
 endelse

;  Make sure Name is 8 characters long

	nn = string(replicate(32b,8))	;8 char name
	strput,nn,strupcase(name) ;insert name

;  Extract first 8 characters of each line of header, and locate END line

 keywrd = strmid(header,0,8)                 ;Header keywords
 iend = where(keywrd eq 'END     ',nfound)
 if nfound eq 0 then header[0]=ENDLINE	    ;no end, insert at beginning
 iend = iend[0] > 0                          ;Make scalar

 if nn eq 'HISTORY ' then begin             ;add history record?
      if iend lt (n-1) then $
		 header[iend+1]=ENDLINE $ ;move end up
      else begin
		 header = [header,replicate(blank,5)] ;yes, add 5.
		 header[n] = ENDLINE
      endelse
      newline = blank
      strput,newline,nn+string(value),0
      header[iend] = newline		 ;add history rec.
      return
 endif		;history

; Find location to insert keyword.   Save the existing comment if user did
; not supply a new one.   Comment starts after column 32 for numeric data,
; after the slash (but at least after column 20) for string data. 

 ipos  = where(keywrd eq nn,nfound)
 if nfound gt 0 then begin
         i = ipos[nfound-1]
         if comment eq '' then begin             ;save comment?
	 if strmid(header[i],10,1) NE "'" then $
		 comment=strmid(header[i],32,48) else begin
		 slash = strpos(header[i],'/', 20)  
		 if slash NE -1 then $
			comment =  strmid(header[i], slash+1, 80) else $
			comment = string(replicate(32B,80))
		endelse
	endif
         goto, REPLACE    
 endif

 if loc ne '' then begin
          iloc =  where(keywrd eq loc,nloc)
          if nloc gt 0 then begin
             i = iloc[0]
             if keyword_set(after) and loc ne 'HISTORY ' then i = i+1 < iend 
	     if i gt 0 then header=[header[0:i-1],blank,header[i:n-1]] $
		        else header=[blank,header[0:n-1]]
             goto, REPLACE  
	  endif
 endif

; At this point keyword and location parameters were not found, so a new
; line is added at the end of the FITS header

	if iend lt (n-1) then begin	;Not found, add more?
		header[iend+1] = ENDLINE	;no, already long enough.
		i = iend		;position to add.
	   endif else begin		;must lengthen.
		header = [header,replicate(blank,5)] ;add an element on the end
		header[n]=ENDLINE		;save "END"
		i =n-1			;add to end
	end

; Now put value into keyword at line i

REPLACE:    
	h=blank			;80 blanks
	strput,h,nn+'= '	;insert name and =.
	apost = "'"	        ;quote a quote
	type = size(value)	;get type of value parameter
	if type[0] ne 0 then $
		message,'Keyword Value (third parameter) must be scalar'

	case type[1] of		;which type?

7:	begin
	  upval = strupcase(value)	;force upper case.
	  if (upval eq 'T') or (upval eq 'F') then begin
		strput,h,upval,29  ;insert logical value.
	    end else begin		;other string?
		if strlen(value) gt 18 then begin	;long string
		    strput, h, apost + strmid(value,0,68) + apost + $
			' /' + comment,10
		    header[i] = h
		    return
		endif
		strput, h, apost + value,10       ;insert string val
		strput, h, apost, 11 + (strlen(value)>8)   ;pad string vals
	  endelse                                          ;to at least 8 chars
	  endcase

5:	BEGIN
	IF (N_ELEMENTS(format) EQ 1) THEN $		; use format keyword
	    v = string(value, FORMAT='('+strupcase(format)+')') $
	ELSE v = STRING(value, FORMAT='(G19.12)')
	s = strlen(v)					; right justify
	strput, h, v, (30-s)>10
	END

 else:  begin
	if (N_elements(format) eq 1) then $            ;use format keyword
	    v = string(value, FORMAT='('+strupcase(format)+')' ) else $
	    v = strtrim(strupcase(value),2)      
	    			      ;convert to string, default format
	s = strlen(v)                 ;right justify
        strput,h,v,(30-s)>10          ;insert
	end
 endcase

 strput,h,' /',30	;add ' /'
 strput, h, comment, 32	;add comment
 header[i] = h		;save line

 return
 end
