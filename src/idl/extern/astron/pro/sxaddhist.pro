pro sxaddhist,history,header,pdu=pdu
;+
; NAME:
;	SXADDHIST                           
; PURPOSE:
;	Procedure to add history line(s) to a FITS header
;
; CALLING SEQUENCE
;	sxaddhist, history, header, [ /PDU ]
;
; INPUTS:
;	history - string or string array containing history line(s)
;		to add to the header
;	header - string array containing the FITS header
;
; KEYWORD INPUTS:
;	/PDU - if specified, the history will be added to the primary
;		data unit header, (before the line beginning BEGIN EXTENSION...)
;		Otherwise, it will be added to the end of the header
; OUTPUTS:
;	header - unpdated header
;
; EXAMPLES:
;	sxaddhist, 'I DID THIS', header
;
;	hist = strarr(3)
;	hist[0] = 'history line number 1'
;	hist[1[ = 'the next history line'
;	hist[2] = 'the last history line'
;	sxaddhist, hist, header
;
; HISTORY:
;	D. Lindler  Feb. 87
;	April 90  Converted to new idl  D. Lindler
;	Put only a single space after HISTORY   W. Landsman  November 1992
;	Aug. 95	  Added PDU keyword parameters
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;--------------------------------------------------------------------
 	On_error,2

 	if N_params() LT 2 then begin
 	   print, ' Syntax - SXADDHIST, hist, header, [ /PDU] '
 	   return
 	endif

; Check input parameters

 	s = size(history) & ndim = s[0] & type = s[ndim+1]
 	if type NE 7 then message, $
 	    'Invalid history lines specified; must be a string or string array'
 
 	nadd = N_elements(history)	      ;Number of lines to add

 	s = size(header) & ndim2 = s[0] & type = s[ndim2+1]
 	if (ndim2 NE 1) or (type NE 7) then message, $
		'Invalid FITS header supplied; header must be a string array'

	nlines = N_elements(header)	      ;Number of lines in header

; Find END statement of FITS header
        
 	endline = where( strtrim(strmid(header,0,8),2) EQ 'END' )
 	n = endline[0]
 	if n LT 0 then message, $
		    'Invalid FITS header array, END keyword not found'

 	blank = string( replicate(32b,80) )
;
; if /PDU find beginning of the extension header and make room for the
; history
;
	n1 = n		;position to insert
	if keyword_set(PDU) then begin
	    extline = where( strtrim(strmid(header,0,8),2) EQ 'BEGIN EX' )
	    n_ext = extline[0]
	    if n_ext gt 1 then n1 = n_ext
	end
;
; make room in the header
;
	if n1 eq 0 then header = [replicate(blank,nadd),header[n1:n]] else $
		header = [header[0:n1-1],replicate(blank,nadd),header[n1:n]]

; Add history records to header starting at position N1

 	for i = 0, nadd-1 do begin
	
		newline = blank
		strput, newline, 'HISTORY ' + history[i]
		header[n1+i] = newline

 	endfor
 return
 end
