Pro KeyParse, line, key, value, Separator=sep
;+
; NAME:
;	KeyParse
; PURPOSE:
;	Breaks a string into a keyword/value pair.
; CALLING SEQUENCE:
;	KeyParse, line, key, value
; INPUTS:
;	line  - The string to break up.
; OUTPUTS:
;	key   - The keyword portion of the string.
;	value - The value portion of the string.
; KEYWORDS:
;	Separator - The character used to separate the two parts.
;	            Defaults to '='.
; COMMENTS:
;	Leading and trailing whitespace is removed from both portions
;	of the string.  If the value has leading and trailing quotes,
;	these are removed (rather, the outermost quote-pair is removed).
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 July 1999.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 3) Then message, 'Syntax: KeyParse, line, key, value'
;
If (n_elements(sep) LE 0) Then sep = '='
;
;		For each line...
;
nl    = n_elements(line)
key   = replicate('', nl)
value = replicate('', nl)
;
quot  = ["'", '"']
nq    = n_elements(quot)
;
For j = 0, (nl - 1) Do Begin
;
;			Look for the separator and break the line up into
;			keyword/value parts.  If the separator is not found
;			then return empty key and value strings for the line.
;
    pos1 = strpos(line[j], sep)
    If (pos1 GT 0) Then Begin
;
        key[j]   = strtrim(strmid(line[j], 0, pos1), 2)
        value[j] = strtrim(strmid(line[j], (pos1 + 1), strlen(line[j])), 2)
;
;			Strip off the leading and trailing quotes.
;
	i    = 0
	done = 0
	While ((done EQ 0) AND (i LT nq)) Do Begin
	    pos1 = strpos(value[j], quot[i], 0);
	    If (pos1 EQ 0) Then Begin
		pos2 = strpos(value[j], quot[i], pos1 + 1)
		If (pos2 EQ (strlen(value[j]) - 1)) Then Begin
		    value[j] = strmid(value[j], 1, (pos2 - 1))
		    done  = 1
		EndIf
	    EndIf
	    i = i + 1
	EndWhile

    EndIf
EndFor
;
Return
End
