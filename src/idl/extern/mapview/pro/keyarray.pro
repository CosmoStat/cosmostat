Function KeyArray, instr
;+
; NAME:
;	KeyArray
; PURPOSE:
;	Breaks a string containing a number of comma-delimited elements
;	into an array of strings.
; CALLING SEQUENCE:
;	array = KeyArray(instr)
; INPUTS:
;	instr - The input string.
; OUTPUTS:
;	array - The output array.
; COMMENTS:
;	If the string contains surrounding parantheses or square brackets,
;	these are discarded.  Leading and trailing whitespace is discarded
;	from each array element.
; EXAMPLE:
;       IDL> print, keyarray('(trucks, cars, planes) ')
;               ===> 'trucks' 'cars' 'planes'
   
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 July 1999.
;       Use STRSPLIT() if since V5.3  W. Landsman  Dec 2002
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  array = KeyArray(instr)'
;
;			Remove the parentheses/brackets used to surround
;			the array of values.
;
str = strtrim(instr, 2)
;
bpar = ["[", "("]
epar = ["]", ")"]
nq   = n_elements(bpar)
i    = 0
done = 0
While ((done EQ 0) AND (i LT nq)) Do Begin
    pos1 = strpos(str, bpar[i], 0);
    If (pos1 EQ 0) Then Begin
	pos2 = strpos(str, epar[i], pos1 + 1)
	If (pos2 EQ (strlen(str) - 1)) Then Begin
	    str  = strmid(str, 1, (pos2 - 1))
	    done = 1
	EndIf
    EndIf
    i = i + 1
EndWhile
;
;			Break up the string into an array.
;
 if !VERSION.RELEASE GE '5.3' then array = strsplit(str,',',/ext) else begin
pos = strpos(str, ',')
If (pos LT 0) Then Begin
    array = [strtrim(str, 2)]
EndIf Else Begin
    array = [strtrim(strmid(str, 0, pos), 2)]
    str = strmid(str, (pos + 1), strlen(str))
EndElse
;
While (pos GE 0) Do Begin
    pos = strpos(str, ',')
    If (pos LT 0) Then Begin
        array = [array, strtrim(str, 2)]
    EndIf Else Begin
        array = [array, strtrim(strmid(str, 0, pos), 2)]
        str = strmid(str, (pos + 1), strlen(str))
    EndElse
EndWhile
EndElse
;
Return, array
End
