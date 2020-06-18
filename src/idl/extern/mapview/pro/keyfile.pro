Pro KeyFile, file, keyvals, Comment=comm, Continue=cont, Separator=sep
;+
; NAME:
;	KeyFile
; PURPOSE:
;	Reads a file containing keyword/value pairs into an array of
;	structures, each element of which contains a keyword/value pair.
; CALLING SEQUENCE:
;	KeyFile, file, keyvals
; INPUTS:
;	file    - The name of the file.
; OUTPUTS:
;	keyvals - The array of keyword/value pairs.  Each element has the form:
;			Key:    The keyword string.
;			Value:  The value string.
; OPTIONAL INPUT KEYWORDS:
;	Comment   - The comment character used in the file.  Defaults to '#'.
;	Continue  - The value line continuation character.  If a value ends in
;	            this character then additional lines are read and appended
;	            to the end of the value until a line is read that does not
;	            end in this character.  Defaults to 7.
;	Separator - The character used to separate the keyword and the value.
;	            Defaults to '='.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 July 1999.
;-
on_error, 2
;
;		Check arguments.
;
If (n_params() LT 2) Then message, 'Syntax: KeyFile, file, keyvals'
;
If (n_elements(comm) LE 0) Then comm = '#'
If (n_elements(cont) LE 0) Then cont = '&'
If (n_elements(sep)  LE 0) Then sep  = '='
;
;		Initialize the output array.
;
ele     = {Key:'', Value:''}
keyvals = [ele];
;
;		Open the file.
;
openr, unit, file, /Get_Lun
;
;		Read and parse the contents of the file.
;
While (Not Eof(unit)) Do Begin
;
;			Read the next line, looking for the next
;			keyword/value pair.
;
    line = TextRead(unit, Comment=comm, /Whitespace)
    If (strlen(line) GT 0) Then Begin
;
;			Parse out the keyword and value.
;
        KeyParse, line, key, value, Separator=sep
        If (strlen(key[0]) GT 0) Then Begin
;
;			Deal with line continuation characters.
;
            value = value[0]
            pos   = strpos(value, cont)
	    While ((pos + 1) EQ strlen(value)) Do Begin
	        value = strmid(value, 0, pos) $
		      + TextRead(unit, Comment=comm, /Whitespace)
		pos   = strpos(value, cont)
	    EndWhile
;
;			Add the keyword and value to the array.
;
            ele.Key   = key[0]
            ele.Value = value
            keyvals   = [keyvals, ele]
        EndIf
    EndIf
EndWhile
;
;			Close the file
;
free_lun, unit
;
;			Drop the first array element.
;
If (n_elements(keyvals) GT 1) Then keyvals = keyvals[1:*]
;
Return
End
