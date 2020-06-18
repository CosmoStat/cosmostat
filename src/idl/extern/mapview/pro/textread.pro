Function TextRead, file, status, Comment=comm, Whitespace=wspace, Zero=zer
;+
; NAME:
;	TextRead
; PURPOSE:
;	Read a non-empty record/line from an open text file.  Lines
;	are read until a nonzero length string is read.
; CALLING SEQUENCE:
;	line = TextRead( file [, status] )
; INPUTS:
;	file   - The unit number for an open text file.
; OUTPUTS:
;	status - A status code: 0=success, -1=EOF.
; RETURNED:
;	line   - The just-read line.
; OPTIONAL INPUT KEYWORDS:
;	Comment    - If supplied, this character/string indicates the start of a
;	             comment in a line.  Text starting at the comment string
;	             is removed from the line before it is returned.
;	/Whitespace - This routine ordinarily does not remove leading and
;	             trailing whitespace, and such white space will lead the
;	             routine to consider a line to not be empty.  If this
;	             keyword is present and nonzero then leading and trailing
;	             whitespace is removed.
;	/Zero       - If present and nonzero then zero-length strings are
;	             acceptable.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 02 July 1999.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax: line = TextRead, file [, status]'
;
If (n_elements(comm) LE 0) Then cchar = '' $
                           Else cchar = comm[0]
;
wflg = keyword_set(wspace)
zflg = keyword_set(zer)
;
;			Read the desired line.
;
status = 0
line   = ''
done   = 0
While (done EQ 0) Do Begin
;
;				Read the next line out of the file.
;
	If (Eof(file)) Then Begin
		status = -1
		Return, line
	EndIf
	readf, file, line
;
;				Strip off any comments.
;
	If (strlen(cchar) GT 0) Then Begin
		pos = strpos(line, cchar)
		If (pos GE 0) Then Begin
			sq1 = strpos(line, "'")
			If (sq1 GE 0) Then sq2 = strpos(line, "'", sq1 + 1) $
			              Else sq2 = -1
			dq1 = strpos(line, "'")
			If (dq1 GE 0) Then dq2 = strpos(line, "'", dq1 + 1) $
			              Else dq2 = -1
			If ((sq1 LT pos) AND (pos LT sq2)) Then pos = -1
			If ((dq1 LT pos) AND (pos LT dq2)) Then pos = -1
		EndIf
		If (pos GE 0) Then Begin
			If (pos EQ 0) Then line = ''		$
			              Else line = strmid(line, 0, pos)
		EndIf
	EndIf
;
;				If a zero-length string is allowed or if a nonzero
;				length line was read then we are done.
;
	If (wflg) Then line = strtrim(line, 2)
	If ((zflg) OR (strlen(line) GT 0)) Then done = 1
;
EndWhile
;
Return, line
End
