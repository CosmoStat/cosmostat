Function Jul2TS, juls, RefGMT=refgmt
;+
; NAME:
;	Jul2TS
; PURPOSE:
;	Converts an array of Julian dates into time stamps.
; CALLING SEQUENCE:
;	ts = Jul2TS(juls)
; INPUTS:
;	juls - The array of MAP reduced Julian dates (standard Julian date 
;	       minus 2,450,000).
; OUTPUTS:
;	ts   - The array of time stamps.
; KEYWORDS:
;	RefGMT - The reference GMT, defining timestamp [0,0].  If not
;	         supplied the DefTSRef is used
; PROCEDURE CALLED:
;       GMT2JUL(), DefTSRef()                      ;MAP Library.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, SSAI, 30 September 2002.
;-
on_error, 2
;
;			Check arguments.
;
If (n_elements(juls) LE 0) Then message, 'Syntax:  ts = Jul2TS(juls)'
;
If (n_elements(refgmt) LE 0) Then refgmt = DefTSRef()
;
;			Perform the conversion.
;
jr = gmt2jul(refgmt)
j = double(juls[*]) - jr
d = long(j)
t = (j - double(d)) * 8.64d+08
ts = [[d], [long(t)]]
;
Return, transpose(ts)
End
