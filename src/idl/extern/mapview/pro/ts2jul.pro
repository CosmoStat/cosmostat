Function TS2Jul, ts, RefGMT=refgmt
;+
; NAME:
;	TS2Jul
; PURPOSE:
;	Converts an array of time stamps into Julian dates.
; CALLING SEQUENCE:
;	juls = TS2Jul(ts)
; INPUTS:
;	ts   - The array of data.
; OUTPUTS:
;	juls - The array of WMAP reduced Julian dates (standard Julian date 
;	       minus 2,450,000).
; OPTIONAL INPUT KEYWORD:
;	RefGMT - The reference GMT, defining timestamp [0,0].  If not
;	         supplied the DefTSRef is used.
; EXAMPLE:
;      IDL> print,TS2Jul([ [513, 677980204],[512, 67880204] ] )
;           2604.2847       2602.5786
; PROCEDURES USED:
;       DefTSRef(), GMT2JUL()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon ITSS, 06 October 1999
;-
On_error, 2
;
;			Check arguments.
;
If (n_elements(ts) LE 0) Then message, 'Syntax:  juls = TS2Juls(ts)'
;
If (n_elements(refgmt) LE 0) Then refgmt = DefTSRef()
;
;			Perform the conversion.
;
juls = gmt2jul(refgmt)				$
     + double(ts[0,*])				$
     + (double(ts[1,*]) / 8.64d+08)
;
Return, reform(juls)
End
