Pro Jul2Tel, jul, tel, RefGMT=refgmt
;+
; NAME:
;	Jul2Tel
; PURPOSE:
;	Converts a MAP Reduced Julian day into a telemetry structure timestamp.
; CALLING SEQUENCE:
;	Jul2Tel, jul, tel
; INPUTS:
;	jul - The MAP Reduced Julian day(s) to convert.
; OUTPUTS:
;	tel - The telemetry structure(s) to contain the results.  These 
;	      structures must already exist and have Day and Time tags
; OPTIONAL INPUT KEYWORD:
;	RefGMT - The reference GMT, defining timestamp [0,0].  If not
;	         supplied the DefTSRef is used.
; PROCEDURES USED:
;       DefTSRef(), GMT2JUL()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon ITSS, 21 April 1999
;-
On_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, 'Syntax:  Jul2Tel, jul, tel'
;
If (n_elements(refgmt) LE 0) Then refgmt = DefTSRef()
;
;			Perform the conversion.
;
j = jul - gmt2jul(refgmt)
tel.day  = long(j)
tel.time = long(0.49d0 + ((j - tel.day) * 8.64d+08))
;
neg = where(tel.time LT 0)
If (neg[0] GE 0) Then Begin
       tel[neg].day  = tel[neg].day  - 1L
       tel[neg].time = tel[neg].time + long(8.64d+08)
EndIf
;
Return
End
