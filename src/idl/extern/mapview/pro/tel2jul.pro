Function Tel2Jul, tel, RefGMT=refgmt
;+
; NAME:
;	Tel2Jul
; PURPOSE:
;	Converts the time stamps in an array of data (structure array with
;	Day and Time members) into Julian dates.
; CALLING SEQUENCE:
;	juls = Tel2Jul(tel, [ RefGMT = ])
; INPUTS:
;	tel  - The array of data.
; OUTPUTS:
;	juls - The array of MAP reduced Julian dates (standard Julian date 
;	       minus 2,450,000).
; OPTIONAL INPUT KEYWORD:
;	RefGMT - The reference GMT, defining timestamp [0,0].  If not
;	         supplied the DefTSRef is used.
; EXAMPLE:
;       Translates the time in the Day 222 Year 2002 time-ordered data into
;       into Julian dates 
;
;       IDL> file ='MAP_tod_20022202357_20022212357.fits'
;       IDL>  fits_read_tod, file, tod
;       IDL> jd = tel2jul( tod.sci) 
; PROCEDURES USED:
;       DefTSRef(), GMT2JUL()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon ITSS, 01 April 1999
;-
on_error, 2
;
;			Check arguments.
;
n = n_elements(tel)
If (n LE 0) Then message, 'Syntax:  juls = Tel2Juls(tel)'
;
If (n_elements(refgmt) LE 0) Then refgmt = DefTSRef()
;
;			Perform the conversion.
;
juls = gmt2jul(refgmt)				$
     + double(tel.day)				$
     + (double(tel.time) / 8.64d+08)
;
Return, juls
End
