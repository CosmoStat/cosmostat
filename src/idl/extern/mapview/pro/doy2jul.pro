Pro DOY2Jul, day, year, jul, Full=ful, Verbose=vrb
;+
; NAME:
;	DOY2Jul
; PURPOSE:
;	Determines the MAP Reduced Julian day from day-of-year and year.
; CALLING SEQUENCE:
;	DOY2Jul, day, year [, jul, /FULL, /VERBOSE]
; INPUTS:
;	day  - The day-of-year, scalar or vector
;	year - The year, same number of elements as day
; OUTPUTS:
;	jul  - The MAP Reduced Julian day (Julian day - 2450000).
; OPTIONAL INPUT KEYWORDS:
;	/Full    - If present and nonzero, the full Julian day is returned
;	          instead of the MAP Reduced Julian day.
;	/Verbose - If present and nonzero, the result is written
;	          to the screen.  If the date argument is not
;	          supplied, this keyword is assumed.
; EXAMPLE:
;       IDL> DOY2JUL,[222,232],[2002,2002],jul,/verbose
;                2496.50 --  222.00/2002
;                2506.50 --  232.00/2002
;
; PROCEDURES USED:
;       JULDATE
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 29 April 1999.
;       Loop loop variable.  RSH, RITSS, 2 Feb 2001.
;       Fix output display loop  W. Landsman Nov 2002
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, 'Syntax:  DOY2Jul, day, year [, jul]'
;
If ((keyword_set(vrb)) OR (n_params() LT 3)) Then vflg = 1 $
                                             Else vflg = 0
;
If (keyword_set(ful)) Then rdconv = 2400000.0d0 $
                      Else rdconv =  -50000.0d0
;
;			Do the calculation.
;
n = n_elements(day)
If (n_elements(year) LT n) Then Begin
	yr = replicate(year[0], n)
	yr[0] = year[0:*]
EndIf Else Begin
	yr = year
EndElse
;
jul = dblarr(n)

For i = 0L, (n - 1) Do Begin
	juldate, [yr[i], 1, 1], js
	jul[i] = js + day[i] + (rdconv - 1.0d0)
	If (vflg NE 0) Then Begin
		fmt = "(F12.2,' -- ',F7.2,'/',I4)"
		print, FORMAT=fmt, jul[i], day[i], year[i]
	EndIf
EndFor
;
Return
End
