PRO DayOfYear, date, day, Verbose=vrb
;+
; NAME:
;	DayOfYear
; PURPOSE:
;	Determines the day-of-year for a date.
; CALLING SEQUENCE:
;	DayOfYear, date [, day]
; INPUTS:
;	date - A 3-5 element vector containing the date:
;			year, month, day, hour, minute
;	       The hour and minute are set to zero if not supplied.
; OUTPUTS:
;	day  - The real day-of-year.  If not supplied then the result
;	       is written regardless of the VERBOSE keyword.
; INPUT KEYWORDS:
;	/Verbose - If present and nonzero, the result is written
;	          to the screen.
; COMMENTS:
;	The date is converted to Julian date and compared to
;	the Julian date of the start of the year.
; EXAMPLE: 
;       Display the day of year of November 25, 2002
;       IDL> DAYOFYEAR, [2002,11,25],doy,/verb
;            11/25/2002:00:00 --  329.00 
; PROCEDURES USED:
;       JULDATE
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon STX, 15 September 1998.
;	Juldate cannot handle a 6+ element input vector.  This routine
;	   ensures that no more than 5 elements are passed.
;	   MRG, RITSS, 11 May 1999.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  DayOfYear, date [, day]'
;
sz = n_elements(date)
If (sz LT 3) Then message, 'Date argument must have at least 3 elements.'
If (sz GT 5) Then indate = date[0:4] $
             Else indate = date
;
If ((keyword_set(vrb)) OR (n_params() LT 2)) Then vflg = 1 $
                                             Else vflg = 0
;
;			Do the calculation.
;
juldate, indate, jd
juldate, [indate[0], 1, 1], js
day = (jd - js) + 1
If (vflg NE 0) Then Begin
	fmt = "(I2.2,'/',I2.2,'/',I4.4,':',I2.2,':',I2.2,' -- ',F7.2)"
	d = replicate(0, 5) + indate
	print, FORMAT=fmt, d[1], d[2], d[0], d[3], d[4], day
EndIf
;
Return
End
