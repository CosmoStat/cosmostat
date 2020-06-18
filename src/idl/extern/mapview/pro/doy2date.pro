Pro DOY2Date, day, year, date, Verbose=vrb
;+
; NAME:
;	DOY2Date
; PURPOSE:
;	Determines the date from day-of-year and year.
; CALLING SEQUENCE:
;	DOY2Date, day, year [, date, /VERBOSE]
; INPUTS:
;	day  - The day-of-year, scalar integer (1-366)
;	year - The year, scalar integer
; OUTPUTS:
;	date - A 5 element vector containing the date:
;			year, month, day, hour, minute
; OPTIONAL INPUT KEYWORD:
;	/Verbose - If present and nonzero, the result is written
;	          to the screen.  If the date argument is not
;	          supplied, this keyword is assumed.
; EXAMPLE:
;      Print the date corresponding to day 221 of the year 2002 
;      IDL> doy2date,221,2002,date,/verbose
;           ===>08/09/2002:00:00 --  221.00
; COMMENTS:
;	The day and year are converted to Julian date, which is then
;	converted into the date vector.
; PROCEDURES USED:
;       DAYCNV, JULDATE
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 27 April 1999.
;-
 On_error, 2
;
;			Check arguments.
;
 If (N_params() LT 2) Then begin 
    print, 'Syntax:  DOY2Date, day, year [, date,/VERBOSE]'
    return
 EndIf
;
If ((keyword_set(vrb)) OR (n_params() LT 3)) Then vflg = 1 $
                                             Else vflg = 0
;
;			Do the calculation.
;
juldate, [year, 1, 1], js
jd = js + day - 1 + 2400000.d0
daycnv, jd, yr, mn, dy, hr
date = [fix(yr), fix(mn), fix(dy), fix(hr), fix((hr - fix(hr)) * 60.)]
If (vflg NE 0) Then Begin
	fmt = "(I2.2,'/',I2.2,'/',I4.4,':',I2.2,':',I2.2,' -- ',F7.2)"
	print, FORMAT=fmt, date[1], date[2], date[0], date[3], date[4], day
EndIf
;
Return
End
