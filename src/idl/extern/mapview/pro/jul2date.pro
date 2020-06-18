Pro Jul2Date, jul, yr, mon, day, hr, mn, sc
;+
; NAME:
;	Jul2Date
; PURPOSE:
;	Converts a reduced Julian date into a Gregorian date and time.
;	Unlike a classic reduced Julian date, the full Julian date is
;	recovered by adding 2450000.
; CALLING SEQUENCE:
;	Jul2Date, jul, yr, mon, day, hr,[ mn, sc ]
; INPUTS:
;	jul - The reduced Julian day, scalar
; OUTPUTS:
;	yr  - The Gregorian year (must be 4 digits!).  Integer.
;	mon - The month (1-12).  Integer.
;	day - The day of month (1-31).  Integer.
;	hr  - The time of day in hours.
;	mn  - The time of day in minutes.
;	sc  - The time of day in seconds.
; EXAMPLE:
;       IDL> jul2date,2405.6,yr,mon,day,hr,mn,sc & print,yr,mon,day,hr,mn,sc
;          2002           5          11       2      24       8.4375000
;            
; COMMENTS:
;	This routine is based upon the Fortran routine in the MAP library.
;
;	If mn is supplied then hr is treated as an integer.  If sc is supplied
;	then mon is treated as an integer.
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 07 August 1998
;	mn, sc added.  MRG, RITSS, 18 October 1999.
;-
on_error, 2
;
IF (n_params() LT 5) Then begin 
    print,'Syntax:  Jul2Date, jul, yr, mon, day, hr, [mn, sc]'
    return
 endif   
;
;			Compute components.
;
z = double(long(jul + 0.5))
f = double((jul + 0.5) - z)
z = double(z + 2450000L)
If (z LT 2299161L) Then Begin
	a = z
Endif Else Begin
	a = double(long((z - 1867216.25d0) / 36524.25d0))
	a = z + 1 + a - long(a / 4)
Endelse
b = a + 1524
c = double(long((b - 122.1) / 365.25))
d = double(long(365.25 * c))
e = double(long((b - d) / 30.6001))
dfr = b - d - long(30.6001 * e) + f
;
;			Extract date.
;
If (e LT 14) Then Begin
	mon = long(e - 1)
Endif Else Begin
	mon = long(e - 13)
Endelse
If (mon GT 2) Then Begin
	yr = long(c - 4716)
Endif Else Begin
	yr = long(c - 4715)
Endelse
day = long(dfr)
hr = (dfr * 24.0d0) MOD 24.0d0
;
;			Break the hour into hour:min:sec.
;
Case (n_params()) Of
     7 : Begin
	 mn = (hr * 60.0d0) mod 60.0d0
	 hr = fix(hr)
	 sc = (mn * 60.0d0) mod 60.0d0
	 mn = fix(mn)
	 End
     6 : Begin
	 mn = (hr * 60.0d0) mod 60.0d0
	 hr = fix(hr)
	 End
  Else : dfr = hr
EndCase
;
Return
End
