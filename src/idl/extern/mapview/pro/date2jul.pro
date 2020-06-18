Function Date2Jul, yr, mon, day, hr, mn, sc
;+
; NAME:
;	Date2Jul
; PURPOSE:
;	Converts a Gregorian date and time into a reduced Julian date.
;	Unlike a classic reducted Julian date, the full Julian date is
;	recovered by adding 2450000.
; CALLING SEQUENCE:
;	jul = Date2Jul(yr, mon, day, hr [, mn, sc])
; INPUTS:
;	yr  - The Gregorian year (must be 4 digits!).  Integer.
;	mon - The month (1-12).  Integer.
;	day - The day of month (1-31).  Integer.
;	hr  - The time of day in hours.
;	mn  - The time of day in minutes.
;	sc  - The time of day in seconds.
; RETURNED:
;	jul - The reduced Julian day.
; COMMENTS:
;	This routine is based upon the Fortran routine in the MAP library.
;
;	If mn is supplied then hr is treated as an integer.  If sc is supplied
;	then mon is treated as an integer.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon STX, 07 August 1998.
;	mn, sc added.  MRG, RITSS, 18 October 1999.
;       accepts arrays of dates, JW, 09 November 1999.
;	Negative value "roundoff" corrected.  MRG, RITSS, 29 September 2000.
;-
on_error, 2
;
If (n_params() LT 4) Then message, 'Syntax:  jul = Date2Jul(yr, mon, day, hr [, mn, sc])'
;
Case (n_params()) Of
     4 : hour = double(hr)
     5 : hour = double(fix(hr)) + (double(mn) / 60.d0)
  Else : hour = double(fix(hr)) + (double(fix(mn)) / 60.d0) + (double(sc) / 3600.d0)
EndCase
;
m = long(mon)

mle2 = where(m le 2)
mgt2 = where(m gt 2)

y = lonarr(n_elements(yr))
if (mle2[0] ne -1) then begin
  y[mle2] = long(yr[mle2] - 1L)
  m[mle2] = m[mle2] + 12L
endif
if (mgt2[0] ne -1) then begin
  y[mgt2] = long(yr[mgt2])
endif

d = double(day) + (hour / 24.0D0)
a = floor(y / 100)
b = 2 - a + floor(a / 4)
jul = double(floor(365.25 * (y - 1992)) + long(30.6001 * (m + 1)) + d + b - 1427.5)
;

if (n_elements(jul) eq 1) then jul = jul[0]

Return, jul
end
