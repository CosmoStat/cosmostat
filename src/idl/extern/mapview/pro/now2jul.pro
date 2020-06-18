Function Now2Jul, dum
;+
; NAME:
;	Now2Jul()
; PURPOSE:
;	Returns the current system time as a MAP Reduced Julian day.
; CALLING SEQUENCE:
;	jul = Now2Jul()
; RETURNED:
;	jul - The MAP Reduced Julian day corresponding to the current
;	      system time.
; COMMENTS:
;	"systime" is used to get the system time; this string is then
;	parsed into its components, which are then passed to Date2Jul.
;       In IDL V5.4 or later one can simply use 
;       systime(/JULIAN,UTC) - 2450000.0d 
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 28 February 2000.
;-
on_error, 2
;
if !VERSION.RELEASE GE '5.4' then $
     return,systime(/JULIAN,/UTC)-2450000.0d
     
mon = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', $
       'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
;
;			Get the current time and parse it.
;
str = strlowcase(strtrim(systime(), 2))
;
n   = strlen(str)
p   = strpos(str, ' ') + 1
str = strtrim(strmid(str, p, n), 2)
;
;				Month.
;
p     = strpos(str, ' ')
tmp   = strtrim(strmid(str, 0, p), 2)
w     = where(mon eq tmp)
If (w[0] LT 0) Then message, 'No month'
mon_n = w[0] + 1
;
;				Day.
;
str   = strtrim(strmid(str, (p + 1), n), 2)
p     = strpos(str, ' ')
tmp   = strtrim(strmid(str, 0, p), 2)
day_n = long(tmp)
;
;				Hour.
;
str  = strtrim(strmid(str, (p + 1), n), 2)
p    = strpos(str, ':')
tmp  = strtrim(strmid(str, 0, p), 2)
hr_n = long(tmp)
;
;				Minute.
;
str   = strtrim(strmid(str, (p + 1), n), 2)
p     = strpos(str, ':')
tmp   = strtrim(strmid(str, 0, p), 2)
min_n = long(tmp)
;
;				Second.
;
str   = strtrim(strmid(str, (p + 1), n), 2)
p     = strpos(str, ' ')
tmp   = strtrim(strmid(str, 0, p), 2)
sec_n = long(tmp)
;
;				Year.
;
str    = strtrim(strmid(str, (p + 1), n), 2)
year_n = long(str)
;
;			Compute and return the MAP Reduced Julian day.
;
Return, Date2Jul(year_n, mon_n, day_n, hr_n, min_n, sec_n)
End
