Function Jul2GMT, jul
;+
; NAME:
;	Jul2GMT
; PURPOSE:
;	Converts  Julian date(s) into a MAP GMT date/time string(s).
; CALLING SEQUENCE:
;	gmt = Jul2GMT(jul)
; INPUTS:
;	jul - The reduced Julian date--the Julian date with 2450000.D0
;	      subtracted off, scalar or vector
; RETURNED:
;	gmt - The GMT:  YYYYDDDHHMMSStttt000
;	      where:	YYYY - 4 digit year,
;			DDD  - 3 digit day of year,
;			HH   - 2 digit hour of day,
;			MM   - 2 digit minute of hour,
;			SS   - 2 digit second of minute,
;			tttt - 7 digit fraction of a second.
; COMMENTS:
;	This routine extracts the Gregorian date from the Julian date
;	using Jul2Date and then builds the GMT string.
; PROCEDURES USED:
;       JUL2DATE
; EXAMPLE:
;       Express the reduced Julian dates 2203.0 and 2604.2025 in GMT
;       IDL> print,jul2gmt( [ [2203.0d, 2604.2025d])
; 
;            ==> 20012931200000000000 20023291651359999999
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon STX, 07 August 1998.
;	Support for an array of GMTs.  MRG, RITSS, 24 July 2000.
;	Fraction of a second expanded to 7 digits. MRG, RITSS, 04 October 2000.
;	Documentation correction.  MRG, RITSS, 30 March 2001.
;-
 On_error, 2
;
;			Check arguments.
;
 If (N_params() LT 1) Then begin 
       print, 'Syntax:  gmt = Jul2GMT(jul)'
       return,''
 endif
;
 n = n_elements(jul)
 If (n LE 0) Then message, 'You must supply a Julian day.'
 gmt = strarr(n)
;
 n = n - 1
 For j = 0L, n Do Begin
    jele = jul[j]
;
;			Extract gregorian date.
;
    Jul2Date, jele, yr, mn, day, hr
    minutes = (hr * 60.D0) MOD 60.D0
    seconds = (minutes * 60.D0) MOD 60.D0
    tms     = long(seconds * 10000000.D0) MOD 10000000L
;
    dom = [31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L]
    If ((yr MOD 4) EQ 0) Then dom[1] = 29L
    doy = replicate(0L, 12)
    For i = 0, 10 Do doy[i+1] = doy[i] + dom[i]
    day = day + doy[long(mn) - 1]
;
;			Convert to GMT string.
;
    gmt[j] = string(Format="(I4.4,I3.3,I2.2,I2.2,I2.2,I7.7)", $
	    long(yr), long(day), long(hr), long(minutes), long(seconds), tms)
EndFor
;
 If (n LE 0) Then gmt = gmt[0]
;
 RETURN, gmt
 END
