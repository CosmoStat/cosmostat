Function GMT2Jul, gmt, Verbose=vrb
;+
; NAME:
;       GMT2Jul
; PURPOSE:
;       Converts a MAP GMT date/time string into a Julian date.
; CALLING SEQUENCE:
;       jul = GMT2Jul(gmt, [/Verbose] )
; INPUTS:
;       gmt - The GMT:  YYYYDDDHHMMSStttt000
;             where:    YYYY - 4 digit year,
;                       DDD  - 3 digit day of year,
;                       HH   - 2 digit hour of day,
;                       MM   - 2 digit minute of hour,
;                       SS   - 2 digit second of minute,
;                       tttt - 7 digit fraction of a second.
; RETURNED:
;       jul - The reduced Julian date--the Julian date with 2450000.D0
;             subtracted off, double precision
; INPUT KEYWORDS:
;       /Verbose - Problems are reported (and bandaged over).
; COMMENTS:
;       This routine extracts the Gregorian date from the GMT string
;       and calls Date2Jul to convert it to reduced Julian date.
; PROCEDURES CALLED:
;       Date2Jul()                         ;Map Library
; COMMON BLOCKS:
;       GMT2JUL_PERSIST_TABS_1 - stores the day/month information 
;          Should not be coded into any other routine.
; MODIFICATION HISTORY:
;       Written by Michael R. Greason, Raytheon STX, 07 August 1998.
;       Support for an array of GMTs.  If the GMT contains day 0, day 1
;          is used.  MRG, RITSS, 24 July 2000.
;       Fraction of a second expanded to 7 digits. MRG, RITSS, 04 October 2000.
;   Day of year lookup tables, vectorization.  RSH, RITSS, 3 Apr 2001.
;   Pads short GMT strings with zeros.  MRG, JP & RSH, RITSS, 6 April 2001
;-
 On_error, 2

;
;   This common block is here purely to make this set of lookup
;   tables persist from one call to the next.  DO NOT code it
;   into any other routine.
;
common gmt2jul_persist_tabs_1, monthtab, daytab, monthtab_ly, daytab_ly

If n_elements(monthtab) NE 365 $
    OR n_elements(daytab) NE 365 $
    OR n_elements(monthtab_ly) NE 366 $
    OR n_elements(daytab_ly) NE 366 Then Begin

    monthtab = bytarr(365)
    daytab = bytarr(365)
    monthtab_ly = bytarr(366)
    daytab_ly = bytarr(366)

    dom = [31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L]

    tot = 0 
    For i=0,11 Do Begin
        monthtab[tot:tot+dom[i]-1] = i + 1
        daytab[tot] = bindgen(dom[i]) + 1
        tot = tot + dom[i]
    Endfor

    dom[1] = 29
    tot = 0 
    For i=0,11 Do Begin
        monthtab_ly[tot:tot+dom[i]-1] = i + 1
        daytab_ly[tot] = bindgen(dom[i]) + 1
        tot = tot + dom[i]
    Endfor
Endif

;
;                       Check arguments.
;
 If (n_params() LT 1) Then message, 'Syntax:  jul = GMT2Jul(gmt)'
;
 n = n_elements(gmt)
 If (n LE 0) Then message, 'At least one GMT must be supplied.'

 jul = dblarr(n)

 gmtpad = strmid(gmt+'00000000000000000000', 0, 20)
;
;
;                       Extract gregorian date.
;
;                               Year.
;
 yr = long(strmid(gmtpad, 0, 4))
;
;                               Month and day.
;
 doy = long(strmid(gmtpad, 4, 3)) - 1
 day = doy
 mn = doy

 leap = where(yr EQ 4*floor(yr/4), nleap, comple=nonleap, ncomple=nnon)

 If nleap GT 0 Then Begin
    d = (doy[leap]>0)<365
    day[leap] = daytab_ly[d]
    mn[leap] = monthtab_ly[d]
 Endif

 If nnon GT 0 Then Begin
    d = (doy[nonleap]>0)<364
    day[nonleap] = daytab[d]
    mn[nonleap] = monthtab[d]
 Endif

;    
;                               Time, in hours.
;
 t1 = long(strmid(gmtpad,7,6))
 t2 = long(strmid(gmtpad,13,7))

 h = t1/10000       & hj = 10000*h
 m = (t1 - hj)/100  & mj = 100*m
 s = t1 - hj - mj

 hr = double(h) + m/60.0d0 + s/3600.0d0 + t2/3600.0d7

;
;                       Convert to julian date.
;
 jul = Date2Jul(yr, mn, day, hr)

 If (n LE 0) Then jul = jul[0]

 Return, jul
 End
