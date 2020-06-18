Pro TimeTransform, input, output,                                       $
        Date2Jul=d2j,  DayOfYear=doy, DispGMT=dspg,  DispTS=dspt,       $
        DOY2Date=dy2d, DOY2Jul=dy2j,  GMT2Jul=g2j,   GMT2TS=g2t,        $
        GMT2YMD=g2ymd, Jul2Date=j2d,  Jul2GMT=j2g,   Jul2Tel=j2t,       $
        Tel2Jul=t2j,   TS2GMT=t2g,    TS2Jul=ts2j,   Reference=ref,     $
        Verbose=vrb, _EXTRA=ext
;+
; NAME:
;       TimeTransform
; PURPOSE:
;       Converts one time format into another.
; CALLING SEQUENCE:
;       TimeTransform, input [, output], TransKeyword
; INPUTS:
;       input  - The input time.  (See KEYWORDS for details.)
; OUTPUTS:
;       output - The converted time.  (See KEYWORDS for details.)
; INPUT KEYWORDS:
;       /Date2Jul  - Converts a Gregorian date into reduced Julian day.
;                   Input should be a 3-6 element array containing the date:
;                   [year, month, day, hour, minute, second].  Time of day
;                   elements are optional.
;       /DayOfYear - Returns the day of the year for a given date.  Input
;                   should be a 3-5 element array containing the date:
;                   [year, month, day, hour, minute]
;       /DispGMT   - Convert a GMT string into a more readable string.
;       /DispTS    - Convert a WMAP timestamp into a more human readable format.
;       /DOY2Date  - Determines the date from day-of-year and year.  Input
;                   should be a 2 element array containing the day of year
;                   and year.  Output will be a 5 element array containing
;                   the date: [year, month, day, hour, minute]
;       /DOY2Jul   - Determines the reduced Julian day from day-of-year and
;                   year.  Input should be a 2 element array containing the 
;                   day of year and year.
;       /GMT2Jul   - Converts a WMAP GMT date/time string into a Julian day.
;       /GMT2TS    - Converts a WMAP GMT into a MAP timestamp.
;       /GMT2YMD   - Converts the YYYYDDD portion of a GMT string to the format
;                   YYYY:MM:DD.
;       /Jul2Date  - Converts a reduced Julian day into a Gregorian date
;                   and time.  Output consists of a 6 element array containing
;                   the date: [year, month, day, hour, minute, second].
;       /Jul2GMT   - Converts a WMAP GMT date/time string into a reduced
;                   Julian day.
;       /Jul2Tel   - Converts a reduced Julian day into a telemetry structure
;                   timestamp.  The telemetry structure must already exist!
;       /Tel2Jul   - Converts a packet time stamp into a reduced Julian day.
;       /TS2GMT    - Converts a WMAP timestamp into a WMAP GMT.
;       /TS2Jul    - Converts an array of time stamps into reduced Julian days.
;
;       Reference - The timestamp reference time in GMT format.
;       /Verbose   - If present and nonzero, the output value is written
;                   to the screen.
;       _EXTRA    - IDL keyword inheritance.  Keywords required by the
;                   conversion routines can be passed directly to them
;                   simply by specifying them in the call to this routine.
; COMMENTS:
;       A conversion keyword indicating the type of conversion MUST be
;       specified.
;
;       Reduced Julian days referred to here are WMAP Reduced Julian days:
;         Full Julian day - 2450000
;
;       See the individual conversion routines for keyword information.
;       Reference GMTs are always passed through the Reference keyword
;       regardless of the _EXTRA inheritance mechanism.
; MODIFICATION HISTORY:
;       Written by Michael R. Greason, Raytheon ITSS, 11 May 1999.
;       Adapted to new Date2Jul and Jul2Date.  MRG, RITSS, 18 October 1999.
;-
on_error, 2
;
;                       Check arguments.
;
np = n_params()
If (np LT 1) Then message, 'Syntax: TimeTransform, input [, output]'
;
;                       Determine the type of transformation.
;
type = -1                               ; No default.
If ((type LT 0) AND (keyword_set(g2j)))   Then type =  0
If ((type LT 0) AND (keyword_set(j2g)))   Then type =  1
If ((type LT 0) AND (keyword_set(g2t)))   Then type =  2
If ((type LT 0) AND (keyword_set(t2g)))   Then type =  3
If ((type LT 0) AND (keyword_set(g2ymd))) Then type =  4
;
If ((type LT 0) AND (keyword_set(d2j)))   Then type = 10
If ((type LT 0) AND (keyword_set(j2d)))   Then type = 11
;
If ((type LT 0) AND (keyword_set(j2t)))   Then type = 20
If ((type LT 0) AND (keyword_set(t2j)))   Then type = 21
If ((type LT 0) AND (keyword_set(ts2j)))  Then type = 22
;
If ((type LT 0) AND (keyword_set(doy)))   Then type = 30
If ((type LT 0) AND (keyword_set(dy2d)))  Then type = 31
If ((type LT 0) AND (keyword_set(dy2j)))  Then type = 32
;
If ((type LT 0) AND (keyword_set(dspg)))  Then type = 90
If ((type LT 0) AND (keyword_set(dspt)))  Then type = 91
;
;                       Perform the transformation.
;
Case (type) Of
;
 0 : output = GMT2Jul(input)
     
 1 : output = Jul2GMT(input)
    
 2 : If (n_elements(ref) LE 0) Then GMT2TS, output, input $
                                  Else GMT2TS, output, input, ref
    
 3 : If (n_elements(ref) LE 0) Then TS2GMT, output, input $
                               Else TS2GMT, output, input, ref
     
 4 : GMT2YMD, input, output

10 : Begin
        n = n_elements(input)
        Case (1) Of
           (n LE 3) : output = Date2Jul(input[0], input[1], input[2], 0.0d0)
           (n EQ 4) : output = Date2Jul(input[0], input[1], input[2], input[3])
           (n EQ 5) : output = Date2Jul(input[0], input[1], input[2], input[3], input[4])
           Else     : output = Date2Jul(input[0], input[1], input[2], input[3], input[4], input[5])
        EndCase
     End
11 : Begin
        Jul2Date, input, yr, mon, day, hr, mn, sc
        output = [yr, mon, day, hr, mn, fix(sc)]
     End
;
20 : Jul2Tel, input, output, REFGMT=ref
    
21 : output = Tel2Jul(input, REFGMT=ref)
    
22 : output = TS2Jul(input, REFGMT=ref)
     
30 : DayOfYear, input, output
     
31 : DOY2Date, input[0], input[1], output, _EXTRA=ext
     
32 : DOY2Jul, input[0], input[1], output, _EXTRA=ext
    

90 : output = DispGMT(input, _EXTRA=ext)
     
91 : If N_elements(ref) LE 0 Then output = DispTS(input, _EXTRA=ext) $
                             Else output = DispTS(input, ref, _EXTRA=ext)
    
Else : message, 'Invalid transformation requested!'
EndCase
;
;                       Report it.
;
If ((keyword_set(vrb)) OR (np LT 2)) Then print, output
;
Return
End
