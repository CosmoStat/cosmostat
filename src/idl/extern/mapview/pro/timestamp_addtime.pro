PRO TimeStamp_AddTime, OutTS, InTS, Dtime, UNITS=units
;+
; NAME:
;	TimeStamp_AddTime
; PURPOSE:
;	Adds to a MAP Omega time stamp.
; CALLING SEQUENCE:
;	TimeStamp_AddTime, OutTS, InTS, Dtime
; INPUTS:
;	InTS  - The input time stamp.  This is a two element array containing
;	        the number of days in the first element and the number of
;	        tenths of milliseconds into the day in the second element.
;	Dtime - The amount of time to add.  See the UNITS description for
;	        the units of this argument.
; OUTPUTS:
;	OutTS - The output time stamp.  This is a two element array containing
;	        the number of days in the first element and the number of
;	        tenths of milliseconds into the day in the second element.
; KEYWORDS:
;	UNITS - Specifies the units of Dtime:
;		  0 - tenths of milliseconds.  Scalar.  Default.
;		  1 - seconds.  Scalar.
;		  2 - days.  Scalar.
;		  3 - another time stamp.  Two element array.
; COMMENTS:
;	A CASE construct performs the arithmetic for the various units.  The
;	output time stamp's time-of-day element is forced to fall within
;	limits, with the day element adjusted as necessary to keep time
;	constant.
; EXAMPLE:
;       Add 863000000 milliseconds to a time stamp of 34 days and 
;           1111322 millliseconds
; 
;       IDL> timestamp_addtime,out,[34,1111322],863000000
;            ===> out = [ 35, 111322 ]
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 22 July 1998
;-
on_error, 2
;
;			Check arguments.
;
IF (n_params() LT 3) THEN $
	message, 'Syntax:  TimeStamp_AddTime, OutTS, InTS, Dtime'
;
IF (n_elements(units) LE 0) THEN units = 0
;
;			Constants.
;
Day2Sec = 86400L
Sec2TMS = 10000L
Day2TMS = Day2Sec * Sec2TMS
;
;			Perform the arithmetic.
;
OutTS = long(InTS)      ;Make sure at least long array
CASE (units) OF
;
;				Add another timestamp.
;
    3 :	BEGIN
	OutTS = OutTS + Dtime
	END
;
;				Add days.
;
    2 :	BEGIN
	OutTS[1] = OutTS[1] + long(Dtime * Day2TMS)
	END
;
;				Add seconds.
;
    1 :	BEGIN
	OutTS[1] = OutTS[1] + long(Dtime * Sec2TMS)
	END
;
;				Add tenths of milliseconds.
;
 ELSE :	BEGIN
	OutTS[1] = OutTS[1] + long(Dtime)
	END
;
ENDCASE
;
;			Force the output TMS component to fall within limits.
;
WHILE (OutTS[1] GE Day2TMS) DO BEGIN
	OutTS[1] = OutTS[1] - Day2TMS
	OutTS[0] = OutTS[0] + 1L
ENDWHILE
;
WHILE (OutTS[1] LT 0L) DO BEGIN
	OutTS[1] = OutTS[1] + Day2TMS
	OutTS[0] = OutTS[0] - 1L
ENDWHILE
;
RETURN
END
