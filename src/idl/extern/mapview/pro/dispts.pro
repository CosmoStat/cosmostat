Function DispTS, ts, refgmt, Date=dt, Local=loc
;+
; NAME:
;	DispTS
; PURPOSE:
;	Converts a MAP timestamp into a more human readable format.
; CALLING SEQUENCE:
;	disp = DispTS(ts [, refgmt, /Date, LOCAL= ])
; INPUTS:
;	ts     - The timestamp.  This is a two element long array.  The first
;	         element contains the number of days between gmt and refgmt.
;	         The second element contains the fraction of a day in tenths
;	         of a millisecond.
; OPTIONAL INPUT:
;	refgmt - The reference GMT, defining timestamp [0,0].  If not supplied
;	         the DefTSRef is used.
; RETURNED:
;	disp  - The string containing the converted date, in GMT epoch:
;		    dd-mon-year hh:mm:ss
; KEYWORDS:
;	/Date  - If present and nonzero, only the date section is returned.
;	Local - If present and nonzero, the time is returned as local
;	        time instead of GMT.  The value of local indicates the
;	        timezone to display:
;			 1 = Eastern Standard Time		(EST)
;			-1 = Eastern Daylight Savings Time	(EDT)
;			 2 = Central Standard Time		(CST)
;			-2 = Central Daylight Savings Time	(CDT)
;			 3 = Mountain Standard Time		(MST)
;			-3 = Mountain Daylight Savings Time	(MDT)
;			 4 = Pacific Standard Time		(PST)
;			-4 = Pacific Daylight Savings Time	(PDT)
;	        The timezone abbreviation will be appended to the output string.
; EXAMPLE:
;       IDL> print,dispts([513,677980204])
;            ==> 25-Nov-2002 18:49:58 GMT
; COMMENTS:
;	The timestamp is converted into GMT using TS2GMT; this, in turn,
;	is converted into the display string using DispGMT.
; PROCEDURES USED:
;       DefTSRef(), DispGMT, TS2GMT
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 12 January 1999.
;-
 On_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  disp = DispTS(ts [, refgmt])'
;
If (n_params() GE 2) Then rgmt = refgmt $
                     Else rgmt = DefTSRef()
;
;			Go for it.
;
ts2gmt, gmt, ts, rgmt
disp = DispGMT(gmt, Date=dt, Local=loc)
;
Return, disp
End
