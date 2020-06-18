Pro GMT2TS, ts, gmt, refgmt
;+
; NAME:
;	GMT2TS
; PURPOSE:
;	Converts a MAP GMT into a MAP timestamp.
; CALLING SEQUENCE:
;	GMT2TS, ts, gmt [, refgmt]
; INPUT:
;	gmt    - The GMT to convert.  This is a string with the format:
;	            YYYYDDDHHMMSStttt000
;	         where:	YYYY - the year,
;			DDD  - the day of year,
;			HH   - the time of day, hours,
;			MM   - the time of day, minutes,
;			SS   - the time of day, seconds,
;			tttt - the time of day, tenths of a millisecond.
;			000  - padding out to 20 characters.
; OPTIONAL INPUT:
;	refgmt - The reference GMT, defining timestamp [0,0].  If not supplied
;	         the DefTSRef is used.
; OUTPUTS:
;	ts     - The timestamp.  This is a two element long array.  The first
;	         element contains the number of days between gmt and refgmt.
;	         The second element contains the fraction of a day in tenths
;	         of a millisecond.
; EXAMPLE:
;       IDL> GMT2TS,ts,'20023291849580204716' & print,ts
;            513   677980204
; COMMENTS:
;	The two GMT's are converted to reduced Julian date using GMT2JUL.
;	The difference between the two Julian dates are used to compute
;	the timestamp.
; PROCEDURES USED:
;       DefTSRef(), GMT2Jul() 
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 10 August 1998
;	DefTSRef returns the default reference date.  MRG, RITSS, 
;	  12 January 1999.
;-
 On_error, 2
;
;			Check arguments.
;
 If (N_params() LT 2) Then begin 
     print,'Syntax:  GMT2TS, ts, gmt [, refgmt]'
     return
 endif
;
If (n_params() GE 3) Then rgmt = refgmt $
                     Else rgmt = DefTSRef()
;
;			Compute the difference in Julian date between
;			the input GMT and reference GMT.
;
jcnv = gmt2jul(gmt)
jref = gmt2jul(rgmt)
djul = jcnv - jref
;
;			Convert this difference into a timestamp.
;
ts = lonarr(2)
ts[0] = long(djul)
ts[1] = long((djul - ts[0]) * 8.64d+08)
;
RETURN
END
