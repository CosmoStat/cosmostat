Pro TS2GMT, gmt, ts, refgmt
;+
; NAME:
;	TS2GMT
; PURPOSE:
;	Converts a WMAP timestamp into a WMAP GMT.
; CALLING SEQUENCE:
;	TS2GMT, gmt, ts [, refgmt]
; INPUTS:
;	ts     - The timestamp.  This is a two element long array.  The first
;	         element contains the number of days between gmt and refgmt.
;	         The second element contains the fraction of a day in tenths
;	         of a millisecond.
;	refgmt - The reference GMT, defining timestamp [0,0].  If not supplied
;	         the DefTSRef is used.
; OUTPUTS:
;	gmt    - The GMT to convert.  This is a string with the format:
;	            YYYYDDDHHMMSStttt000
;	         where:	YYYY - the year,
;			DDD  - the day of year,
;			HH   - the time of day, hours,
;			MM   - the time of day, minutes,
;			SS   - the time of day, seconds,
;			tttt - the time of day, tenths of a millisecond.
;			000  - padding out to 20 characters.
; COMMENTS:
;	The timestamp is used to determine the difference between the target
;	GMT and the reference GMT, in days.  GMT2JUL is used to convert the
;	reference GMT into a Julian date; the difference is then used to
;	determine the target Julian date.  Finally, JUL2GMT is used to find
;	the target GMT.
; EXAMPLE:
;      IDL> TS2GMT, gmt, [513, 677980204] & print,gmt
;         ===> '20023291849580204490'
; PROCEDURES USED:
;       DefTSRef(), GMT2JUL(), JUL2GMT()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 10 August 1998
;	DefTSRef returns the default reference date.  MRG, RITSS, 
;	  12 January 1999.
;-
on_error, 2
;
;			Check arguments.
;
 If (n_params() LT 2) Then begin 
       print, 'Syntax:  TS2GMT, gmt, ts [, refgmt]'
       return
 endif
;
If (n_params() GE 3) Then rgmt = refgmt $
                     Else rgmt = DefTSRef()
;
;			Compute the difference in Julian date from the
;			timestamp.
;
djul = double(ts[0]) + ((double(ts[1]) + 0.49d0) / 8.64d+08)
;
;			Convert this difference into a GMT using the
;			reference GMT.
;
jul = gmt2jul(rgmt) + djul
gmt = jul2gmt(jul)
;
RETURN
END
