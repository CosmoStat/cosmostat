FUNCTION timestamp_diff, ts1, ts2, DAYS=days, TS=tsk
;+
; NAME:
;	TIMESTAMP_DIFF
; PURPOSE:
;	To determine the diffence between two time stamps.
; CALLING SEQUENCE:
;	diff = timestamp_diff(ts1, ts2)
; INPUTS:
;	ts1  - A two element integer array giving the first time stamp,
;	       with the day in the first element, and the time of day in
;	       tenths of milliseconds in the second.
;	ts2  - A two element integer array giving the second time stamp,
;	       with the day in the first element, and the time of day in
;	       tenths of milliseconds in the second.
; RETURNED:
;	diff - The difference between the two, in tenths of milliseconds.
;	       Double precision.
; OPTIONAL INPUT KEYWORDS:
;	/DAYS - If present and nonzero, the difference is returned in days
;	       instead of tms.
;	/TS   - If present and nonzero, the difference is returned as a
;	       timestamp array.  This keyword takes precedence over DAYS.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Hughes STX, 19 September 1997.
;	TS keyword added.  MRG, RSTX, 10 August 1998.
;-
on_error, 2
;
;			Check arguments.
;
IF (n_params() LT 2) THEN message, 'Syntax:  diff = timestamp_diff(ts1, ts2)'
;
;			Constants.
;
d2tmsL = 864000000L		; Days to tenths-of-milliseconds.
d2tms  = double(d2tmsL)
;
;			Determine the difference.
;
IF (keyword_set(tsk)) THEN BEGIN
	diff = long(ts1) - long(ts2)
	WHILE (diff[1] LT 0) DO BEGIN
		diff[1] = diff[1] + d2tmsL
		diff[0] = diff[0] - 1L
	ENDWHILE
	WHILE (diff[1] GE d2tmsL) DO BEGIN
		diff[1] = diff[1] - d2tmsL
		diff[0] = diff[0] + 1L
	ENDWHILE
	RETURN, diff[0:1]
ENDIF
;
IF (keyword_set(days)) THEN BEGIN
	diff = double(ts1[0] - ts2[0]) + (double(ts1[1] - ts2[1]) / d2tms)
ENDIF ELSE BEGIN
	diff = (double(ts1[0] - ts2[0]) * d2tms) + double(ts1[1] - ts2[1])
ENDELSE
;
RETURN, diff
END
