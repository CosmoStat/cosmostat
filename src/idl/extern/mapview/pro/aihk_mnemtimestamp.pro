Pro AIHK_MnemTimeStamp, aihk, mnem, ts, Julian=jul, RefGMT=refgmt, $
	Status=status
;+
; NAME:
;	AIHK_MnemTimeStamp
; PURPOSE:
;	Returns an array of time stamps for a given mnemonic in
;	an array of analog instrument housekeeping (AIHK) elements.
; CALLING SEQUENCE:
;	AIHK_MnemTimeStamp, aihk, mnem, ts
; INPUTS:
;	aihk - The array of AIHK elements.  At least one must be supplied.
;	mnem - The mnemonic to search for.
; OUTPUTS:
;	ts   - A 4xN array of time stamps.  N is the number of aihk elements.
;	       For each aihk element/row, the first two elements/columns
;	       is the time stamp for the first sweep and the second two
;	       elements/columns is the time stamp for the second sweep.
;	       A time stamp is a two element array where the first element
;	       is the number of days since the reference date and the
;	       second element is the number of tenths of a millisecond in
;	       that day.
; OPTIONAL INPUT KEYWORDS:
;	/Julian - If present and nonzero, ts is a 2xN array of MAP Reduced
;	         Julian days.
;	RefGMT - The reference GMT, defining timestamp [0,0].  If not
;	         supplied the DefTSRef is used.
; OPTIONAL OUTPUT KEYWORDS:
;	Status - A status code: 1=success, 0=not found.
; EXAMPLE:
;      IDL> fits_read_tod,'MAP_tod_20022162357_20022172357.fits',tod
;      IDL> AIHK_MnemTimeStamp,tod.aeu,'DFW222B5DNI',ts,/julian
;
;      The output array ts will be dblarr(2,1875) containing the MAP reduced
;      Julian date for the DFW222B5DNI' mnemonic for both sweeps of the 1875 
;      AIHK elements
; COMMENTS:
;	AIHK_MnemTime is used to get the time offsets into each element.
;	TimeStamp_AddTime is used to apply the offsets to each element.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 21 September 1999.
;	Julian day support.  MRG, RITSS, 09 January 2001.
;       Set loop limits to LONG.  JP, 01 Apr 2002.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 3) Then message, 'Syntax: AIHK_MnemTimeStamp, aihk, mnem, ts'
;
n = n_elements(aihk)
If (n LE 0) Then message, 'At least one AIHK element must be supplied.'
;
;			Get the time offsets for the mnemonic.
;
AIHK_MnemTime, mnem, tms1, tms2, Status=status
If (status EQ 0) Then Return
;
;			Apply the offsets to the packet times.
;
ts = replicate(0L, 4, n)
n = n - 1L
For i = 0L, n Do Begin
	tsin  = [aihk[i].Day, aihk[i].Time]
;
	TimeStamp_AddTime, tsout, tsin, tms1, Units=0
	ts[0,i] = tsout[0:1]
;
	TimeStamp_AddTime, tsout, tsin, tms2, Units=0
	ts[2,i] = tsout[0:1]
EndFor
;
;			Convert to Julian day.
;
If (keyword_set(jul)) Then Begin
	j = transpose([[ts2jul(ts[0:1,*])], [ts2jul(ts[2:3,*])]])
	ts = j
EndIf
;
Return
End
