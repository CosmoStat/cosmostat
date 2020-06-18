Pro AIHK_MnemTime, mnem, tms1, tms2, Status=status
;+
; NAME:
;	AIHK_MnemTime
; PURPOSE:
;	Returns the time offset into a AIHK packet for the measurement
;	of a given mnemonic.
; CALLING SEQUENCE:
;	AIHK_MnemTime, mnem, tms1 [, tms2]
; INPUTS:
;	mnem - The mnemonic to search for.
; OUTPUTS:
;	tms1 - The time offset in tenths of a millisecond for the mnemonic
;	       in the first sweep.  Required.
;	tms2 - The time offset in tenths of a millisecond for the mnemonic
;	       in the second sweep.  Optional.
; OPTIONAL OUTPUT KEYWORD:
;	Status - A status code: 1=success, 0=not found.
; COMMENTS:
;	The offset into the packet data arrays for the mnemonic is used
;	to compute the time offset.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 21 September 1999.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, $
	'Syntax: AIHK_MnemTime, mnem, tms1 [, tms2]'

;
;			Get the mnemonic array index and type.
;
AIHK_Mnemonic, mnem, index, arr, Status=status
If (status EQ 0) Then Return
;
;			Convert the index into the time offsets.
;
ticktime = 256L				; tenths of a millisecond.
;
If (arr EQ 1) Then Begin		; PDU array.
	ticks =  9L
	narr  = 97L
EndIf Else Begin			; AEU array (either 1 or 2).
	ticks = 15L
	narr  = 57L
EndElse
;
tms1 = ticks * ticktime * index
tms2 = ticks * ticktime * (index + narr)
;
Return
End
