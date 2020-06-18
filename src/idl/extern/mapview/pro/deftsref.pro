Function DefTSRef, dum
;+
; NAME:
;	DefTSRef
; PURPOSE:
;	Returns the default WMAP timestamp reference GMT.
; CALLING SEQUENCE:
;	refgmt = DefTSRef()
; RETURNED:
;	refgmt - The default WMAP timestamp reference GMT.
; COMMENTS:
;	This function allows this reference to be maintained in
;	one location, simplifying maintenance.
;
;	The reference date is set to the start of the LAUNCH DAY
;	instead of the precise MET time definition; this prevents
;	negative time stamps.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon ITSS, 12 January 1999.
;	Flight update.  MRG, RITSS, 03 July 2001.
;-
on_error, 2
;
Return, '20011810000000000000'
End
