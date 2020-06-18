Function DIHK_Pckt2Mnemonic, DIHK, Mnem, Status, Counts=cnts
;+
; NAME:
;	DIHK_Pckt2Mnemonic
; PURPOSE:
;	Returns the physical value associated with a instrument
;	digital telemetry mnemonic, extracting the data from a 
;	telemetry packet.
; CALLING SEQUENCE:
;	value = DIHK_Pckt2Mnemonic (DIHK, Mnem [, Status])
; INPUTS:
;	DIHK   - The packet(s) of DIHK telemetry.
;	Mnem   - The name of the mnemonic.
; OUTPUTS:
;	Status - An optional argument returning a status code:
;	         0=success, 1=undefined mnemonic.
; RETURNED:
;	value  - The value(s) associated with the mnemonic, with the conversion
;	         applied.
; OPTIONAL INPUT KEYWORD:
;	Counts - If present and nonzero, the value is in integer counts
;	         instead of double precision physical units.
; COMMENTS:
;	This routine simply sets up the call to DIHK_GetMnemonic, while
;	providing support for an array of packet data.
; EXAMPLE:
;       Get the DINT1T temperatures from the file 
;        MAP_tod_20022182357_20022192357.fits'
;
;       IDL> fits_read_tod,'MAP_tod_20022182357_20022192357.fits',tod
;       IDL> print,dihk_pckt2mnemonic(tod.deu, 'DINT1T')
;            ==> 23.9877      23.5471      23.9877      23.9877
; PROCEDURES CALLED:
;       DIHK_GetMnemonic()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 12 November 1998.
;	Support for an array of packet data added.  Counts support added.
;	   MRG, RITSS, 11 February 1999
;       Remove packet loop, similar to aihk_pckt2mnemonic. JW, 31 July 2000.
;	
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, $
	'Syntax:  value = DIHK_Pckt2Mnemonic (DIHK, Mnem [, Status])'
c = keyword_set(cnts)

Val = DIHK_GetMnemonic(DIHK.Data, Mnem, Status, Counts=c)
If (n_elements(val) EQ 1) Then Val = Val[0]
;
Return, Val
End
