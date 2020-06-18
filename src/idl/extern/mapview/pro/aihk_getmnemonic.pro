Function AIHK_GetMnemonic, AIHK, Mnem, Status, Counts=cnts, OHMS=Ohms, $
                           Window=Window, Sim_PRT=sim_prt
;+
; NAME:
;	AIHK_GetMnemonic
; PURPOSE:
;	Returns a physical value associated with a mnemonic, extracting
;	the data out of a sweep of analog instrument housekeeping (AIHK)
;	telemetry data.
; CALLING SEQUENCE:
;	value = AIHK_GetMnemonic(AIHK, Mnem [, Status])
; INPUTS:
;	AIHK   - The packet data associated with a single sweep of AIHK telemetry.
;	Mnem   - The name of the AIHK mnemonic.
; OUTPUTS:
;	Status - An optional status value: 0=success, 1=undefined mnemonic.
; RETURNED:
;	value  - The converted value associated with the mnemonic, in
;	         the appropriate physical units.
; OPTIONAL INPUT KEYWORDS:
;	/Counts - If present and nonzero, the counts are returned instead
;	         of physical units.
;	/Ohms   - If present and nonzero, ohms are returned instead Kelvin
;	         for the high resolution PRT signals
;       /Sim_PRT- If present and nonzero, a simulated PRT ohms to K conversion
;                curve is used (fake_ihk).
; OPTIONAL OUTPUT KEYWORD:
;       Window - scalar window value for the supplied PRT
; COMMENTS:
;	AIHK_Mnemonic converts the mnemonic name into an index into the
;	array; this routine applies the appropriate conversion.
; EXAMPLE:
;       Return the values of 'DFW222B5DNI' for the first sweep of file
;      (1975 elements) in the file MAP_tod_20022162357_20022172357.fits
;
;       IDL> fits_read_tod,'MAP_tod_20022162357_20022172357.fits',tod
;       IDL> tt = aihk_getmnemonic(tod.aeu.sweep1,'DFW222B5DNI')
; PROCEDURES USED:
;      AIHK_Mnemonic, AIHK_Mnem_coefs, AIHK_Mnem2Serial,  Get_PRT_Temp
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 16 November 1998.
;	Corrected a logic and an arithmetic error in Windowed mnemonic
;	   calculations.  Changed the Windows mnemonic coefficients to
;	   deal with differences between the two Window arrays.  
;	   MRG, RITSS, 08 January 1999.
;	Counts keyword added.  MRG, RITSS, 11 February 1999.
;       Sim_PRT keyword added.  JW, June 2000.
;       Window keyword added as per GH version.  JW, July 2000.
;       Added bitmask support, JW, August 2000.
;
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, $
	'Syntax: value = AIHK_GetMnemonic(AIHK, Mnem [, Status])'
;
;			Get the index of the mnemonic in the array.  If not
;			found, set the status appropriately and return.
;
AIHK_Mnemonic, Mnem, Index, Arr
If (Index LT 0) Then Begin
	Status = 1
	Return, 0.0d
EndIf Else Begin
	Status = 0
EndElse
;
;			Get the raw telemetry value.
;
Win  = 0
Wflg = 0
Case (Arr) Of
	   3 :	Begin
		Tmp = AIHK.AEU_Analog2[Index]
		If ((Index GE 1) AND (Index LE 31)) Then Begin
		   Wflg = 1
		   Win = AIHK.AEU_Window2[Index / 2]
		   If ((Index MOD 2) NE 0) Then Win = (Win AND '00FF'xL) $
		                           Else Win = ISHFT((Win AND 'FF00'xL), -8)
		   Slp  = 254.968244D-06
		   YInt = 319.5004D0
		   WInt = 129L
		   WSlp = 256.0D0
		   Rmax = 650.25838D0
		EndIf
		End
;
;
	   2 :	Begin
		Tmp = AIHK.AEU_Analog1[Index]
		If ((Index GE 1) AND (Index LE 31)) Then Begin
		   Wflg = 2
		   Win = AIHK.AEU_Window1[Index / 2]
		   If ((Index MOD 2) NE 0) Then Win = (Win AND '00FF'xL) $
		                           Else Win = ISHFT((Win AND 'FF00'xL), -8)
		   Slp  = 255.381467D-06
		   YInt = 319.5197D0
		   WInt = 129L
		   WSlp = 256.0D0
		   Rmax = 650.58226D0
		EndIf
		End
;
;
	Else :	Begin
		Tmp = AIHK.PDU_Analog[Index]
		End
EndCase
;
;			Convert the value into physical units.
;

Window = Win
If (keyword_set(cnts)) Then  Return, Tmp		; Return counts.

If (Wflg EQ 0) Then Begin
;
;				A straightforward polynomial conversion.
;
	AIHK_Mnem_Coefs, Mnem, Coefs, BitMsk
        If (BitMsk ne '0'xL) then Tmp = double(long(Tmp) AND BitMsk)
	n = n_elements(coefs) - 1
	T = double(Tmp)
	Val = 0.0D
	For i = 0, n Do Val = Val + (Coefs[i] * (T ^ i))
;
EndIf Else Begin
;
;				Window'ed telemetry.
;
	Res = (DOUBLE(Tmp) * Slp) + YInt + ((DOUBLE(Win - WInt) / WSlp) * Rmax)


        IF( KEYWORD_SET(Ohms) )THEN RETURN, Res		; Return ohms.
        AIHK_Mnem2Serial, Mnem, Serial_No
        Get_PRT_Temp, Serial_No, Res, Val, sim_prt=sim_prt
;
EndElse
;
Return, Val
End
