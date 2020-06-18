Function AIHK_Pckt2Mnemonic, AIHK, Mnem, Swp, Status, Counts=cnts,      $
                             Juls=JulTime, Ohms=ohms, Window=Window,    $
                             Sim_PRT=sim_prt
;+
; NAME:
;	AIHK_Pckt2Mnemonic
; PURPOSE:
;	Returns the physical value(s) associated with an analog
;	instrument housekeeping (AIHK) mnemonic, extracting 
;	the data from a AIHK telemetry packet(s).
; CALLING SEQUENCE:
;	value = AIHK_Pckt2Mnemonic (AIHK, Mnem, Swp [, Status])
; INPUTS:
;	AIHK   - The packet(s) of AIHK telemetry.
;	Mnem   - The name of the mnemonic.
;	Swp    - Which of the two AIHK sweeps to extract from.  If <=0,
;	         both sweeps are extracted (and folded together in time
;	         order [1,2,1,2,1,2,...].  If 1, sweep 1 is processed; 
;	         else sweep 2 is processed.
; OUTPUTS:
;	Status - An optional argument returning a status code:
;	         0=success, 1=undefined mnemonic.
; RETURNED:
;	value  - The value(s) associated with the mnemonic, with the conversion
;	         applied.
; OPTIONAL INPUT KEYWORDS:
;	/Counts - If present and nonzero, the counts are returned instead
;	         of physical units.
;       /Ohms   - If present and nonzero, prt resistance in ohms is returned
;                rather than temperature in Kelvin.  Note that the use of
;                Ohms and Counts is mutually exclusive.
;       Sim_PRT- If present and nonzero, conversion from ohms to K is performed
;                using the fake_ihk conversion curve, which is independent of
;                serial no.
; OPTIONAL OUTPUT KEYOWRD:
;       Juls -  vector of Julian times of the mnemonic data points
;       Window - scalar window value for the PRT
; COMMENTS:
;	This routine simply sets up the call to AIHK_GetMnemonic, while 
;	providing support for an array of telemetry packets.
; EXAMPLE:
;       Return the values of 'DFW222B5DNI' for the first sweep of file
;      (1975 elements) in the file MAP_tod_20022162357_20022172357.fits
;
;       IDL> fits_read_tod,'MAP_tod_20022162357_20022172357.fits',tod
;       IDL> tt = aihk_pckt2mnemonic(tod.aeu,'DFW222B5DNI',1)
; PROCEDURES CALLED:
;       AIHK_GetMnemonic(), TEL2JUL()
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 16 November 1998
;	Counts keyword added.  Array support.  MRG, RITSS, 11 February 1999.
;       Ohms keyword added for consistency with aihk_getmnemonic. 
;                                                            JW, 25 Oct 1999
;       Sim_PRT keyword added. JW, June 2000
;       Added Window keyword and commented out unneeded loop over packets 
;           as per GH code.  JW, July 2000.
;       Juls keyword added.  MRG, SSAI, 12 December 2001.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, $
	'Syntax:  value = AIHK_Pckt2Mnemonic (AIHK, Mnem, Swp [, Status])'
;
c = keyword_set(cnts)
oh = keyword_set(ohms)
fake_prt = keyword_set(sim_prt)
;
;			Get the data.
;
If (Swp LE 1) Then Begin
    Val1 = AIHK_GetMnemonic(AIHK.Sweep1, Mnem, Status, $
                            Counts=c, Ohms=oh, Window=w1, sim_prt=fake_prt)
    Jul1 = tel2jul(AIHK)
EndIf
If (Swp NE 1) Then Begin
    Val2 = AIHK_GetMnemonic(AIHK.Sweep2, Mnem, Status, $
                            Counts=c, Ohms=oh, Window=w2, sim_prt=fake_prt)
    Jul2 = tel2jul(AIHK) + (23.04D0 / (24.d * 3600.d))
EndIf
;
;			Return the requested data.
;
If (Swp EQ 1) Then Begin
	Val = Val1
        Window = w1
        JulTime = Jul1
EndIf Else If (Swp GT 1) Then Begin
	Val = Val2
        Window = w2
        JulTime = Jul2
EndIf Else Begin
	Val = transpose([[Val1], [Val2]])
	Val = Val[*]
        Window = transpose([[w1],[w2]])
        Window = Window[*]
        JulTime = transpose([[Jul1], [Jul2]])
        JulTime = JulTime[*]
EndElse
;
If (n_elements(Val) EQ 1L) Then begin
    Val = Val[0]
    Window = Window[0]
    JulTime = JulTime[0]
Endif
;
Return, Val
End
