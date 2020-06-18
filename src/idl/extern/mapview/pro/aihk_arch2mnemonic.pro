Function AIHK_Arch2Mnemonic, TOD, Mnem, Swp, Status, Counts=cnts
;+
; NAME:
;       AIHK_Arch2Mnemonic
; PURPOSE:
;       Returns the physical value(s) associated with an analog
;       instrument housekeeping (AIHK) mnemonic, extracting 
;       the data from a time-ordered archive record(s).
; CALLING SEQUENCE:
;       value = AIHK_Arch2Mnemonic (TOD, Mnem, Swp [, Status])
; INPUTS:
;       TOD    - The time-ordered archive record(s).
;       Mnem   - The name of the mnemonic.
;       Swp    - Which of the two AIHK sweeps to extract from.  If <=0,
;                both sweeps are extracted (and folded together in time
;                order [1,2,1,2,1,2,...].  If 1, sweep 1 is processed; 
;                else sweep 2 is processed.
; OUTPUTS:
;       Status - An optional argument returning a status code:
;                0=success, 1=undefined mnemonic.
; RETURNED:
;       value  - The value(s) associated with the mnemonic, with the conversion
;                applied.
; OPTIONAL INPUT KEYWORDS:
;       /Counts - If present and nonzero, the counts are returned instead
;                of physical units.
; COMMENTS:
;       This routine simply sets up the call to AIHK_GetMnemonic, while 
;       providing support for an array of telemetry packets.
; EXAMPLE:
;       Print the 3750 values of 'DFQ211B10DNI' for both sweeps of day 220 of 
;      2002
;
;       IDL> file ='MAP_tod_20022202357_20022212357.fits'
;       IDL> fits_read_tod,file,arch     ;Read TOD structure
;       IDL> print,AIHK_arch2Mnemonic (arch, 'DFQ211B10DNI',-1)
; MODIFICATION HISTORY:
;       Michael R. Greason, Raytheon STX, 16 November 1998
;       Counts keyword added.  Array support.  MRG, RITSS, 11 February 1999.
;       Vectorized. MRG & JP, RITSS, 18 October 2000.
;-
on_error, 2
;
;                       Check arguments.
;
If (n_params() LT 3) Then message, $
        'Syntax:  value = AIHK_Arch2Mnemonic (TOD, Mnem, Swp [, Status])'
c = keyword_set(cnts)
;
;                       Initialize the two sweep output arrays.
;
n = n_elements(TOD)
If (c) Then Begin
        Val1 = lonarr(n)
        Val2 = Val1
EndIf Else Begin
        Val1 = dblarr(n)
        Val2 = Val1
EndElse
;
;                       Get the data.
;
n = n - 1
If (Swp LE 1) Then Begin
        Val1 = AIHK_GetMnemonic(TOD.aeu.Sweep1, Mnem, Status, Counts=c)
Endif
If (Swp NE 1) Then Begin
        Val2 = AIHK_GetMnemonic(TOD.aeu.Sweep2, Mnem, Status, Counts=c)
Endif
;
;                       Return the requested data.
;
If (Swp EQ 1) Then Begin
        Val = Val1
EndIf Else If (Swp GT 1) Then Begin
        Val = Val2
EndIf Else Begin
        Val = transpose([[Val1], [Val2]])
        Val = Val[*]
EndElse
;
Return, Val
End
