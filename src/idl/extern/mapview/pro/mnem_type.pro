Function Mnem_Type, Mnemonic
;+
; NAME:
;	Mnem_Type
; PURPOSE:
; 	Function to return the mnemonic type associated with Mnemonic.  
;	The result is either 'AIHK', 'DIHK', 'SCI' or 'UNDEFINED'.
; CALLING SEQUENCE:
;	Type = Mnem_Type(Mnemonic)
; INPUTS:
;	Mnemonic - The mnemonic(s) to test.
; ReturnED:
;	Type     - The type associated with each mnemonic.
; PROCEDURES CALLED:
;       AIHK_MNEMONIC, DIHK_MNEMONIC()
; MODIfICATION HISTORY:
;	Written by Al Kogut, NASA/GSFC
;	Updated documentation.  Multiple mnemonic support.  MRG, RITSS,
;	   09 January 2001.
;       Remove spacecraft and GPIB mnemonics   W. Landsman January 2003
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  Type = Mnem_Type(Mnemonic)'
;
; 			Create the output array, initialized with "UNDEFINED"
;
n = n_elements(Mnemonic)
If (n LE 0) Then message, 'At least one mnemonic must be supplied!'
Type = replicate('UNDEFINED', n)
;
;			Interrogate each mnemonic.
;
n = n - 1
For i = 0, n Do Begin
;
;				Analog instrument housekeeping.
;
  AIHK_Mnemonic, Mnemonic[i], Index, Array, Status=Stat
  If (Stat GT 0) Then Begin
    Type[i] = 'AIHK'
    Continue
  EndIf
;
;				Digital instrument housekeeping.
;
  If (DIHK_Mnemonic(Mnemonic[i]) GE 0) Then Begin
    Type[i] = 'DIHK' 
    Continue
  EndIf
;
;				Science.
;
  If (Sci_Mnemonic(Mnemonic[i]) GE 0) Then Begin
    Type[i] = 'SCI' 
    Continue
  EndIf
;
;
EndFor
;
Return, Type
End
