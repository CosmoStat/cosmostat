FUNCTION Sci_Mnemonic, Mnemonic, Status
;+
; NAME:
;	Sci_Mnemonic
; PURPOSE:
;	Tests whether or not Mnemonic is a valid science mnemonic for use
;	with Sci_GetMnemonic.
; CALLING SEQUENCE:
;	Test = Sci_Mnemonic(Mnemonic [,Status])
; INPUTS:
;	Mnem	The name of the science mnemonic.  The mnemonic naming
;		conventions are as follows: let 'chan' be a radiometer 
;		channel string, such as 'K113', etc.  then:
;		'DK113'  returns all observations in time-order.
;		'QDK113' returns the first observation in each packet
;		'ADK113' returns the major frame (packet) average
;		'RDK113' returns the major frame (packet) rms
;	Status	Optional status return argument, equivalent to the returned
;		function value.
; OUTPUTS:
;	Status	Test result: 0=defined mnemonic, -1=undefined mnemonic.
; MODIFICATION HISTORY:
;	Gary Hinshaw, NASA/GSFC, 4 October 1999
;-
ON_ERROR, 2
;
IF( N_PARAMS() LT 1 )THEN MESSAGE, $
  'Syntax:  Test = Sci_Mnemonic(Mnemonic [,Status])

Status = -1
Mnem = STRUPCASE(STRTRIM(Mnemonic,2))
Mnem_Type = STRMID(Mnem,0,1)

CASE Mnem_Type OF
 'D' : Channel = STRMID(Mnemonic,1,100)
 'Q' : Channel = STRMID(Mnemonic,2,100)
 'A' : Channel = STRMID(Mnemonic,2,100)
 'R' : Channel = STRMID(Mnemonic,2,100)
 ELSE : RETURN, Status
ENDCASE

; Extract the band and index from the supplied channel
Extract_Band_Index, Channel, Band, Index, Status

RETURN, Status
END
