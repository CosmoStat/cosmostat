Pro AIHK_Mnem2Serial, Mnem, Serial_No, Status=status
;+
; NAME:
;       AIHK_Mnem2Serial
; PURPOSE:
;       Returns the serial number of an analog instrument housekeeping (AIHK) 
;       PRT.
; CALLING SEQUENCE:
;       AIHK_Mnem2Serial, Mnem, Serial_No
; INPUTS:
;       Mnem  - The name of the AIHK mnemonic.
; OUTPUTS:
;       Serial_No - The serial number of the PRT with this mnemonic
; OPTIONAL OUTPUT KEYWORD:
;       Status - Returns a status code: 0=not found, 1=found.
; COMMENTS:
;       A giant case statement does the work here.
; PROCEDURES USED:
;       AIHK_Mnem2Serial_List
; MODIFICATION HISTORY:
;       Written by Gary Hinshaw NASA/GSFC, 03/30/1999 14:05
;       Assigned serial number to DRV111RXBAMPT.
;         MRG, RITSS, 27 Aug 1999.
;       Modified to use AIHK_Mnem2Serial_List to get the Mnemonic/Serial
;         number information.  Michael R. Greason, Raytheon ITSS,
;         02 September 1999
;-
ON_ERROR, 2
;
;                       Check arguments
;
If (N_PARAMS() LT 2) Then MESSAGE, 'Syntax: AIHK_Mnem2Serial, Mnem, Serial_No'
;
;                       Determine the serial number
;
AIHK_Mnem2Serial_List, info
w = where(strupcase(strtrim(info.Mnemonic, 2)) EQ strupcase(strtrim(Mnem, 2)))
If (w[0] LT 0) Then Begin
        Status = 0
EndIf Else Begin
        Status = 1
        Serial_No = info[w[0]].SerialNum
EndElse
;
Return
End
