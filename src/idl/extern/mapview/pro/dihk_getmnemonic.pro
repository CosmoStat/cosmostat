Function DIHK_GetMnemonic, DIHK, Mnem, Status, Counts=cnts
;+
; NAME:
;	DIHK_GetMnemonic
; PURPOSE:
;	Returns a physical value associated with a mnemonic, extracting
;	the data out of an array of digital instrument housekeeping
;	(DIHK) telemetry data.
; CALLING SEQUENCE:
;	value = DIHK_GetMnemonic(DIHK, Mnem [, Status, /COUNTS])
; INPUTS:
;	DIHK   - The packet data array of DIHK telemetry.
;	Mnem   - The name of the DIHK mnemonic.
; OUTPUTS:
;	Status - An optional status value: 0=success, 1=undefined mnemonic.
; RETURNED:
;	value  - The converted value associated with the mnemonic, in
;	         the appropriate physical units.
; OPTIONAL INPUT KEYWORD:
;	/Counts - If present and nonzero, the value is in integer counts
;	         instead of double precision physical units.
; COMMENTS:
;	DIHK_Mnemonic converts the mnemonic name into an index into the
;	array; this routine applies the appropriate conversion.
; EXAMPLE:
;       Get the DINT1T temperatures from the time-ordered file 
;        MAP_tod_20022182357_20022192357.fits'
;
;       IDL> fits_read_tod,'MAP_tod_20022182357_20022192357.fits',tod
;       IDL> temp = dihk_getmnemonic(tod.deu.data, 'DINT1T')
; PROCEDURES CALLED:
;       DIHK_Mnemonic(), DIHK_Mnem_Coefs GET_DEU_INT_TEMP
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon STX, 12 November 1998
;	Counts added.  MRG, RITSS, 11 February 1999.
;       poly conversion added.  JW, 3 June 1999.
;       special case for DEU internal temps.  JW, 01 Oct 1999
;       added bitmask support.  JW, August 2000.
;       Add Poly() call, work with TOD structure  W. Landsman January 2003
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 2) Then message, $
	'Syntax: value = DIHK_GetMnemonic(DIHK, Mnem [, Status, /COUNTS])'
;
;			Get the index of the mnemonic in the array.  If not
;			found, set the status appropriately and return.
;
Index = DIHK_Mnemonic(Mnem)
If (Index LT 0) Then Begin
	Status = 1
	Return, 0.0d
EndIf Else Begin
	Status = 0
EndElse
;
;			Convert the value into physical units.
;
Tmp = reform(DIHK[Index,*,*])
Case (strtrim(strupcase(Mnem), 2)) Of

;
'DAEUGAINST' : Begin
	Tmp =  (Tmp And '0001'xL)
	End
'DA1WINALGST' : Begin
	Tmp = ISHFT((Tmp And '0002'xL), -1)
	End
'DA2WINALGST' : Begin
	Tmp = ISHFT((Tmp And '0004'xL), -2)
	End
'DSCIPKTTRNST' : Begin
	Tmp = ISHFT((Tmp And '0008'xL), -3)
	End
'DHKSTWDSPARE' : Begin
	Tmp = ISHFT((Tmp And 'FFF0'xL), -4)
	End
;
'DPWRONRESET' : Begin
	Tmp = (Tmp And '0001'xL)
	End
'DWCHDOGRESET' : Begin
	Tmp = ISHFT((Tmp And '0002'xL),  -1)
	End
'DHDWARERESET' : Begin
	Tmp = ISHFT((Tmp And '0004'xL),  -2)
	End
'DSFTWARERESET' : Begin
	Tmp = ISHFT((Tmp And '0008'xL),  -3)
	End
'DHKSTWD4DSPARE' : Begin
	Tmp = ISHFT((Tmp And '00F0'xL),  -4)
	End
'DFLTSWSUBVER' : Begin
	Tmp = ISHFT((Tmp And '0F00'xL),  -8)
	End
'DFLTSWVER' : Begin
	Tmp = ISHFT((Tmp And 'F000'xL), -12)
	End
;
'DBADAPIDERR' : Begin
	Tmp = (Tmp And '0001'xL)
	End
'DASRAMLDCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0002'xL),  -1)
	End
'DAWINCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0004'xL),  -2)
	End
'DPHEMTCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0008'xL),  -3)
	End
'DASETGNCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0010'xL),  -4)
	End
'DAHKTBCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0020'xL),  -5)
	End
'DPHEMTVLDCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0040'xL),  -6)
	End
'DSCIPKTCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0080'xL),  -7)
	End
'DAOVLPLDCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0100'xL),  -8)
	End
'DAHFWINLDCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0200'xL),  -9)
	End
'DCHKSMCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0400'xL), -10)
	End
'DFUNCODECMDERR' : Begin
	Tmp = ISHFT((Tmp And '0800'xL), -11)
	End
'DCMDPKTBUFFERR' : Begin
	Tmp = ISHFT((Tmp And '1000'xL), -12)
	End
'DPCMDQERR' : Begin
	Tmp = ISHFT((Tmp And '2000'xL), -13)
	End
'DA1CMDQERR' : Begin
	Tmp = ISHFT((Tmp And '4000'xL), -14)
	End
'DA2CMDQERR' : Begin
	Tmp = ISHFT((Tmp And '8000'xL), -15)
	End
;
'DPEPRMLDCMDERR' : Begin
	Tmp = (Tmp And '0001'xL)
	End
'DPBIASLDWOENCE' : Begin
	Tmp = ISHFT((Tmp And '0002'xL),  -1)
	End
'DPHEMTONOFWOENCE' : Begin
	Tmp = ISHFT((Tmp And '0004'xL),  -2)
	End
'DPHEMTVEDLDWOENCE' : Begin
	Tmp = ISHFT((Tmp And '0008'xL),  -3)
	End
'DPEPRMLDWOENCE' : Begin
	Tmp = ISHFT((Tmp And '0010'xL),  -4)
	End
'DLSPDRPCMDERR' : Begin
	Tmp = ISHFT((Tmp And '0020'xL),  -5)
	End
'DHKERRWD1SPARE' : Begin
	Tmp = ISHFT((Tmp And 'FFC0'xL),  -6)
	End
;
'DROSAHKBUFFERR' : Begin
	Tmp = (Tmp And '0001'xL)
	End
'DROSAHKRTERR' : Begin
	Tmp = ISHFT((Tmp And '0002'xL),  -1)
	End
'DROSSCI1BUFFERR' : Begin
	Tmp = ISHFT((Tmp And '0004'xL),  -2)
	End
'DROSSCI1RTERR' : Begin
	Tmp = ISHFT((Tmp And '0008'xL),  -3)
	End
'DROSSCI2BUFFERR' : Begin
	Tmp = ISHFT((Tmp And '0010'xL),  -4)
	End
'DROSSCI2RTERR' : Begin
	Tmp = ISHFT((Tmp And '0020'xL),  -5)
	End
'DROSDHKBUFFERR' : Begin
	Tmp = ISHFT((Tmp And '0040'xL),  -6)
	End
'DROSDHKRTERR' : Begin
	Tmp = ISHFT((Tmp And '0080'xL),  -7)
	End
'DHKERRWD4SPARE' : Begin
	Tmp = ISHFT((Tmp And 'FF00'xL),  -8)
	End
;
Else :
EndCase
;
If (keyword_set(cnts)) Then Return, Tmp		; Return counts.

;
; Special case for DEU internal temperatures conversion from cts to C
;
Case (strtrim(strupcase(Mnem), 2)) Of

'DINT1T' : get_deu_int_temp,tmp,val
'DINT2T' : get_deu_int_temp,tmp,val
 Else :  Begin
;
; otherwise, A straightforward polynomial conversion.
;
	DIHK_Mnem_Coefs, Mnem, Coefs, BitMsk
        If (BitMsk ne '0'xL) then Tmp = double(long(Tmp) AND BitMsk)
        val = poly(double(tmp), coefs)
        End

EndCase

Return, Val
End



