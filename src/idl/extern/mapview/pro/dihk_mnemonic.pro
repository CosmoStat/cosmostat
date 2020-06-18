Function DIHK_Mnemonic, Mnem, Status=status
;+
; NAME:
;       DIHK_Mnemonic
; PURPOSE:
;       Returns the array index of a digital instrument housekeeping
;	(DIHK) telemetry mnemonic.
; CALLING SEQUENCE:
;       index = DIHK_Mnemonic(Mnem)
; INPUTS:
;       Mnem  - The name of the DIHK mnemonic.
; RETURNED:
;       Index - The array index.  -1 indicates not found.
; KEYWORD:
;	Status - Returns a status code: 0=not found, 1=found.
; COMMENTS:
;       A giant case statement does the work here.
; MODIFICATION HISTORY:
;       Written automatically from a telemetry list.  08/22/2000 16:38
;-
on_error, 2
;
;                       Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  index = DIHK_Mnemonic(Mnem)'
;
;                       Determine the Fortran-style array index.  This
;                       allows the case statement to perform precisely
;                       like the Fortran version of this routine, simplifying
;                       maintenance.
;
status = 1
Case (strtrim(strupcase(Mnem), 2)) Of
	'ARRAY3B'            : Index =   1
	'DHKSTWD0'           : Index =   1
	'DAEUGAINST'         : Index =   1
	'DA1WINALGST'        : Index =   1
	'DA2WINALGST'        : Index =   1
	'DSCIPKTTRNST'       : Index =   1
	'DHKSTWDSPARE'       : Index =   1
	'DHKSTWD1'           : Index =   2
	'DHKSTWD2'           : Index =   3
	'DHKSTWD3'           : Index =   4
	'DHKSTWD4'           : Index =   5
	'DPWRONRESET'        : Index =   5
	'DWCHDOGRESET'       : Index =   5
	'DSFTWARERESET'      : Index =   5
	'DHKSTWD4DSPARE'     : Index =   5
	'DFLTSWSUBVER'       : Index =   5
	'DFLTSWVER'          : Index =   5
	'DHKSYSERRWD0'       : Index =   6
	'DBADAPIDERR'        : Index =   6
	'DASRAMLDCMDERR'     : Index =   6
	'DAWINCMDERR'        : Index =   6
	'DPHEMTCMDERR'       : Index =   6
	'DASETGNCMDERR'      : Index =   6
	'DAHKTBCMDERR'       : Index =   6
	'DPHEMTVLDCMDERR'    : Index =   6
	'DSCIPKTCMDERR'      : Index =   6
	'DAOVLPLDCMDERR'     : Index =   6
	'DAHFWINLDCMDERR'    : Index =   6
	'DCHKSMCMDERR'       : Index =   6
	'DFUNCODECMDERR'     : Index =   6
	'DCMDPKTBUFFERR'     : Index =   6
	'DPCMDQERR'          : Index =   6
	'DA1CMDQERR'         : Index =   6
	'DA2CMDQERR'         : Index =   6
	'DHKSYSERRWD1'       : Index =   7
	'DHKERRWD1SPARE'     : Index =   7
	'DPEPRMLDCMDERR'     : Index =   7
	'DPBIASLDWOENCE'     : Index =   7
	'DPHEMTONOFWOENCE'   : Index =   7
	'DPHEMTVEDLDWOENCE'  : Index =   7
	'DPEPRMLDWOENCE'     : Index =   7
	'DLSBDRPCMDERR'      : Index =   7
	'DPDUTLMCMDERR'      : Index =   7
	'DBLKCHCMDERR'       : Index =   7
	'DBLKWOENCE'         : Index =   7
	'DBLKDLYCMDERR'      : Index =   7
	'DHKERRWD2SPARE'     : Index =   7
	'DPEPRMRDERR'        : Index =   8
	'DSRAMLOCERR'        : Index =   9
	'DHKSYSERRWD4'       : Index =  10
	'DROSAHKBUFFERR'     : Index =  10
	'DROSAHKRTERR'       : Index =  10
	'DROSSCI1BUFFERR'    : Index =  10
	'DROSSCI1RTERR'      : Index =  10
	'DROSSCI2BUFFERR'    : Index =  10
	'DROSSCI2RTERR'      : Index =  10
	'DROSDHKBUFFERR'     : Index =  10
	'DROSDHKRTERR'       : Index =  10
	'DHKERRWD4SPARE'     : Index =  10
	'DHKSYSERRWD5'       : Index =  11
	'DHKSYSERRWD6'       : Index =  12
	'DHKSYSERRWD7'       : Index =  13
	'DA1WINOVLP'         : Index =  14
	'DA2WINOVLP'         : Index =  15
	'DCMDRCVDCTR'        : Index =  16
	'DCMDEXECTR'         : Index =  17
	'DCMDRJTCTR'         : Index =  18
	'DINT1T'             : Index =  19
	'DINT2T'             : Index =  20
	'DINTP5V'            : Index =  21
	'DINTP15V'           : Index =  22
	'DINTN15V'           : Index =  23
	'DA1HFWINOFFSET'     : Index =  24
	'DA2HFWINOFFSET'     : Index =  25
	'DDIGHKSPARE1'       : Index =  26
	'DDIGHKSPARE2'       : Index =  27
	'DDIGHKSPARE3'       : Index =  28
	'DDIGHKSPARE4'       : Index =  29
	'DDIGHKSPARE5'       : Index =  30
	'DDIGHKSPARE6'       : Index =  31
	'DDIGHKSPARE7'       : Index =  32
	'DDIGHKSPARE8'       : Index =  33
	'DDIGHKSPARE9'       : Index =  34
	'DCALREF'            : Index =  35
	'D2_5V'              : Index =  36
	'DDIGHKSPARE12'      : Index =  37
	'DDIGHKSPARE13'      : Index =  38
	'DDIGHKSPARE14'      : Index =  39
	'DDIGHKSPARE15'      : Index =  40
	'DDIGHKSPARE16'      : Index =  41
	'DDIGHKSPARE17'      : Index =  42
	'DDHKCHKSUM'         : Index =  43
;
	Else : Begin & Index = 0 & status = 0 & End
EndCase
;
;                       Convert the index into an IDL-style array 
;                       index (starting at 0) and return.
;
Return, Index - 1
End
