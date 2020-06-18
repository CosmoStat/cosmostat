Pro DIHK_Mnem_Coefs, Mnem, Coefs, BitMsk, Status=status
;+
; NAME:
;       DIHK_Mnem_Coefs
; PURPOSE:
;       Returns the conversion coefficients associated with an
;	digital instrument housekeeping (DIHK) mnemonic.
; CALLING SEQUENCE:
;       DIHK_Mnem_Coefs, Mnem, Coefs [, BitMsk, STATUS = ]
; INPUTS:
;       Mnem  - The name of the DIHK mnemonic.
; OUTPUTS:
;       Coefs - The coefficients array:
;               convval = sum( Coefs(ind) * X ^ ind )
;	BitMsk- The bitmask for the mnemonic.  Optional
; COMMENTS:
;       A giant case statement does the work here.
; OPTIONAL OUTPUT KEYWORD:
;	Status - Returns a status code: 0=not found, 1=found.
; MODIFICATION HISTORY:
;       Written automatically from a telemetry list.  08/22/2000 16:38
;-
on_error, 2
;
;                       Check arguments.
;
If (n_params() LT 2) Then message, 'Syntax:  DIHK_Mnem_Coefs, Mnem, Coefs [, BitMsk]'
;
;                       Determine the coefficients.
;
status = 1
Case (strtrim(strupcase(Mnem), 2)) Of
;
;			PDU Analog array
;
	'P03BPVNO'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BPCKT'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1000'xL & End
	'P03BSHDF'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '800'xL & End
	'P03BID'             : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BSEGF'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BSCNT'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BPLEN'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BSTIME'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'P03BMTIME'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'ARRAY3B'            : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSTWD0'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAEUGAINST'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1'xL & End
	'DA1WINALGST'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2'xL & End
	'DA2WINALGST'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '4'xL & End
	'DSCIPKTTRNST'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8'xL & End
	'DHKSTWDSPARE'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSTWD1'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSTWD2'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSTWD3'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSTWD4'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPWRONRESET'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1'xL & End
	'DWCHDOGRESET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2'xL & End
	'DSFTWARERESET'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8'xL & End
	'DHKSTWD4DSPARE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFLTSWSUBVER'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFLTSWVER'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSYSERRWD0'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DBADAPIDERR'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1'xL & End
	'DASRAMLDCMDERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2'xL & End
	'DAWINCMDERR'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '4'xL & End
	'DPHEMTCMDERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8'xL & End
	'DASETGNCMDERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '10'xL & End
	'DAHKTBCMDERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '20'xL & End
	'DPHEMTVLDCMDERR'    : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '40'xL & End
	'DSCIPKTCMDERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '80'xL & End
	'DAOVLPLDCMDERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '100'xL & End
	'DAHFWINLDCMDERR'    : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '200'xL & End
	'DCHKSMCMDERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '400'xL & End
	'DFUNCODECMDERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '800'xL & End
	'DCMDPKTBUFFERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1000'xL & End
	'DPCMDQERR'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2000'xL & End
	'DA1CMDQERR'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '4000'xL & End
	'DA2CMDQERR'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8000'xL & End
	'DHKSYSERRWD1'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKERRWD1SPARE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1'xL & End
	'DPEPRMLDCMDERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2'xL & End
	'DPBIASLDWOENCE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '4'xL & End
	'DPHEMTONOFWOENCE'   : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8'xL & End
	'DPHEMTVEDLDWOENCE'  : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '10'xL & End
	'DPEPRMLDWOENCE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '20'xL & End
	'DLSBDRPCMDERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '40'xL & End
	'DPDUTLMCMDERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '80'xL & End
	'DBLKCHCMDERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '100'xL & End
	'DBLKWOENCE'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '200'xL & End
	'DBLKDLYCMDERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '400'xL & End
	'DHKERRWD2SPARE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPEPRMRDERR'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DSRAMLOCERR'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSYSERRWD4'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DROSAHKBUFFERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '1'xL & End
	'DROSAHKRTERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '2'xL & End
	'DROSSCI1BUFFERR'    : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '4'xL & End
	'DROSSCI1RTERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '8'xL & End
	'DROSSCI2BUFFERR'    : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '10'xL & End
	'DROSSCI2RTERR'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '20'xL & End
	'DROSDHKBUFFERR'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '40'xL & End
	'DROSDHKRTERR'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '80'xL & End
	'DHKERRWD4SPARE'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSYSERRWD5'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSYSERRWD6'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DHKSYSERRWD7'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DA1WINOVLP'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DA2WINOVLP'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DCMDRCVDCTR'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DCMDEXECTR'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DCMDRJTCTR'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DINT1T'             : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DINT2T'             : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DINTP5V'            : Begin & Coefs = [ -1.99885000d+01,   9.76600000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DINTP15V'           : Begin & Coefs = [ -5.25496000d+01,   2.56590000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DINTN15V'           : Begin & Coefs = [ -5.25496000d+01,   2.56590000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DA1HFWINOFFSET'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DA2HFWINOFFSET'     : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE1'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE2'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE3'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE4'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE5'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE6'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE7'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE8'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE9'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DCALREF'            : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'D2_5V'              : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE12'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE13'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE14'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE15'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE16'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDIGHKSPARE17'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DDHKCHKSUM'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
;
	Else : Begin & Coefs = [0d, 1d0, 0d0] & BitMsk = 0L & status = 0 & End
EndCase
;
Return
End
