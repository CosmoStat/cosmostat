Pro AIHK_Mnem_Coefs, Mnem, Coefs, BitMsk, Status=status
;+
; NAME:
;       AIHK_Mnem_Coefs
; PURPOSE:
;       Returns the conversion coefficients associated with an
;	analog instrument housekeeping (AIHK) mnemonic.
; CALLING SEQUENCE:
;       AIHK_Mnem_Coefs, Mnem, Coefs [, BitMsk]
; INPUTS:
;       Mnem  - The name of the AIHK mnemonic.
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
If (n_params() LT 2) Then message, 'Syntax:  AIHK_Mnem_Coefs, Mnem, Coefs [, BitMsk]'
;
;                       Determine the coefficients.
;
status = 1
Case (strtrim(strupcase(Mnem), 2)) Of
;
;			PDU Analog array
;
	'DFQ111B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV111B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA111B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW111B4DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW211B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW311B6DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW411B7DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK111B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV211B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ211B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPPINTT1'           : Begin & Coefs = [ -2.31500000d+01,   4.06901000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPPINTT2'           : Begin & Coefs = [ -2.31500000d+01,   4.06901000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPPINTT3'           : Begin & Coefs = [ -2.31500000d+01,   4.06901000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPRM7_2V'           : Begin & Coefs = [ -9.54545000d+00,   4.66087000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFLDAM6_2V'        : Begin & Coefs = [ -8.33333000d+00,   4.06901000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ112B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV112B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA112B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW112B4DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW212B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW312B6DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW412B7DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK112B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV212B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ212B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPRP7_2V'           : Begin & Coefs = [ -9.54545000d+00,   4.66087000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFLDBM6_2V'        : Begin & Coefs = [ -8.33333000d+00,   4.06901000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ121B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV121B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA121B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW121B4DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW221B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW321B6DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW421B7DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK121B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV221B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ221B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFLEDP10V'         : Begin & Coefs = [ -1.32500000d+01,   6.46973000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPPHSWCONVP9V'      : Begin & Coefs = [ -1.20000000d+01,   5.85938000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ122B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV122B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA122B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW122B4DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW222B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW322B6DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW422B7DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK122B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV222B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ222B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFM7_2V'           : Begin & Coefs = [ -9.54545000d+00,   4.66087000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPPHSWCONVM9V'      : Begin & Coefs = [ -1.20000000d+01,   5.85938000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ111B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV111B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA111B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW111B4DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW211B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW311B6DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW411B7DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK111B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV211B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ211B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFP7_2V'           : Begin & Coefs = [ -9.54545000d+00,   4.66087000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPHKM15V'           : Begin & Coefs = [ -1.71242000d+01,   8.36145000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ112B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV112B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA112B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW112B4DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW212B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW312B6DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW412B7DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK112B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV212B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ212B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFLDAP6_2V'        : Begin & Coefs = [ -8.33333000d+00,   4.06901000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPHKP15V'           : Begin & Coefs = [ -1.71242000d+01,   8.36145000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ121B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV121B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA121B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW121B4DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW221B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW321B6DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW421B7DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK121B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV221B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ221B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPHKP5V'            : Begin & Coefs = [ -5.56112000d+00,   2.71539000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ122B1DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV122B2DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA122B3DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW122B4DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW222B5DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW322B6DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW422B7DNI'        : Begin & Coefs = [ -6.41026000d+01,   3.13001000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK122B8DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV222B9DNI'        : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ222B10DNI'       : Begin & Coefs = [ -4.16667000d+01,   2.03451000d-02,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPFLDBP6_2V'        : Begin & Coefs = [ -8.33333000d+00,   4.06901000d-03,   0.00000000d+00] & BitMsk = '0'xL & End
;
;			First AEU Analog sweep array
;
	'DABD1V'             : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DARREF1BD1'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAW323_4AMPT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAV113_4ADT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAW2_14_23AMP_ADT'  : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAIHK1BDT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DACONVBDT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPV221_2RXBT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPW111_2RXBT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPW321_2FPAT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV222RXBAMPT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW111RXBAMPT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW221RXBAMPT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA12RXBRIBT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ2RXBRIBT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRPYPSHPRTKT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DARREF2BD1'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK1AFEEDT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ1AFEEDT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW3BFEEDT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFK1BOMTT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ1BOMTT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW3AOMTT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV22FPATEET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW32FPATEET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTATOPSECT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTABOTSECT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBMIDSECT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBTOPPRIT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTAMIDPRIT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTAPXMIDRADT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBMXBOTRADT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK113RFBI0'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK114RFBI1'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK123RFBI2'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK124RFBI3'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW413RFBI8'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW414RFBI9'        : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW423RFBI10'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW424RFBI11'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW313RFBI16'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW314RFBI17'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW323RFBI18'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW324RFBI19'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW213RFBI24'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW214RFBI25'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW223RFBI26'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW224RFBI27'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV113RFBI32'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV114RFBI33'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV123RFBI34'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV124RFBI35'       : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP15VBD1'          : Begin & Coefs = [ -2.99513000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAM15VBD1'          : Begin & Coefs = [ -2.99513000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP12VBD1'          : Begin & Coefs = [ -2.99513000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAM12VBD1'          : Begin & Coefs = [ -2.99513000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP5VBD1'           : Begin & Coefs = [ -2.99513000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
;
;			Second AEU Analog sweep array
;
	'DABD2V'             : Begin & Coefs = [ -4.99188000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DARREF1BD2'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAW113_4ADT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAV223_4AMPT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAQ113_4ADT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAIHK2BDT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DASPARE1'           : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPV111_2FPAT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPW321_2RXBT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DPW221_2FPAT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV111RXBAMPT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW321RXBAMPT'      : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRK12RXBRIBT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ1RXBRIBT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW3RXBRIBT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRMYPSHPRTKT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DARREF2BD2'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA1BFEEDT'        : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ2BFEEDT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW3AFEEDT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFKA1AOMTT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFQ2AOMTT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW3BOMTT'          : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFV11FPATEET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW11FPATEET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DFW22FPATEET'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTAMIDSECT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBTOPSECT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTATOPPRIT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBMIDPRIT'         : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTBPXMIDRADT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DTAMXTOPRADT'       : Begin & Coefs = [  0.00000000d+00,   1.00000000d+00,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW113RFBI4'        : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW114RFBI5'        : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW123RFBI2'        : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRW124RFBI3'        : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV213RFBI12'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV214RFBI13'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV223RFBI14'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRV224RFBI15'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ113RFBI20'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ114RFBI21'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ123RFBI22'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ124RFBI23'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ213RFBI28'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ214RFBI29'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ223RFBI30'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRQ224RFBI31'       : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA113RFBI36'      : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA114RFBI37'      : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA123RFBI38'      : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DRKA124RFBI39'      : Begin & Coefs = [ -4.98047000d+00,   1.56250000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP15VBD2'          : Begin & Coefs = [ -2.98828000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAM15VBD2'          : Begin & Coefs = [ -2.98828000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP12VBD2'          : Begin & Coefs = [ -2.98828000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAM12VBD2'          : Begin & Coefs = [ -2.98828000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
	'DAP5VBD2'           : Begin & Coefs = [ -2.98828000d+01,   9.37500000d-04,   0.00000000d+00] & BitMsk = '0'xL & End
;
	Else : Begin & Coefs = [0d, 1d0, 0d0] & BitMsk = 0L & status = 0 & End
EndCase
;
Return
End
