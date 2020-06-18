Pro AIHK_Mnemonic, Mnem, Index, Arr, Status=status
;+
; NAME:
;       AIHK_Mnemonic
; PURPOSE:
;       Returns the sweep/array index of an analog instrument
;	housekeeping (AIHK) mnemonic.
; CALLING SEQUENCE:
;       AIHK_Mnemonic, Mnem, Index, Arr
; INPUTS:
;       Mnem  - The name of the AIHK mnemonic.
; OUTPUTS:
;       Index - The array index.  -1 indicates not found.
;       Arr   - A flag indicating the sweep element:
;                       1 = PDU array,
;                       2 = AEU 1 array,
;                       3 = AEU 2 array.
; OTPIONAL OUTPUT KEYWORD:
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
If (n_params() LT 3) Then message, 'Syntax:  AIHK_Mnemonic, Mnem, Index, Arr'
;
;                       Determine the Fortran-style array index.  This
;                       allows the case statement to perform precisely
;                       like the Fortran version of this routine, simplifying
;                       maintenance.
;
status = 1
Case (strtrim(strupcase(Mnem), 2)) Of
;
;			PDU Analog array
;
	'DFQ111B1DNI'        : Begin & Index =   1 & Arr =  1 & End
	'DFV111B2DNI'        : Begin & Index =   2 & Arr =  1 & End
	'DFKA111B3DNI'       : Begin & Index =   3 & Arr =  1 & End
	'DFW111B4DNI'        : Begin & Index =   4 & Arr =  1 & End
	'DFW211B5DNI'        : Begin & Index =   5 & Arr =  1 & End
	'DFW311B6DNI'        : Begin & Index =   6 & Arr =  1 & End
	'DFW411B7DNI'        : Begin & Index =   7 & Arr =  1 & End
	'DFK111B8DNI'        : Begin & Index =   8 & Arr =  1 & End
	'DFV211B9DNI'        : Begin & Index =   9 & Arr =  1 & End
	'DFQ211B10DNI'       : Begin & Index =  10 & Arr =  1 & End
	'DPPINTT1'           : Begin & Index =  11 & Arr =  1 & End
	'DPPINTT2'           : Begin & Index =  12 & Arr =  1 & End
	'DPPINTT3'           : Begin & Index =  13 & Arr =  1 & End
	'DPRM7_2V'           : Begin & Index =  14 & Arr =  1 & End
	'DPFLDAM6_2V'        : Begin & Index =  15 & Arr =  1 & End
	'DFQ112B1DNI'        : Begin & Index =  16 & Arr =  1 & End
	'DFV112B2DNI'        : Begin & Index =  17 & Arr =  1 & End
	'DFKA112B3DNI'       : Begin & Index =  18 & Arr =  1 & End
	'DFW112B4DNI'        : Begin & Index =  19 & Arr =  1 & End
	'DFW212B5DNI'        : Begin & Index =  20 & Arr =  1 & End
	'DFW312B6DNI'        : Begin & Index =  21 & Arr =  1 & End
	'DFW412B7DNI'        : Begin & Index =  22 & Arr =  1 & End
	'DFK112B8DNI'        : Begin & Index =  23 & Arr =  1 & End
	'DFV212B9DNI'        : Begin & Index =  24 & Arr =  1 & End
	'DFQ212B10DNI'       : Begin & Index =  25 & Arr =  1 & End
	'DPRP7_2V'           : Begin & Index =  26 & Arr =  1 & End
	'DPFLDBM6_2V'        : Begin & Index =  27 & Arr =  1 & End
	'DFQ121B1DNI'        : Begin & Index =  28 & Arr =  1 & End
	'DFV121B2DNI'        : Begin & Index =  29 & Arr =  1 & End
	'DFKA121B3DNI'       : Begin & Index =  30 & Arr =  1 & End
	'DFW121B4DNI'        : Begin & Index =  31 & Arr =  1 & End
	'DFW221B5DNI'        : Begin & Index =  32 & Arr =  1 & End
	'DFW321B6DNI'        : Begin & Index =  33 & Arr =  1 & End
	'DFW421B7DNI'        : Begin & Index =  34 & Arr =  1 & End
	'DFK121B8DNI'        : Begin & Index =  35 & Arr =  1 & End
	'DFV221B9DNI'        : Begin & Index =  36 & Arr =  1 & End
	'DFQ221B10DNI'       : Begin & Index =  37 & Arr =  1 & End
	'DPFLEDP10V'         : Begin & Index =  38 & Arr =  1 & End
	'DPPHSWCONVP9V'      : Begin & Index =  39 & Arr =  1 & End
	'DFQ122B1DNI'        : Begin & Index =  40 & Arr =  1 & End
	'DFV122B2DNI'        : Begin & Index =  41 & Arr =  1 & End
	'DFKA122B3DNI'       : Begin & Index =  42 & Arr =  1 & End
	'DFW122B4DNI'        : Begin & Index =  43 & Arr =  1 & End
	'DFW222B5DNI'        : Begin & Index =  44 & Arr =  1 & End
	'DFW322B6DNI'        : Begin & Index =  45 & Arr =  1 & End
	'DFW422B7DNI'        : Begin & Index =  46 & Arr =  1 & End
	'DFK122B8DNI'        : Begin & Index =  47 & Arr =  1 & End
	'DFV222B9DNI'        : Begin & Index =  48 & Arr =  1 & End
	'DFQ222B10DNI'       : Begin & Index =  49 & Arr =  1 & End
	'DPFM7_2V'           : Begin & Index =  50 & Arr =  1 & End
	'DPPHSWCONVM9V'      : Begin & Index =  51 & Arr =  1 & End
	'DRQ111B1DNI'        : Begin & Index =  52 & Arr =  1 & End
	'DRV111B2DNI'        : Begin & Index =  53 & Arr =  1 & End
	'DRKA111B3DNI'       : Begin & Index =  54 & Arr =  1 & End
	'DRW111B4DNI'        : Begin & Index =  55 & Arr =  1 & End
	'DRW211B5DNI'        : Begin & Index =  56 & Arr =  1 & End
	'DRW311B6DNI'        : Begin & Index =  57 & Arr =  1 & End
	'DRW411B7DNI'        : Begin & Index =  58 & Arr =  1 & End
	'DRK111B8DNI'        : Begin & Index =  59 & Arr =  1 & End
	'DRV211B9DNI'        : Begin & Index =  60 & Arr =  1 & End
	'DRQ211B10DNI'       : Begin & Index =  61 & Arr =  1 & End
	'DPFP7_2V'           : Begin & Index =  62 & Arr =  1 & End
	'DPHKM15V'           : Begin & Index =  63 & Arr =  1 & End
	'DRQ112B1DNI'        : Begin & Index =  64 & Arr =  1 & End
	'DRV112B2DNI'        : Begin & Index =  65 & Arr =  1 & End
	'DRKA112B3DNI'       : Begin & Index =  66 & Arr =  1 & End
	'DRW112B4DNI'        : Begin & Index =  67 & Arr =  1 & End
	'DRW212B5DNI'        : Begin & Index =  68 & Arr =  1 & End
	'DRW312B6DNI'        : Begin & Index =  69 & Arr =  1 & End
	'DRW412B7DNI'        : Begin & Index =  70 & Arr =  1 & End
	'DRK112B8DNI'        : Begin & Index =  71 & Arr =  1 & End
	'DRV212B9DNI'        : Begin & Index =  72 & Arr =  1 & End
	'DRQ212B10DNI'       : Begin & Index =  73 & Arr =  1 & End
	'DPFLDAP6_2V'        : Begin & Index =  74 & Arr =  1 & End
	'DPHKP15V'           : Begin & Index =  75 & Arr =  1 & End
	'DRQ121B1DNI'        : Begin & Index =  76 & Arr =  1 & End
	'DRV121B2DNI'        : Begin & Index =  77 & Arr =  1 & End
	'DRKA121B3DNI'       : Begin & Index =  78 & Arr =  1 & End
	'DRW121B4DNI'        : Begin & Index =  79 & Arr =  1 & End
	'DRW221B5DNI'        : Begin & Index =  80 & Arr =  1 & End
	'DRW321B6DNI'        : Begin & Index =  81 & Arr =  1 & End
	'DRW421B7DNI'        : Begin & Index =  82 & Arr =  1 & End
	'DRK121B8DNI'        : Begin & Index =  83 & Arr =  1 & End
	'DRV221B9DNI'        : Begin & Index =  84 & Arr =  1 & End
	'DRQ221B10DNI'       : Begin & Index =  85 & Arr =  1 & End
	'DPHKP5V'            : Begin & Index =  86 & Arr =  1 & End
	'DRQ122B1DNI'        : Begin & Index =  87 & Arr =  1 & End
	'DRV122B2DNI'        : Begin & Index =  88 & Arr =  1 & End
	'DRKA122B3DNI'       : Begin & Index =  89 & Arr =  1 & End
	'DRW122B4DNI'        : Begin & Index =  90 & Arr =  1 & End
	'DRW222B5DNI'        : Begin & Index =  91 & Arr =  1 & End
	'DRW322B6DNI'        : Begin & Index =  92 & Arr =  1 & End
	'DRW422B7DNI'        : Begin & Index =  93 & Arr =  1 & End
	'DRK122B8DNI'        : Begin & Index =  94 & Arr =  1 & End
	'DRV222B9DNI'        : Begin & Index =  95 & Arr =  1 & End
	'DRQ222B10DNI'       : Begin & Index =  96 & Arr =  1 & End
	'DPFLDBP6_2V'        : Begin & Index =  97 & Arr =  1 & End
;
;			First AEU Analog sweep array
;
	'DABD1V'             : Begin & Index =   1 & Arr =  2 & End
	'DARREF1BD1'         : Begin & Index =   2 & Arr =  2 & End
	'DAW323_4AMPT'       : Begin & Index =   3 & Arr =  2 & End
	'DAV113_4ADT'        : Begin & Index =   4 & Arr =  2 & End
	'DAW2_14_23AMP_ADT'  : Begin & Index =   5 & Arr =  2 & End
	'DAIHK1BDT'          : Begin & Index =   6 & Arr =  2 & End
	'DACONVBDT'          : Begin & Index =   7 & Arr =  2 & End
	'DPV221_2RXBT'       : Begin & Index =   8 & Arr =  2 & End
	'DPW111_2RXBT'       : Begin & Index =   9 & Arr =  2 & End
	'DPW321_2FPAT'       : Begin & Index =  10 & Arr =  2 & End
	'DRV222RXBAMPT'      : Begin & Index =  11 & Arr =  2 & End
	'DRW111RXBAMPT'      : Begin & Index =  12 & Arr =  2 & End
	'DRW221RXBAMPT'      : Begin & Index =  13 & Arr =  2 & End
	'DRKA12RXBRIBT'      : Begin & Index =  14 & Arr =  2 & End
	'DRQ2RXBRIBT'        : Begin & Index =  15 & Arr =  2 & End
	'DRPYPSHPRTKT'       : Begin & Index =  16 & Arr =  2 & End
	'DARREF2BD1'         : Begin & Index =  17 & Arr =  2 & End
	'DFK1AFEEDT'         : Begin & Index =  18 & Arr =  2 & End
	'DFQ1AFEEDT'         : Begin & Index =  19 & Arr =  2 & End
	'DFW3BFEEDT'         : Begin & Index =  20 & Arr =  2 & End
	'DFK1BOMTT'          : Begin & Index =  21 & Arr =  2 & End
	'DFQ1BOMTT'          : Begin & Index =  22 & Arr =  2 & End
	'DFW3AOMTT'          : Begin & Index =  23 & Arr =  2 & End
	'DFV22FPATEET'       : Begin & Index =  24 & Arr =  2 & End
	'DFW32FPATEET'       : Begin & Index =  25 & Arr =  2 & End
	'DTATOPSECT'         : Begin & Index =  26 & Arr =  2 & End
	'DTABOTSECT'         : Begin & Index =  27 & Arr =  2 & End
	'DTBMIDSECT'         : Begin & Index =  28 & Arr =  2 & End
	'DTBTOPPRIT'         : Begin & Index =  29 & Arr =  2 & End
	'DTAMIDPRIT'         : Begin & Index =  30 & Arr =  2 & End
	'DTAPXMIDRADT'       : Begin & Index =  31 & Arr =  2 & End
	'DTBMXBOTRADT'       : Begin & Index =  32 & Arr =  2 & End
	'DRK113RFBI0'        : Begin & Index =  33 & Arr =  2 & End
	'DRK114RFBI1'        : Begin & Index =  34 & Arr =  2 & End
	'DRK123RFBI2'        : Begin & Index =  35 & Arr =  2 & End
	'DRK124RFBI3'        : Begin & Index =  36 & Arr =  2 & End
	'DRW413RFBI8'        : Begin & Index =  37 & Arr =  2 & End
	'DRW414RFBI9'        : Begin & Index =  38 & Arr =  2 & End
	'DRW423RFBI10'       : Begin & Index =  39 & Arr =  2 & End
	'DRW424RFBI11'       : Begin & Index =  40 & Arr =  2 & End
	'DRW313RFBI16'       : Begin & Index =  41 & Arr =  2 & End
	'DRW314RFBI17'       : Begin & Index =  42 & Arr =  2 & End
	'DRW323RFBI18'       : Begin & Index =  43 & Arr =  2 & End
	'DRW324RFBI19'       : Begin & Index =  44 & Arr =  2 & End
	'DRW213RFBI24'       : Begin & Index =  45 & Arr =  2 & End
	'DRW214RFBI25'       : Begin & Index =  46 & Arr =  2 & End
	'DRW223RFBI26'       : Begin & Index =  47 & Arr =  2 & End
	'DRW224RFBI27'       : Begin & Index =  48 & Arr =  2 & End
	'DRV113RFBI32'       : Begin & Index =  49 & Arr =  2 & End
	'DRV114RFBI33'       : Begin & Index =  50 & Arr =  2 & End
	'DRV123RFBI34'       : Begin & Index =  51 & Arr =  2 & End
	'DRV124RFBI35'       : Begin & Index =  52 & Arr =  2 & End
	'DAP15VBD1'          : Begin & Index =  53 & Arr =  2 & End
	'DAM15VBD1'          : Begin & Index =  54 & Arr =  2 & End
	'DAP12VBD1'          : Begin & Index =  55 & Arr =  2 & End
	'DAM12VBD1'          : Begin & Index =  56 & Arr =  2 & End
	'DAP5VBD1'           : Begin & Index =  57 & Arr =  2 & End
;
;			Second AEU Analog sweep array
;
	'DABD2V'             : Begin & Index =   1 & Arr =  3 & End
	'DARREF1BD2'         : Begin & Index =   2 & Arr =  3 & End
	'DAW113_4ADT'        : Begin & Index =   3 & Arr =  3 & End
	'DAV223_4AMPT'       : Begin & Index =   4 & Arr =  3 & End
	'DAQ113_4ADT'        : Begin & Index =   5 & Arr =  3 & End
	'DAIHK2BDT'          : Begin & Index =   6 & Arr =  3 & End
	'DASPARE1'           : Begin & Index =   7 & Arr =  3 & End
	'DPV111_2FPAT'       : Begin & Index =   8 & Arr =  3 & End
	'DPW321_2RXBT'       : Begin & Index =   9 & Arr =  3 & End
	'DPW221_2FPAT'       : Begin & Index =  10 & Arr =  3 & End
	'DRV111RXBAMPT'      : Begin & Index =  11 & Arr =  3 & End
	'DRW321RXBAMPT'      : Begin & Index =  12 & Arr =  3 & End
	'DRK12RXBRIBT'       : Begin & Index =  13 & Arr =  3 & End
	'DRQ1RXBRIBT'        : Begin & Index =  14 & Arr =  3 & End
	'DRW3RXBRIBT'        : Begin & Index =  15 & Arr =  3 & End
	'DRMYPSHPRTKT'       : Begin & Index =  16 & Arr =  3 & End
	'DARREF2BD2'         : Begin & Index =  17 & Arr =  3 & End
	'DFKA1BFEEDT'        : Begin & Index =  18 & Arr =  3 & End
	'DFQ2BFEEDT'         : Begin & Index =  19 & Arr =  3 & End
	'DFW3AFEEDT'         : Begin & Index =  20 & Arr =  3 & End
	'DFKA1AOMTT'         : Begin & Index =  21 & Arr =  3 & End
	'DFQ2AOMTT'          : Begin & Index =  22 & Arr =  3 & End
	'DFW3BOMTT'          : Begin & Index =  23 & Arr =  3 & End
	'DFV11FPATEET'       : Begin & Index =  24 & Arr =  3 & End
	'DFW11FPATEET'       : Begin & Index =  25 & Arr =  3 & End
	'DFW22FPATEET'       : Begin & Index =  26 & Arr =  3 & End
	'DTAMIDSECT'         : Begin & Index =  27 & Arr =  3 & End
	'DTBTOPSECT'         : Begin & Index =  28 & Arr =  3 & End
	'DTATOPPRIT'         : Begin & Index =  29 & Arr =  3 & End
	'DTBMIDPRIT'         : Begin & Index =  30 & Arr =  3 & End
	'DTBPXMIDRADT'       : Begin & Index =  31 & Arr =  3 & End
	'DTAMXTOPRADT'       : Begin & Index =  32 & Arr =  3 & End
	'DRW113RFBI4'        : Begin & Index =  33 & Arr =  3 & End
	'DRW114RFBI5'        : Begin & Index =  34 & Arr =  3 & End
	'DRW123RFBI2'        : Begin & Index =  35 & Arr =  3 & End
	'DRW124RFBI3'        : Begin & Index =  36 & Arr =  3 & End
	'DRV213RFBI12'       : Begin & Index =  37 & Arr =  3 & End
	'DRV214RFBI13'       : Begin & Index =  38 & Arr =  3 & End
	'DRV223RFBI14'       : Begin & Index =  39 & Arr =  3 & End
	'DRV224RFBI15'       : Begin & Index =  40 & Arr =  3 & End
	'DRQ113RFBI20'       : Begin & Index =  41 & Arr =  3 & End
	'DRQ114RFBI21'       : Begin & Index =  42 & Arr =  3 & End
	'DRQ123RFBI22'       : Begin & Index =  43 & Arr =  3 & End
	'DRQ124RFBI23'       : Begin & Index =  44 & Arr =  3 & End
	'DRQ213RFBI28'       : Begin & Index =  45 & Arr =  3 & End
	'DRQ214RFBI29'       : Begin & Index =  46 & Arr =  3 & End
	'DRQ223RFBI30'       : Begin & Index =  47 & Arr =  3 & End
	'DRQ224RFBI31'       : Begin & Index =  48 & Arr =  3 & End
	'DRKA113RFBI36'      : Begin & Index =  49 & Arr =  3 & End
	'DRKA114RFBI37'      : Begin & Index =  50 & Arr =  3 & End
	'DRKA123RFBI38'      : Begin & Index =  51 & Arr =  3 & End
	'DRKA124RFBI39'      : Begin & Index =  52 & Arr =  3 & End
	'DAP15VBD2'          : Begin & Index =  53 & Arr =  3 & End
	'DAM15VBD2'          : Begin & Index =  54 & Arr =  3 & End
	'DAP12VBD2'          : Begin & Index =  55 & Arr =  3 & End
	'DAM12VBD2'          : Begin & Index =  56 & Arr =  3 & End
	'DAP5VBD2'           : Begin & Index =  57 & Arr =  3 & End
;
	Else : Begin & Index = 0 & Arr = 0 & status = 0 & End
EndCase
;
;                       Convert the index into an IDL-style array 
;                       index (starting at 0) and return.
;
Index = Index - 1
;
Return
End
