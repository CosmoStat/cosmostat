Pro AIHK_Mnem2Serial_List, info
;+
; NAME:
;    AIHK_Mnem2Serial_List
; PURPOSE:
;    Returns an array of structures that associates mnemonics
;    with serial numbers for sensor ids.
; CALLING SEQUENCE:
;    AIHK_Mnem2Serial_List, info
; OUTPUTS:
;    info - The array of structures.
; MODIFICATION HISTORY:
;    Automatically generated.
;    Manual edit for A side secondary prt id changes, JW, 22 Feb 2000.
;    Manual edit for B side secondary prt id changes, JW, 08 Mar 2000.
;    Manual edit of prt serial nos. for W3BOMT & W3AFEED, JW, 24 May 2000.
;-
on_error, 2
;
;            Check arguments.
;
If (n_params() LT 1) Then message, $
     'Syntax: AIHK_Mnem2Serial_List, info'
;
;            Define the array of structures.
;
info = replicate({Mnemonic:'', SerialNum:'', SensorId:''}, 64)
;
;            And fill the array.
;
i = -1
i = i + 1
info[i].Mnemonic  = 'DABD1V'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK1_AC_0'
;
i = i + 1
info[i].Mnemonic  = 'DARREF1BD1'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK1_AC_1'
;
i = i + 1
info[i].Mnemonic  = 'DAW323_4AMPT'
info[i].SerialNum = 'UG85'
info[i].SensorId  = 'IHK1_AC_2'
;
i = i + 1
info[i].Mnemonic  = 'DAV113_4ADT'
info[i].SerialNum = 'UH25'
info[i].SensorId  = 'IHK1_AC_3'
;
i = i + 1
info[i].Mnemonic  = 'DAW2_14_23AMP_ADT'
info[i].SerialNum = 'UH26'
info[i].SensorId  = 'IHK1_AC_4'
;
i = i + 1
info[i].Mnemonic  = 'DAIHK1BDT'
info[i].SerialNum = 'UG55'
info[i].SensorId  = 'IHK1_AC_5'
;
i = i + 1
info[i].Mnemonic  = 'DACONVBDT'
info[i].SerialNum = 'UH96'
info[i].SensorId  = 'IHK1_AC_6'
;
i = i + 1
info[i].Mnemonic  = 'DPV221_2RXBT'
info[i].SerialNum = 'UH90'
info[i].SensorId  = 'IHK1_AC_7'
;
i = i + 1
info[i].Mnemonic  = 'DPW111_2RXBT'
info[i].SerialNum = 'UG59'
info[i].SensorId  = 'IHK1_AC_8'
;
i = i + 1
info[i].Mnemonic  = 'DPW321_2FPAT'
info[i].SerialNum = 'UH95'
info[i].SensorId  = 'IHK1_AC_9'
;
i = i + 1
info[i].Mnemonic  = 'DRV222RXBAMPT'
info[i].SerialNum = 'UH03'
info[i].SensorId  = 'IHK1_AC_10'
;
i = i + 1
info[i].Mnemonic  = 'DRW111RXBAMPT'
info[i].SerialNum = 'UH13'
info[i].SensorId  = 'IHK1_AC_11'
;
i = i + 1
info[i].Mnemonic  = 'DRW221RXBAMPT'
info[i].SerialNum = 'UG93'
info[i].SensorId  = 'IHK1_AC_12'
;
i = i + 1
info[i].Mnemonic  = 'DRKA12RXBRIBT'
info[i].SerialNum = 'UG64'
info[i].SensorId  = 'IHK1_AC_13'
;
i = i + 1
info[i].Mnemonic  = 'DRQ2RXBRIBT'
info[i].SerialNum = 'UH24'
info[i].SensorId  = 'IHK1_AC_14'
;
i = i + 1
info[i].Mnemonic  = 'DRPYPSHPRTKT'
info[i].SerialNum = 'UG74'
info[i].SensorId  = 'IHK1_AC_15'
;
i = i + 1
info[i].Mnemonic  = 'DARREF2BD1'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK1_AC_16'
;
i = i + 1
info[i].Mnemonic  = 'DFK1AFEEDT'
info[i].SerialNum = 'UG89'
info[i].SensorId  = 'IHK1_AC_17'
;
i = i + 1
info[i].Mnemonic  = 'DFQ1AFEEDT'
info[i].SerialNum = 'UG81'
info[i].SensorId  = 'IHK1_AC_18'
;
i = i + 1
info[i].Mnemonic  = 'DFW3BFEEDT'
info[i].SerialNum = 'UG60'
info[i].SensorId  = 'IHK1_AC_19'
;
i = i + 1
info[i].Mnemonic  = 'DFK1BOMTT'
info[i].SerialNum = 'UG77'
info[i].SensorId  = 'IHK1_AC_20'
;
i = i + 1
info[i].Mnemonic  = 'DFQ1BOMTT'
info[i].SerialNum = 'UG82'
info[i].SensorId  = 'IHK1_AC_21'
;
i = i + 1
info[i].Mnemonic  = 'DFW3AOMTT'
info[i].SerialNum = 'UG54'
info[i].SensorId  = 'IHK1_AC_22'
;
i = i + 1
info[i].Mnemonic  = 'DFV22FPATEET'
info[i].SerialNum = 'UH23'
info[i].SensorId  = 'IHK1_AC_23'
;
i = i + 1
info[i].Mnemonic  = 'DFW32FPATEET'
info[i].SerialNum = 'UG65'
info[i].SensorId  = 'IHK1_AC_24'
;
i = i + 1
info[i].Mnemonic  = 'DTATOPSECT'
info[i].SerialNum = 'UH21'         ; was 'UG95' prior to 19 Feb 2000
info[i].SensorId  = 'IHK1_AC_25'
;
i = i + 1
info[i].Mnemonic  = 'DTABOTSECT'
info[i].SerialNum = 'UG57'         ; was 'UH02' prior to 19 Feb 2000
info[i].SensorId  = 'IHK1_AC_26'
;
i = i + 1
info[i].Mnemonic  = 'DTBMIDSECT'
info[i].SerialNum = 'UG94'         ; was 'UK17' prior to 05 Mar 2000
info[i].SensorId  = 'IHK1_AC_27'
;
i = i + 1
info[i].Mnemonic  = 'DTBTOPPRIT'
info[i].SerialNum = 'UG90'
info[i].SensorId  = 'IHK1_AC_28'
;
i = i + 1
info[i].Mnemonic  = 'DTAMIDPRIT'
info[i].SerialNum = 'UG79'
info[i].SensorId  = 'IHK1_AC_29'
;
i = i + 1
info[i].Mnemonic  = 'DTAPXMIDRADT'
info[i].SerialNum = 'UK19'
info[i].SensorId  = 'IHK1_AC_30'
;
i = i + 1
info[i].Mnemonic  = 'DTBMXBOTRADT'
info[i].SerialNum = 'UG61'
info[i].SensorId  = 'IHK1_AC_31'
;
i = i + 1
info[i].Mnemonic  = 'DABD2V'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK2_AC_0'
;
i = i + 1
info[i].Mnemonic  = 'DARREF1BD2'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK2_AC_1'
;
i = i + 1
info[i].Mnemonic  = 'DAW113_4ADT'
info[i].SerialNum = 'UH07'
info[i].SensorId  = 'IHK2_AC_2'
;
i = i + 1
info[i].Mnemonic  = 'DAV223_4AMPT'
info[i].SerialNum = 'UG78'
info[i].SensorId  = 'IHK2_AC_3'
;
i = i + 1
info[i].Mnemonic  = 'DAQ113_4ADT'
info[i].SerialNum = 'UI88'
info[i].SensorId  = 'IHK2_AC_4'
;
i = i + 1
info[i].Mnemonic  = 'DAIHK2BDT'
info[i].SerialNum = 'UG62'
info[i].SensorId  = 'IHK2_AC_5'
;
i = i + 1
info[i].Mnemonic  = 'DASPARE1'
info[i].SerialNum = 'not used'
info[i].SensorId  = 'IHK2_AC_6'
;
i = i + 1
info[i].Mnemonic  = 'DPV111_2FPAT'
info[i].SerialNum = 'UG83'
info[i].SensorId  = 'IHK2_AC_7'
;
i = i + 1
info[i].Mnemonic  = 'DPW321_2RXBT'
info[i].SerialNum = 'UG66'
info[i].SensorId  = 'IHK2_AC_8'
;
i = i + 1
info[i].Mnemonic  = 'DPW221_2FPAT'
info[i].SerialNum = 'UG70'
info[i].SensorId  = 'IHK2_AC_9'
;
i = i + 1
info[i].Mnemonic  = 'DRV111RXBAMPT'
info[i].SerialNum = 'UG97'
info[i].SensorId  = 'IHK2_AC_10'
;
i = i + 1
info[i].Mnemonic  = 'DRW321RXBAMPT'
info[i].SerialNum = 'UG84'
info[i].SensorId  = 'IHK2_AC_11'
;
i = i + 1
info[i].Mnemonic  = 'DRK12RXBRIBT'
info[i].SerialNum = 'UG86'
info[i].SensorId  = 'IHK2_AC_12'
;
i = i + 1
info[i].Mnemonic  = 'DRQ1RXBRIBT'
info[i].SerialNum = 'UG92'
info[i].SensorId  = 'IHK2_AC_13'
;
i = i + 1
info[i].Mnemonic  = 'DRW3RXBRIBT'
info[i].SerialNum = 'UG56'
info[i].SensorId  = 'IHK2_AC_14'
;
i = i + 1
info[i].Mnemonic  = 'DRMYPSHPRTKT'
info[i].SerialNum = 'UH97'
info[i].SensorId  = 'IHK2_AC_15'
;
i = i + 1
info[i].Mnemonic  = 'DARREF2BD2'
info[i].SerialNum = 'n.a.'
info[i].SensorId  = 'IHK2_AC_16'
;
i = i + 1
info[i].Mnemonic  = 'DFKA1BFEEDT'
info[i].SerialNum = 'UH14'
info[i].SensorId  = 'IHK2_AC_17'
;
i = i + 1
info[i].Mnemonic  = 'DFQ2BFEEDT'
info[i].SerialNum = 'UK18'
info[i].SensorId  = 'IHK2_AC_18'
;
i = i + 1
info[i].Mnemonic  = 'DFW3AFEEDT'
info[i].SerialNum = 'UH22'         ; was 'UG71' prior to 24 May 2000
info[i].SensorId  = 'IHK2_AC_19'
;
i = i + 1
info[i].Mnemonic  = 'DFKA1AOMTT'
info[i].SerialNum = 'UH15'
info[i].SensorId  = 'IHK2_AC_20'
;
i = i + 1
info[i].Mnemonic  = 'DFQ2AOMTT'
info[i].SerialNum = 'UG96'
info[i].SensorId  = 'IHK2_AC_21'
;
i = i + 1
info[i].Mnemonic  = 'DFW3BOMTT'
info[i].SerialNum = 'UH94'         ; was 'UH06' prior to 24 May 2000
info[i].SensorId  = 'IHK2_AC_22'
;
i = i + 1
info[i].Mnemonic  = 'DFV11FPATEET'
info[i].SerialNum = 'UG91'
info[i].SensorId  = 'IHK2_AC_23'
;
i = i + 1
info[i].Mnemonic  = 'DFW11FPATEET'
info[i].SerialNum = 'UG73'
info[i].SensorId  = 'IHK2_AC_24'
;
i = i + 1
info[i].Mnemonic  = 'DFW22FPATEET'
info[i].SerialNum = 'UH28'
info[i].SensorId  = 'IHK2_AC_25'
;
i = i + 1
info[i].Mnemonic  = 'DTAMIDSECT'
info[i].SerialNum = 'UG63'         ; was 'UH00' prior to 19 Feb 2000
info[i].SensorId  = 'IHK2_AC_26'
;
i = i + 1
info[i].Mnemonic  = 'DTBTOPSECT'
info[i].SerialNum = 'UG75'         ; was 'UG76' prior to 05 Mar 2000
info[i].SensorId  = 'IHK2_AC_27'
;
i = i + 1
info[i].Mnemonic  = 'DTATOPPRIT'
info[i].SerialNum = 'UG99'
info[i].SensorId  = 'IHK2_AC_28'
;
i = i + 1
info[i].Mnemonic  = 'DTBMIDPRIT'
info[i].SerialNum = 'UH18'
info[i].SensorId  = 'IHK2_AC_29'
;
i = i + 1
info[i].Mnemonic  = 'DTBPXMIDRADT'
info[i].SerialNum = 'UG67'
info[i].SensorId  = 'IHK2_AC_30'
;
i = i + 1
info[i].Mnemonic  = 'DTAMXTOPRADT'
info[i].SerialNum = 'UH05'
info[i].SensorId  = 'IHK2_AC_31'
;
;           Done.
;
Return
End
