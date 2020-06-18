FUNCTION Pckt2Mnemonic, Pckt, Mnemonic, Status, TIME=Time, TJ=Tj, $
         NFRAMES=Nframes, _EXTRA=Extra
;+
; NAME:
;       Pckt2Mnemonic
; PURPOSE:
;       Procedure to extract data for a given mnemonic from a packet.  This 
;       is a wrapper to the various pckt2mnemonic routines available for 
;       each packet type.
; CALLING SEQUENCE:
;       data = Pckt2Mnemonic(Pckt, Mnemonic [, Status])
; INPUTS:
;       Pckt     - IDL structure giving Time-ordered data, e.g. as read by
;                  FITS_READ_TOD
;       Mnemonic - scalar string giving the mnemonic to extract, e.g. 'DK113'
; OUTPUTS:
;       Status   - A status code: 0=success, -1=error.
; RETURNED:
;       data     - The extracted data.
; OPTIONAL INPUT KEYWORDS:
;       /TIME    - If present and nonzero, the data returned as DLBARR[N,2]
;                 where data[*,0] gives the MAP reduced Julian date, and 
;                 data[*,1] gives the extracted data.   It is usually preferable
;                 to use the TJ output keyword to obtain the Julian dates. 
;       NFRAMES - The number of major frames over which to average detector data
;                 or produce an rms of detector data. Only applies to science
;                 data mnemonics beginning with 'A' (average) ro 'R' (rms). 
;                 Default = 1 MF.
;       Other keywords are used by the various data extraction routines.
; OPTIONAL OUTPUT KEYWORD:
;       TJ - vector of MAP-reduced Julian dates associated with each data 
;            point
; COMMENTS:
;       This routine acts as a wrapper for the data extraction routines (either
;       AIHK_Pckt2Mnemonic(), DIHK_Pckt2Mnemonic() or Sci_GetMnemonic()).  The 
;       appropriate time-ordered data structure  must be supplied.
; PROCEDURES CALLED:
;       AIHK_MnemTimeStamp, Extract_Band_Index, Mnem_Type, TimeTransform
;       AIHK_Pckt2Mnemonic(), DIHK_Pckt2Mnemonic(), Sci_GetMnemonic()
; MODIFICATION HISTORY:
;       Written by Gary Hinshaw, GSFC.
;       Extended the documentation.  MRG, RITSS, 03 March 2000.
;       SCFT Tj computation changed for supercom/pseudo tm compatibility. 
;            08 May 2001, JPW.
;       Remove spacecraft and GPIB mnemonics W. Landsman  January 2003
;-
ON_ERROR, 2
;
; Check arguments.
;
Status = -1
If (n_params() LT 2) Then message, $
   'Syntax:  data = Pckt2Mnemonic(Pckt, Mnemonic [, Status, TJ=, NFRAMES=])'
;
N_pckts = N_ELEMENTS(Pckt)
If (N_pckts LE 0) Then message, 'Data must be supplied.'
;
; Extract the raw data from the supplied packets
;
M_Type = Mnem_Type(Mnemonic)
M_Type = M_Type[0]              ;ensure scalar, not array
IF( M_Type EQ 'SCI' )THEN S_Type = STRUPCASE(STRMID(Mnemonic,0,1))
CASE M_Type OF
 'AIHK' : Data = AIHK_Pckt2Mnemonic(Pckt.aeu, Mnemonic, -1, $
                                    Status, _Extra=Extra)
 'DIHK' : Data = DIHK_Pckt2Mnemonic(Pckt.deu, Mnemonic, Status, _Extra=Extra)
 'SCI'  : Data = Sci_GetMnemonic(Pckt.sci, Mnemonic, Status,  $
                                  Nframes=Nframes, _Extra=Extra)
 'UNDEFINED' : RETURN, Status
ENDCASE
;
; Extract the julian time, if requested
;
 if keyword_set(time) or arg_present(Tj) then begin

CASE M_Type OF
 'AIHK' : BEGIN
           AIHK_MnemTimeStamp, Pckt.aeu, Mnemonic, TS
           TS = REFORM(TS,2,2*N_pckts)
           TimeTransform, TS, Tj, /TS2Jul
          END
 'DIHK' : TimeTransform, Pckt.deu, Tj, /Tel2Jul
 'SCI'  : BEGIN
           TimeTransform, Pckt.sci, Tj, /Tel2Jul
           Tj = reform(Tj, N_elements(Tj))
           IF( S_Type EQ 'D' )THEN BEGIN
             Channel = STRMID(Mnemonic,1,10)
             Extract_Band_Index, Channel, Band       
             Tj = Spread_Pckt_Tjul(Tj,Band)
           ENDIF
           IF( KEYWORD_SET(Nframes) )THEN BEGIN
             Index = Nframes*LINDGEN(N_ELEMENTS(Data))
             Tj = Tj[Index]
           ENDIF
          END
ENDCASE
 endif

; Maintain compatibilty with older calling sequence
IF  KEYWORD_SET(Time)  THEN Data = [[Tj],[Data]]
;
Status = 0
;
RETURN, Data
END
