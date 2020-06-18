Function Sci_GetMnemonic, Sci, Mnem, Status, NFRAMES=Nframes, _EXTRA=Extra
;+
; NAME:
;       Sci_GetMnemonic
; PURPOSE:
;       Returns the physical value associated with a science mnemonic,
;       extracting the data from a science telemetry packet.
; CALLING SEQUENCE:
;       Data = Sci_GetMnemonic(Sci, Mnem [,Status] [,Nframes= ] 
;             
; INPUTS:
;       Sci     The packet of science telemetry data.
;       Mnem    The name of the science mnemonic.  The mnemonic naming
;               conventions are as follows: let 'chan' be a radiometer 
;               channel string, such as 'K113', etc.  then:
;               'DK113'  returns all observations in time-order.
;               'QDK113' returns the first observation in each packet
;               'ADK113' returns the major frame (packet) average
;               'RDK113' returns the major frame (packet) rms
; OPTIONAL INPUT KEYWORD:
;       Nframes - Optional keyword which may be used to specify the number
;               of major frames over which to average detector data
;               or produce an rms of detector data. Default = 1 MF.
; OUTPUTS:
;       Status - An optional status value: 0=success, 1=undefined mnemonic.
; RETURNED:
;       Data   - The raw or calibrated/corrected telemetry value.
; EXAMPLE:
;       Return the first observation in each packet from a time-ordered
;       data file
;
;       IDL> fits_read_tod,'MAP_tod_20022162357_20022172357.fits',tod
;       IDL> data = sci_getmnemonic(tod.sci,'QDK113')
; MODIFICATION HISTORY:
;       Gary Hinshaw, NASA/GSFC, 29 December 1998
;       J. Weiland, 8 June 1999:  added Nframes keyword
;       BASELINE in pars structure.  RSH, 10 June 2002.
;       Remove calibration keywords  W. Landsman January 2003
;-
;ON_ERROR, 2
;
IF( N_PARAMS() LT 2 )THEN MESSAGE, $
  'Syntax:  Data = Sci_GetMnemonic(Sci, Mnem [,Status] [,Nframes= ] )' 

Status = 0
N_sci = N_ELEMENTS(Sci)
Mnemonic = STRUPCASE(STRTRIM(Mnem,2))
Mnem_Type = STRMID(Mnemonic,0,1)

IF KEYWORD_SET(Nframes) THEN Nframes = LONG(Nframes) ELSE  Nframes = 1L


CASE Mnem_Type OF
 'D' : Channel = STRMID(Mnemonic,1,10)
 ELSE :Channel = STRMID(Mnemonic,2,10)
ENDCASE
; Extract the band and index from the supplied channel
Extract_Band_Index, Channel, Band, Index

; Extract the data for the requested channel
CASE Mnem_Type OF
 'Q'  : BEGIN
        CASE Band OF
         'K'  : Data = REFORM(Sci.K[Index,0])
         'KA' : Data = REFORM(Sci.Ka[Index,0])
         'Q'  : Data = REFORM(Sci.Q[Index,0])
         'V'  : Data = REFORM(Sci.V[Index,0])
         'W'  : Data = REFORM(Sci.W[Index,0])
        ENDCASE
        END
 ELSE : BEGIN
        CASE Band OF
         'K'  : Data = reform(Sci.K[Index,*])
         'KA' : Data = reform(Sci.Ka[Index,*])
         'Q'  : Data = reform(Sci.Q[Index,*])
         'V'  : Data = reform(Sci.V[Index,*])
         'W'  : Data = reform(Sci.W[Index,*])
        ENDCASE
        END
ENDCASE

; Complete the requested data processing
CASE Mnem_Type OF
 'D' : Data = REFORM(Data,N_ELEMENTS(Data))
 'Q' : Data = REFORM(Data,N_Sci)
'A' : Data = Frame_Avg(Data,Nframes=Nframes)
'R' : Data = Frame_RMS(Data,Nframes=Nframes)
ENDCASE

RETURN, Data
END
