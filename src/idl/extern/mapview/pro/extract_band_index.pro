PRO Extract_Band_Index, Channel, Band, Index, Status
;+
; NAME:
;    EXTRACT_BAND_INDEX
; PURPOSE:
;    Extract the frequency band ('K','KA','Q','V','W') and the index
;    (an integer from 0-15) for use in extracting science data from IDL data
;    structures.
; CALLING SEQUENCE:
;    Extract_Band_Index, Channel, Band, Index, Status 
; INPUT:
;  Channel - A scalar string containing the radiometer channel to be culled,
;		eg., 'KA113', must be either 4 or 5 characters in length
; OUTPUTS:
;  Band - A string containing the frequency band ('K','KA','Q','V','W').
;  Index -  An integer containing the array index of the specific channel,
;           from 0-3 for K, Ka; from 0-7 for Q, V; from 0-15 for W.
;  Status - return status flag: 0=valid channel, -1=undefined channel.
; EXAMPLE:
;     extract_band_index,'V114',band,index,status
;     ===> band = 'V'  index = 1    status= 0
; REVISION HISTORY:
;     Old    January 2003
; -
 Status = -1
ON_IOError, Abort
Channel = STRTRIM(STRUPCASE(Channel),2)
CASE  STRLEN(Channel) OF
 4 : BEGIN
      Band = STRMID(Channel,0,1)
      Ida = LONG(STRMID(Channel,1,1))
      Irad = LONG(STRMID(Channel,2,1))
      Ich = LONG(STRMID(Channel,3,1))
     END
 5 : BEGIN
      Band = STRMID(Channel,0,2)
      Ida = LONG(STRMID(Channel,2,1))
      Irad = LONG(STRMID(Channel,3,1))
      Ich = LONG(STRMID(Channel,4,1))
     END
 ELSE : RETURN
ENDCASE

; Do some error checking
CASE Band OF
 'K'  : IF( Ida NE 1 )THEN RETURN
 'KA' : IF( Ida NE 1 )THEN RETURN
 'Q'  : IF( Ida LT 1 OR Ida GT 2 )THEN RETURN
 'V'  : IF( Ida LT 1 OR Ida GT 2 )THEN RETURN
 'W'  : IF( Ida LT 1 OR Ida GT 4 )THEN RETURN
 ELSE : RETURN
ENDCASE
IF( Irad LT 1 OR Irad GT 2 )THEN RETURN
IF( Ich LT 3 OR Ich GT 4 )THEN RETURN

Index = 4L*(Ida-1) + 2L*(Irad-1) + (Ich-3)

Status = 0

Abort:
RETURN
END
