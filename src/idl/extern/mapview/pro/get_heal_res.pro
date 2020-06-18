FUNCTION GET_HEAL_RES, N_Pixels
;+
; NAME:
;       GET_HEAL_RES()
; PURPOSE:
;       Compute the resolution of a heal-pixelized sky map
; EXPLANATION:
;       Function to compute the resolution of a heal-pixelized sky map, given 
;       the number of map pixels.
; CALLING SEQUENCE
;       res = GET_HEAL_RES( N_Pixels)
; INPUTS:
;       N_pixels - integer scalar giving the number of pixels in the sky map
;
; OUTPUTS:
;       res - Integer scalar (1-13) giving the resolution factor 
; REVISION HISTORY:
;      Written M. Greason  October 1998
; -

CASE N_Pixels OF
        48L : Resolution = 1
       192L : Resolution = 2
       768L : Resolution = 3
      3072L : Resolution = 4
     12288L : Resolution = 5
     49152L : Resolution = 6
    196608L : Resolution = 7
    786432L : Resolution = 8
   3145728L : Resolution = 9
  12582912L : Resolution = 10
  50331648L : Resolution = 11
 201326592L : Resolution = 12
 805306368L : Resolution = 13
 ELSE       : BEGIN
               PRINT,' We do not support the selected resolution...'
               PRINT,' Number of input pixels:', N_Pixels
               Resolution = -1
               RETURN, Resolution
              END
ENDCASE

RETURN, Resolution
END
