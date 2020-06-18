; $Id: loadct_012.pro,v 1.1 1998/10/02 15:37:44 greason Exp $

PRO LOADCT_012, Color_Table

; Procedure to load a standard IDL color table with the first three colors
; reserved as:
;   0 - white
;   1 - black
;   2 - grey

; Load the user-selected color table
COMMON Colors, R,G,B,R_Curr,G_Curr,B_Curr
LOADCT, Color_Table
; Define 0 to be white and 1 to be black and 2 to be grey
R[0] = 255 & G[0] = 255 & B[0] = 255
R[1] = 0   & G[1] = 0   & B[1] = 0
R[2] = 175 & G[2] = 175 & B[2] = 175
; Re-load the color table
TVLCT, R,G,B

RETURN
END

