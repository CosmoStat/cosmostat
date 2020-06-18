;+
; NAME:
;	TVLUT
;
; PURPOSE:
;      Display an image with a zoom factor, and the LUT is displayed
;      The zoom factor is calculated automatically in order to visualize
;      the image in a window of size 320x320 (for an image 32x32, the
;      zoom factor is 320/32 = 10). An offset is automatically calculated
;      (or is set by keywords) in order to center the image in the IDL window. 
;
; CALLING SEQUENCE:
;	TVLUT, Image 
;
; INPUTS:
;	Image -- two dimensional array: image to visualize
;
; KEYED INPUTS: 
;       Depx -- scalar: offset in x (default is 50)
;       Depy -- scalar: offset in y (default is 80)
;
; EXAMPLE:
; tvlut, my_array
;
; 
; MODIFICATION HISTORY:
; 	Written by:	JL Starck 2/3/95 
;	modification of !order : HA 19/10/95
;       15-Mar-1996 SO, header upgrade
;
;-


PRO tvlut, IMAGE, Depx, Depy

if !window EQ -1 then window, xsize=480, ysize=480

vsize = size(IMAGE)
if vsize[0] NE 2 then BEGIN
     print, 'CALLING SEQUENCE: TVLUT, Image
     goto, CLOSING
     END

Nx = vsize[1]
Ny = vsize[2]

MINIMAGE = MIN(IMAGE)
MAXIMAGE = MAX(IMAGE)

S = 320
ech=S/float(max((size(IMAGE))[1:2]))
imageout=congrid(IMAGE,ech*(size(IMAGE))[1],ech*(size(IMAGE))[2])

LUT=FLOAT(INDGEN(1,S))*(MAXIMAGE - MINIMAGE)/S + MINIMAGE
LUT = REBIN(LUT, 30, S, /SAMPLE)

IF (!order EQ 1) THEN LUT[*,*] = LUT[*,S-INDGEN(S)]

ERASE

if not keyword_set(Depx) then Depx = 50
if not keyword_set(Depy) then Depy = 80

TVSCL, imageout, Depx, Depy

TVSCL, LUT, S+Depx+30, Depy

!P.COLOR=250

XYOUTS, S+Depx+30, Depy+4+s, STRING(STRCOMPRESS(MAXIMAGE)), /DEVICE,CHARSIZE=1.2
XYOUTS, S+Depx+30, Depy-10, STRING(STRCOMPRESS(MINIMAGE)), /DEVICE,CHARSIZE=1.2

CLOSING:

RETURN

END
