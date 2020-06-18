FUNCTION FORMAT_IMAGE, NameImage
;+ 
; NAME: 
;        FORMAT_IMAGE
;
; PURPOSE: 
;        return the format of an image.
;        format can be:
;             ".d" for raw images
;             ".fits" for fits images
;             ".bdf" for midas images
;             ".cub" for isocam SAD images
;
; CALLING SEQUENCE: 
;   output=FORMAT_IMAGE(NameImage)
;
; INPUTS: 
;   NameImage -- string : image file name with extension
; 
;
; MODIFICATION HISTORY: 
;    13-Dec-1995 JL Starck
;-
  
;------------------------------------------------------------
; initialization
;------------------------------------------------------------
 
Format = 'unknown'

;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'output=FORMAT_IMAGE(NameImage)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

if rstrpos(NameImage, '.d') GT 0  then Format = 'raw' else $
if rstrpos(NameImage, '.fit') GT 0 then Format = 'fits' else $
if rstrpos(NameImage, '.bdf') GT 0 then Format = 'midas' else $
if rstrpos(NameImage, '.cub') GT 0 then Format = 'sad' 

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, Format
 
END
