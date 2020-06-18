FUNCTION MR_BACKGROUND, Image, nscale=nscale, border=border
;+ 
; NAME: 
;   MR_BACKGROUND
;
; PURPOSE: 
;   Estimate the background of an image
;
; CALLING SEQUENCE: 
;   output=MR_BACKGROUND(Image, nscale=nscale, border=border)
;
; INPUTS: 
;   Image -- IDL array: background image
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED INPUTS: 
;   nscale -- scalar: if set, fix the number of scales used for
;                     for the multiresolution transform
;   border -- scalar: if set, then the background is computed from the 
;                     border of the image.
;
; OUTPUTS: 
;    output -- IDl 2D array: border image
;
; MODIFICATION HISTORY: 
;    19-Jan-1996 JL Starck written with template_gen 
;-
 
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
Result = -1

 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'output=MR_BACKGROUND(Image, nscale=nscale, border=border)'
   GOTO, CLOSING
 ENDIF
 
vsize = size(image)
if vsize(0) NE 2 then BEGIN
   print, 'Error: first parameter must be an image'
   GOTO, CLOSING
END

;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
Nx = (size(image))(1)
Ny = (size(image))(2)
N = min([Nx,Ny])
if keyword_set(border) then $
BEGIN
   Bordup = image[1:Nx-2, 1:4]
   Borddown = image[1:Nx-2, Ny-5:Ny-2]
   Bordleft = image[1:4, 5:Ny-6]
   Bordright = image[Ny-5:Ny-2, 5:Ny-6]
   total = total(Bordup) + total(Borddown) + $
                                  total(Bordleft) + total(Bordright)
   Nel = N_ELEMENTS(Bordup) + N_ELEMENTS(Borddown) + $
                       N_ELEMENTS(Bordleft) + N_ELEMENTS(Bordright)
   Mean = total / Nel
   Result = fltarr(Nx,Ny)
   Result[*,*] = Mean
END ELSE BEGIN
   if not keyword_set(nscale) then nscale = fix((alog(N/16.)/alog(2.))+1.5)
   if nscale LT 3 then nscale = 3

   NameImag = 'xx_imag.fits'
   NameResult = 'xx_result.fits'
   NameBackgr = 'xx_bgr.fits'

   writefits, NameImag, image
   OPT = '-n '+ strcompress(string(nscale),/REMOVE_ALL) + ' -w '+NameBackgr
   com = 'mr_background ' + OPT + ' ' + NameImag + ' ' + NameResult
   spawn, com
   Result = readfits(NameBackgr, /silent)

   delete, NameImag
   delete, NameResult
   delete, NameBackgr 
   END
 
;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, Result
 
 END
