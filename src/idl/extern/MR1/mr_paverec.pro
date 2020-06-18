FUNCTION MR_PAVEREC, Pave
;+ 
; NAME: 
;      MR_PAVEREC
;
; PURPOSE: 
;      reconstruct an image from its multiresolution transform
;
; CALLING SEQUENCE: 
;      output=MR_PAVEREC(Pave)
;
; INPUTS: 
;   Pave -- 3D IDL array: multiresolution transform
;
; OUTPUTS: 
;    output -- 2d IDL array: reconstructed image
;
; MODIFICATION HISTORY: 
;    19-Jan-1996 JL Starck  
;-
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=MR_PAVEREC(Pave)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

Imag=-1
vsize = size(Pave)
if vsize(0) NE 3 then BEGIN
    print, 'Error: first parameter has to be a cube'
    goto, CLOSING
    END

Nx = vsize(1)
Ny = vsize(2)
Np = vsize(3)
Imag = fltarr(Nx, Ny)
Imag(*,*) = Pave(*,*,Np-1)
for i = 0,Np-2 do Imag(*,*) = Imag(*,*) + Pave(*,*,i)
 
;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, Imag
 
 END
