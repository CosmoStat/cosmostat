PRO MR1D_DEGLITCH, Cube, filter=filter, Nscale=Nscale, Nsigma=Nsigma, mask=mask
;+ 
; NAME: 
;       MR1D_DEGLITCH
;
; PURPOSE: 
;      Deglitch a cube using the multiresolution median.
;
; CALLING SEQUENCE: 
;   MR1D_DEGLITCH, Cube, Nscale=Nscale, Nsigma=Nsigma, mask=mask
;
; INPUT-OUTPUT: 
;   Cube -- IDL 3D array: Cube to degltich
;
;
; KEYED INPUTS: 
;   Nscale -- scalar: number of scale for the multiresolution analysis
;   Nsigma -- float: glitch detection at Nsigma
;   filter -- scalar : if set, then a adaptative temporal filter is applied.
;
; KEYED OUTPUTS: 
;   mask-- IDL 3D array: Mask where glitches have been detected
; 
; MODIFICATION HISTORY: 
;    13-Jan-1996 JL Starck written with template_gen 
;-
 
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'MR1D_DEGLITCH, Cube, filter=filter, Nscale=Nscale, Nsigma=Nsigma, mask=mask'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
vsize = size(cube)
Nx = vsize(1)
Ny = vsize(2)
Nz = vsize(3)
mask = bytarr(Nx,Ny,Nz)

; if not keyword_set(SigmaNoise) then Sigma = get_noise(Cube)
Opt = ' '
if keyword_set(filter) then Opt = Opt + '-f '
if keyword_set(Nscale) then $
            Opt = Opt + ' -n '+  strcompress(string(Nscale),/REMOVE_ALL)
if keyword_set(Nsigma) then $
            Opt = Opt + ' -s '+  strcompress(string(Nsigma),/REMOVE_ALL)

NameCube = 'xx_cube.fits'
NameResult = 'xx_result.fits'

writefits, NameCube, Cube
com = 'mr1d_deglitch ' + OPT + ' ' + NameCube + ' ' + NameResult
spawn, com
OldCube = Cube
Cube = readfits(NameResult, /silent) 

delete, NameCube
delete, NameResult

if not keyword_set(filter) then BEGIN
   ind = where ( abs(OldCube - Cube) gt 1e-7, count)
   if count GT 0 then mask(ind) = 1
   END ELSE BEGIN
   if not keyword_set(Nsigma) then Nsigma = 3.
   if not keyword_set(SigmaNoise) then Sigma = Nsigma*get_noise(OldCube)
   ind = where ( abs(OldCube - Cube) gt Sigma, count)
   if count GT 0 then mask(ind) = 1
   END

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN
 
 END
