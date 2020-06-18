;PRO IM_SMOOTH, Im_in, Im_Out, method=method, Step_Trou=Step_Trou, 
;               Border=Border, WinSize=WinSize
;+ 
; NAME: 
;     IM_SMOOTH
;
; PURPOSE: 
;     Smooth an image, taking into account the border
;
; CATEGORY: 
;     III-3
;
; CALLING SEQUENCE: 
;   IM_SMOOTH, Im_in, Im_Out, method=method, Step_Trou=Step_Trou, Border=Border
;
; INPUTS: 
;   Im_in -- 2D IDL array: image to smooth
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED INPUTS: 
;   method -- string: smoothing method ('linear','bspline',or 'median')
;   Step_Trou -- scalar: 2^Step_Trou = distance between two adjancent 
;                                      pixels (default is 0)
;   Border -- scalar: type of border ('cont', 'mirror', 'zero', 'nobord')
;   WinSize -- scalar: window size (only used if method = 'median')
;
; OUTPUTS: 
;   Im_Out -- 2D IDL array: smoothed image 
; 
; OPTIONAL OUTPUT PARAMETERS: 
;   none
;
; MODIFICATION HISTORY: 
;    20-Dec-1995 JL Starck written with template_gen 
;-
 
;==============================
function test_cont, ind, N
ret = ind
if ind LT 0 then ret = 0 $
else if ind GE N then ret = N-1 
return, ret
end
;==============================
function test_mirror, ind, N
ret = ind
if ind LT 0 then ret = -ind $
else if ind GE N then ret = 2*N-ind-2
return, ret
end
;==============================

pro smooth_bspline, Im_in, Im_Out, Step_trou=Step_trou, Border=Border
   Nl = (size(Im_in))[1]
   Nc = (size(Im_in))[2]

   Step = 2^Step_trou
   Window_Size = Step*4+1
   K1D = fltarr(Window_Size)

   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.

   W2 = Window_Size/2
   K1D[W2] = C3
   K1D[W2-Step] = C2
   K1D[W2+Step] = C2
   K1D[W2-2*Step] = C1
   K1D[W2+2*Step] = C1
   
   Kernel = K1D # K1D
   Im_Out = convol(Im_in, Kernel, /EDGE_WRAP)
   ; border manadgement

END
;==============================

pro smooth_linear, Im_in, Im_Out, Step_trou
   Nl = (size(Im_in))[1]
   Nc = (size(Im_in))[2]
   Step = 2^Step_trou
   Window_Size = Step*2+1
   K1D = fltarr(Window_Size)

   W2 = Window_Size/2
   K1D[W2] = 0.5
   K1D[W2-Step] = 0.25
   K1D[W2+Step] = 0.25
   Kernel = K1D # K1D
   Im_Out = convol(Im_in, Kernel, /EDGE_WRAP)
end

;--------------------------------------------------------------------

pro smooth_median, Im_in, Im_Out, WinSize=WinSize

Im_Out = median(Im_in, WinSize)

Nx = (size(Im_in))[1]
Ny = (size(Im_in))[2]
W=WinSize/2

for i=0,Nx-1 do $
for k=0,W-1 do $
begin
   j = k
   voisin = Im_in[max([i-W,0]):min([i+W,Nx-1]), max([j-W,0]):min([j+W,Ny-1])]
   Im_Out[i,j]=voisin [ (sort(voisin))[ n_elements(voisin)/2 ] ]
   j = Ny-k-1
   voisin = Im_in[max([i-W,0]):min([i+W,Nx-1]), max([j-W,0]):min([j+W,Ny-1])]
   Im_Out[i,j]=voisin [ (sort(voisin))[ n_elements(voisin)/2 ] ]
end

for k=0,W-1 do $
for j=0,Ny-1 do $
begin
   i = k
   voisin = Im_in[max([i-W,0]):min([i+W,Nx-1]), max([j-W,0]):min([j+W,Ny-1])]
   Im_Out[i,j]=voisin [ (sort(voisin))[ n_elements(voisin)/2 ] ]
   i = Nx-k-1
   voisin = Im_in[max([i-W,0]):min([i+W,Nx-1]), max([j-W,0]):min([j+W,Ny-1])]
   Im_Out[i,j]=voisin [ (sort(voisin))[ n_elements(voisin)/2 ] ]
end
end

;--------------------------------------------------------------------

PRO IM_SMOOTH, Im_in, Im_Out, method=method, Step_Trou=Step_Trou, Border=Border, WinSize=WinSize

;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 2 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'IM_SMOOTH, Im_in, Im_Out, method=method, Step_Trou=Step_Trou, Border=Border'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
if not keyword_set(Border) then Border = 'mirror'
if not keyword_set(Step_Trou) then Step_Trou = 0
if not keyword_set(WinSize) then WinSize = 3
if not keyword_set(method) then method = 'median'

case method of
   'median': smooth_median, Im_in, Im_Out, WinSize=WinSize
   'bspline': smooth_bspline, Im_in, Im_Out, Border=Border, Step_Trou=Step_Trou
   'linear': smooth_linear, Im_in, Im_Out, Step_Trou
   end

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN
 
 END
