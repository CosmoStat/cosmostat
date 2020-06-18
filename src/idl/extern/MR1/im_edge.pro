;PRO IM_EDGE, Im_in, Im_Out, opt=opt
;+ 
; NAME: 
;     IM_SMOOTH
;
; PURPOSE: 
;     Caculate the edges of an image
;
; CATEGORY: 
;     III-3
;
; CALLING SEQUENCE: 
;   IM_EDGE, Im_in, Im_Out, opt=opt
;
; INPUTS: 
;   Im_in -- 2D IDL array: input image  
;
; KEYED INPUTS: 
;      Opt: string which contains the differents options. Options are:
; 
;         [-M edge_detection_method]
;              1: Canny (derivative of a Gaussian) 
;              2: Pixel difference 
;              3: Separated pixel difference 
;              4: Sobel 
;              5: Prewitt 
;              6: Roberts 
;              7: Frei Chen 
;              8: Laplacian: filter 1 
;              9: Laplacian: filter 2 
;              10: Laplacian: filter 3 
;              11: Marr and Hildrith: Laplacian of Gaussian (LoG) 
;              12: Prewitt compass gradient 
;              13: Kirsch 
;              14: Robinson 3-level 
;              15: Robinson 5-level 
;             default is Sobel method
;         [-k]
;             Select maxima or zero-crossing.
;         [-t ThresholdValue]
;             Threshold the edge map.
;         [-S Sigma]
;             Scale parameter. Used with DroG and LoG methods.
;             Default is 1/sqrt(3) 
;
; OUTPUTS: 
;   Im_Out -- 2D IDL array: edge image 
; 
; EXTERNAL CALLS:
;           im_edge (C++ program)
;
;  EXAMPLE:
;           edge detection with all default options
;                IM_EDGE, Imag, Edge
;
;          same example, but using Canny menthod and edge thinning.
;                IM_EDGE, Imag, Edge, OPT='-M 1 -k'
;
;
; MODIFICATION HISTORY: 
;    20-oct-2002 JL Starck  
;-
 
;==============================
 
  
PRO IM_EDGE, imag, Result, OPT=OPT

;------------------------------------------------------------
; parameters check
;------------------------------------------------------------

if N_PARAMS() LT 2 then begin 
        spawn, 'im_edge'
        print, 'CALL SEQUENCE: im_edge Imag Edge OPT=Opt'
        goto, DONE
        end


 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

if not keyword_set(Opt) then Opt = ' '  
Nl = (size(imag))[2]
Nc = (size(imag))[1]
NameImag = 'xx_imag.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag
 
com = 'im_edge ' + OPT + ' ' + NameImag + ' ' + NameResult
print, com
spawn, com
Result = readfits(NameResult, /silent)

; delete, NameImag
; delete, NameResult

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 DONE:
 
  RETURN
 
 END
