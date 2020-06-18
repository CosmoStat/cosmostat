;+
; NAME: 
;       IM_DECONV
;
; PURPOSE: 
;       Deconvolve of an image using standard methods
;
;
; CALLING:
;
;       IM_DECONV, Imag, Psf, Result, Opt=Opt
;
; INPUTS:
;       Imag: image to deconvolve
;       Psf:  Point Spread Function
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;           [-d type_of_deconvolution]
;
;               type_of_deconvolution = 
;                 1: Deconvolution by Van Cittert's algorithm 
;                 2: Deconvolution by gradient algorithm 
;                 3: Deconvolution by division in Fourier space 
;                 4: Deconvolution by Lucy's algorithm 
;                 5: Deconvolution by CLEAN algorithm 
;                 6: Deconvolution by the MEM method (Frieden entropy) 
;                 7: Deconvolution by the MEM method (Gull entropy) 
;                 8: Deconvolution using Tikhonov Regularization 
;                 9: Deconvolution using MAP method
; 
;                 default is Deconvolution by gradient algorithm 
;
;           [-i number_of_iterations]
;              Maximum number of iterations
;              default is 500
;
;            [-P]
;             Suppress the positivity constraint.
;
;           [-G RegulParam]
;                Regularization parameter 
;              default is 0.000000
;  
;           [-C ConvergParam]
;                Convergence parameter. 
;                default is 1.000000
;  
;           [-f ICF_Fwhm]
;                Intrinsic correlation function.
;                Fwhm = Full Width at Half Maximum.
;  
;           [-I ICF_FileName]
;                Intrinsic correlation function file.
;  
;           [-F First_Guess]
;                Input solution file.
;  
;           [-O]
;                Optimization.
;  
;           [-M Model_Image]
;                Input model for MEM method.
;  
;           [-r residual_file_name]
;               Residual_file_name = file name
;               write the residual to the disk 
;  
;           [-S]
;               Do not shift automatically the maximum  
;               of the PSF at the center.
;  
 
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           im_deconv (C++ program)
;
; EXAMPLE:
;           deconvolve an image with all default options (Gradient method)
;                IM_DECONV, Imag, Psf, Result
;
;          same example, but impose the number of iterations to be 30
;                IM_DECONV, Imag, Psf, Result, OPT='-i 30 -e 0'
;
;          deconvolution by the Lucy method, without 
;          any regularization, with 30 iterations
;                IM_DECONV, Imag, Psf, Result, OPT='-d 4 -i 30'
;
; HISTORY:
;	Written: Jean-Luc Starck 1996.
;	December, 1996 File creation
;-

PRO im_deconv, imag, Psf, result, OPT=OPT

if N_PARAMS() LT 3 then begin 
        spawn, 'im_deconv'
        print, 'CALL SEQUENCE: im_deconv Imag Psf OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
NamePsf = 'xx_psf.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag
writefits, NamePsf, Psf

com = 'im_deconv ' + OPT + ' ' + NameImag + ' ' + NamePsf + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
delete, NamePsf
DONE:
end
