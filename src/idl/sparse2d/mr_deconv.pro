;+
; NAME: 
;       MR_DECONV
;
; PURPOSE: 
;       Deconvolve of an image by using a multiresolution transform.
;       The deconvolution method is described in:
;       J.L. Starck, A. Bijaoui, and F. Murtagh, "Multiresolution Support
;       Applied to Image Filtering and Deconvolution", in CVIP: Graphical
;       Models and Image Processing, Vol. 57, 5, pp 420-431, Sept. 1995.
;
;
; CALLING:
;
;       MR_DECONV, Imag, Psf, Result, Opt=Opt
;
; INPUTS:
;       Imag: image to deconvolve
;       Psf:  Point Spread Function
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;           [-d type_of_deconvolution]
;               type_of_deconvolution = 
;                  1: Deconvolution by multiresolution Van Cittert's algorithm
;                  2: Deconvolution by multiresolution Gradient's algorithm
;                  3: Deconvolution by multiresolution Lucy's algorithm 
;                  4: Deconvolution by multiresolution MAP algorithm 
;                  5: Deconvolution by the division in Fourier space 
;                                                          + Wavelet filtering 
;                  default is 3
;
;          [-t type_of_multiresolution_transform]
;                  1: linear wavelet transform: a trous algorithm 
;                  2: bspline wavelet transform: a trous algorithm 
;                  3: wavelet transform in Fourier space 
;                  4: morphological median transform 
;                  5: morphological minmax transform 
;                  6: pyramidal linear wavelet transform 
;                  7: pyramidal bspline wavelet transform 
;                  8: pyramidal wavelet transform in Fourier space: 
;                     wavelet =  between two resolutions 
;                  9: pyramidal wavelet transform in Fourier space: 
;                     wavelet = difference  between the square 
;                                                of two resolutions
;                 10: pyramidal median transform 
;                 11: morphological pyramidal minmax transform 
;                 12: pyramidal laplacian 
;                 13: decomposition on scaling function 
;                 14: Mallat's wavelet transform (7/9 filters) 
;                 15: Feauveau's wavelet transform 
;                 16: Feauveau's wavelet transform without undersampling 
;                 17: G transform (morphological min-max algorithm) 
;                 18: Haar's wavelet transform 
;                 19: Half-pyramidal transform 
;                 20: Mixed Half-pyramidal WT and Median method (WT-HPMT) 
;                 21: diadic wavelet transform 
;                 22: Mixed WT and PMT method (WT-PMT) 
;                 23: undecimated Haar transform: a trous algorithm 
;                 24: undecimated mallat transform (three bands per scale)
;                 default is bspline wavelet transform: 
;                                       a trous algorithm
;
;          [-T type_of_filters]
;                 1: Antonini 7/9 filters 
;                 2: Daubechies filter 4 
;                 3: Biorthogonal 2/6 Haar filters 
;                 4: Biorthogonal 2/10 Haar filters 
;                 5: Odegard 7/9 filters 
;                 6: User's filters 
;                 default is Antonini 7/9 filters
;                
;           [-u number_of_undecimated_scales]
;                Number of undecimated scales used in the 
;                Undecimated Wavelet Transform
;
;           [-g sigma]
;                Gaussian noise
;                  sigma = noise standard deviation 
;                by default, the noise is gaussian, and the standard 
;                deviation is automatically estimated. 
;
;           [-c gain,sigma,mean]
;                case of a CCD: noise = Poisson noise + read-out noise
;                  gain = CCD gain 
;                  sigma = standard deviation of the read-out noise
;                  mean = mean of the read-out noise
;                if this option is set, 
;                           Noise = Poisson + Gaussian read-out Noise
;                it is generally the case with the CCD.
;                Attention, these parameters must be separated by a comma 
;                without space. example: -c 0.133,1.733,0.
;                If mean, or sigma and mean are omitted, default values are 0.
;                gain can not be omitted. 
;
;         [-m type_of_noise]
;              1: Gaussian noise 
;              2: Poisson noise 
;              3: Poisson noise + Gaussian noise 
;              4: Multiplicative noise 
;              5: Non-stationary additive noise 
;              6: Non-stationary multiplicative noise 
;              7: Undefined uniform noise 
;              8: Undefined noise 
;              9: Stationary correlated noise 
;             default is Gaussian noise
;
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4
;
;           [-s NSigma]
;                Thresolding at NSigma * SigmaNoise at each scale
;                default is 3
;
;           [-i number_of_iterations]
;              Maximum number of iterations
;              default is 500
;
;           [-e Epsilon]
;                Convergence parameter
;                default is 0.001
;   
;           [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options). 
;
;           [-f ICF_Fwhm]
;                Intrinsic correlation function.
;                Fwhm = Full Width at Half Maximum.
;  
;           [-P]
;               Suppress the positivity constraint.
;
;           [-I ICF_FileName]
;                Intrinsic correlation function file.
;
;           [-F First_Guess]
;                Input solution file.
;
;           [-r residual_file_name]
;               Residual_file_name = file name
;               write the residual to the disk 
; 
;           [-S]
;               Do not shift automatically the maximum  
;               of the PSF at the center.
;
;           [-p]
;               Detect only positive structure. Default is no.
;   
;           [-K]
;               Suppress the last scale. Default is no.
;   
;           [-O]
;               Optimization.
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           mr_deconv (C++ program)
;
; EXAMPLE:
;           deconvolve an image with all default options (Lucy method + 
;           regularization in the wavelet space by using the a-trou 
;           algorithm, ...)
;                MR_DECONV, Imag, Psf, Result
;
;          same example, but impose the number of iterations to be 30
;                MR_DECONV, Imag, Psf, Result, OPT='-i 30 -e 0'
;
;          deconvolution by the one step gradient method, without 
;          any regularization, with 30 iterations
;                MR_DECONV, Imag, Psf, Result, OPT='-d 2 -i 30'
;
; HISTORY:
;	Written: Jean-Luc Starck 1996.
;	December, 1996 File creation
;       October, 1999 Header Update
;-

PRO mr_deconv, imag, Psf, result, OPT=OPT

if N_PARAMS() LT 3 then begin 
        spawn, 'mr_deconv'
        print, 'CALL SEQUENCE: mr_deconv Imag Psf OPT=Opt'
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

com = 'mr_deconv ' + OPT + ' ' + NameImag + ' ' + NamePsf + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
delete, NamePsf
DONE:
end
