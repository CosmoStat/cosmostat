;+
; NAME: 
;       MR_FILTER
;
; PURPOSE: 
;       filter  an image by using a multiresolution transform.
;       
; CALLING:
;
;       MR_FILTER, Imag, Result, Opt
;
; INPUTS:
;       Imag: image to filter
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;         [-f type_of_filtering]
;              1: Multiresolution Hard K-Sigma Thresholding  
;              2: Multiresolution Soft K-Sigma Thresholding
;              3: Iterative multiresolution thresholding 
;              4: Adjoint operator applied to the multiresolution support  
;              5: Hierarchical thresholding 
;              6: Hierarchical Wiener filtering 
;              7: Multiresolution Wiener filtering 
;              8: Median filtering 
;              9: Average filtering 
;              10: Bspline filtering 
;              11: Donoho hard thesholding 
;              12: Donoho soft thesholding 
;              13: SURE Hard Thesholding 
;              14: SURE Soft Thesholding 
;              15: MULTI-SURE Hard Thesholding 
;              16: MULTI-SURE Soft Thesholding 
;              17: Median Absolute Deviation (MAD) Hard Thesholding 
;              18: Median Absolute Deviation (MAD) Soft Thesholding 
;              default is Multiresolution Hard K-Sigma Thresholding.
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
;                 11: pyramidal laplacian 
;                 12: morphological pyramidal minmax transform 
;                 13: decomposition on scaling function 
;                 14: Mallat's wavelet transform 
;                 15: Feauveau's wavelet transform 
;                 16: Feauveau's wavelet transform without undersampling 
;                 17: G transform (morphological min-max algorithm 
;                 18: Haar's wavelet transform 
;                 19: Half-pyramidal transform 
;                 20: Mixed Half-pyramidal WT and Median method (WT-HPMT) 
;                 21: diadic wavelet transform 
;                 22: Mixed WT and PMT method (WT-PMT) 
;                 23: undecimated Haar transform: a trous algorithm 
;                 24: undecimated mallat transform (three bands per scale)
;                 default is 2
;
;          [-T type_of_filters]
;              1: Antonini 7/9 filters 
;              2: Daubechies filter 4 
;              3: Biorthogonal 2/6 Haar filters 
;              4: Biorthogonal 2/10 Haar filters 
;              5: Odegard 7/9 filters 
;              6: User's filters 
;              default is Antonini 7/9 filters
;
;           [-u number_of_undecimated_scales]
;                Number of undecimated scales used in the 
;                Undecimated Wavelet Transform
;
;           [-g sigma]
;                Gaussian noise
;                  sigma = noise standard deviation 
;                by default, the noise is gaussian, and the standard 
;                devaition is automatically estimated. 
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
;           [-m]
;              1: Gaussian Noise 
;              2: Poisson Noise 
;              3: Poisson Noise + Gaussian Noise 
;              4: Multiplicative Noise 
;              5: Non uniform additive noise 
;              6: Non uniform multiplicative noise 
;              7: Undefined uniform Noise 
;              8: Undefined Noise 
;              9: Stationary correlated noise 
;              10: Poisson Noise with few events 
;             default is Gaussian noise
;
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4.
;                default is 6 in case of poisson noise with few events 
;
;           [-s NSigma]
;                Thresolding at NSigma * SigmaNoise at each scale
;                default is 3
;
;           [-i number_of_iterations]
;              Maximum number of iterations
;              default is 10
;              default is 0.000010 in case of poisson noise with few events
;
;           [-e Epsilon]
;                Convergence parameter
;                default is 0.0001
;
;           [-w support_file_name]
;              support_file_name = file name
;              creates an image from the support 
;              and write it on the disk
;
;           [-k]
;              Suppress isolated pixel in the support.
;              Default is no.
;
;           [-K]
;              Suppress the last scale. Default is no.
;
;           [-p]
;              Detect only positive structure 
;              default is no
;
;           [-E Epsilon]
;             Epsilon = precision for computing thresholds
;                       (only used in case of poisson noise with few events)
;             default is 1.000000e-03 
;
;
;           [-S SizeBlock]
;             Size of the  blocks used for local variance estimation.
;             default is 7 
;
;
;           [-N NiterSigmaClip]
;             iteration number used for local variance estimation.
;             default is 1 
;
;
;           [-F first_detection_scale]
;              first scale used for the detection 
;              default is 1
;
;
;           [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options).
;
;           [-W WindowSize]
;               Window size for median and average filtering.
;               default is 5.
;
;           [-b]
;               Add the maximum level constraint.
;               Max value is 255. Default is no.
;
;           [-v]
;              Verbose. Default is no.
;
;           [-P]
;              Suppress the positivity constraint.
;
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           mr_filter (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                MR_FILTER, Imag, Result
;
;          same example, but impose the number of scales  to be 3, and 
;          a thresolding at 5 sigma
;                MR_FILTER, Imag, Psf, Result, OPT='-n 3 -s 5'
;
;
; HISTORY:
;	Written: Jean-Luc Starck 1994.
;	December, 1995 File creation
;       October, 1999 Header Update
;-

PRO mr_filter, imag, result, OPT=OPT

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_filter'
        print, 'CALL SEQUENCE: mr_filter Imag Imag_Out OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag

com = 'mr_filter ' + OPT + ' ' + NameImag + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
DONE:
end
