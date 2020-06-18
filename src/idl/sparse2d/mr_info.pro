;+
; NAME: 
;       MR_INFO
;
; PURPOSE: 
;       Calculate statistical information about the wavelet transform
;       of an image.
;
; CALLING:
;       MR_Info, Imag, TabStat, OPT=Opt, nodel=nodel, NameRes=NameRes
;
; INPUT:
;       Imag: image 
;
; KEYWORDS:
;
;      NameRes: Outfile file name for the statistic table
;
;      nodel: if set, the created files (by the program mr_info)
;             are not deleted.
;
;      Opt: string which contains the differents options. Options are:
;
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
;                 25: Wavelet transform via lifting scheme 
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
;         [-n number_of_scales]
;             number of scales used in the multiresolution transform
;              default is 5
;
;         [-L]
;             Use a L2 normalization. Default is L1
;
;         [-u number_of_undecimated_scales]
;                Number of undecimated scales used in the 
;                Undecimated Wavelet Transform
;
;         [-g sigma]
;             sigma = noise standard deviation
;              default is automatically estimated
;
;         [-c gain,sigma,mean]
;             gain = gain of the CCD
;             sigma = read-out noise standard deviation
;             mean = read-out noise mean
;               noise = poisson + readout noise. default is no (Gaussian)
;
;         [-p]
;            Poisson Noise
;             default is no (Gaussian).
;
;         [-s nsigma]
;              Thresholding at nsigma * SigmaNoise
;              default is 3.000000
;
;         [-k]
;           Suppress isolated pixels in the support. Default is no.
; 
;         [-a]
;              Significant structures analysis.
;              default is no.
;
;         [-l] 
;              Dilate the support
;	     
;         [-v] 
;              Verbose
;
; OUTPUTS:
;           TabStat: fltarr(NbrBand, NbrStat) 
;                    TabStat(j,0) = standard deviation of the jth band
;                    TabStat(j,1) = skewness of the jth band
;                    TabStat(j,2) = kurtosis of the jth band
;                    TabStat(j,3) = minimum of the jth band
;                    TabStat(j,4) = maximum of the jth band
;
; EXTERNAL CALLS:
;           mr_info (C++ program)
;
; EXAMPLE:
;           calculate the statistics with all default options 
;                mr_info, Imag, TabStat
;
; HISTORY:
;	Written: Jean-Luc Starck 2002.
;-

PRO mr_info, imag, result, OPT=OPT, nodel=nodel, NameRes=NameRes

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_info'
        print, 'CALL SEQUENCE: mr_info, Imag, Imag_Out, OPT=Opt, nodel=nodel'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
if keyword_set(NameRes) then NameResult = NameRes $
else NameResult = 'xx_result.fits' 

writefits, NameImag, imag 
com = 'mr_info ' + OPT + ' ' + NameImag + ' ' + NameResult
; print, com
spawn, com

Result = readfits(NameResult, /silent)
;help, result
if not keyword_set(nodel) then begin
delete, NameImag
if not keyword_set(NameRes) then delete, NameResult
end

DONE:
end
