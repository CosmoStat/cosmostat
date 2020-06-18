;+
; NAME:
;        NFIRSTINMAX
;
; PURPOSE:
;	Computes a multiresolution transform of an image. Set to zero
;       the smallest coefficients and reconstruct the image.
;       The number of coefficient we want to preserve from the thresholding
;       is given by the number  NFPER.
;       All coeff with an absolute value larger than MAX-NFPER*MAX are kept
;       Default transform is the orthogonal wavelet transform using five scales,
;       and a L2 normalization. The default value is 10%. 
;
; CALLING:
;
;      NFIRSTINMAX, Imag, Rec, Opt=Opt, nfper 
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want transform
;
;     nfper -- int [0,1]: Percentage of coefficients we want to keep.
;
; OUTPUTS:
;     Rec -- 2D: reconstructed image 
;
; KEYWORDS:
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
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
;          [-l type_of_lifting_transform]
;              1: Lifting scheme: CDF WT 
;              2: Lifting scheme: median prediction 
;              3: Lifting scheme: integer Haar WT 
;              4: Lifting scheme: integer CDF WT 
;              5: Lifting scheme: integer (4,2) interpolating transform 
;              6: Lifting scheme: 7/9 WT 
;              7: Lifting scheme: integer 7/9 WT 
;             default is Lifting scheme: integer Haar WT
;   
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4
;
;           [-x] 
;                write all scales separatly as images with prefix "scale_j"
;                (j being the scale number)
;
;           [-B] 
;                same as x option, but interpolate by block the scales.
;                This option is valid only if the choosen multiresolution 
;                transform is pyramidal (6,7,8,9,10,11,12)
;
;           [-c iter] 
;                iterative transformation. Iter = number of iterations.
;                This option is valid only if the choosen multiresolution 
;                transform is pyramidal (6,7,8,9,10,11). The reconstruction
;                is not exact and we need few iterations. Generally 3
;                iterations are enough.
;
;           [-u number_of_undecimated_scales]
;                Number of undecimated scales used in the 
;                Undecimated Wavelet Transform
;
;           [-L]
;                Use a L2 normalization. Default is L1
;
; EXTERNAL CALLS:
;       mr_transform (C++ program)
;       mr_recons (C++ program)
;
; EXAMPLE:
;
;       Compute the multiresolution of the image I with default options
;       (i.e. OWT with 5 scales,  coefficients  within 10% of the max are not thresholded. 
;               NFIRST, I, Rec, 0.1
;
;-

pro lc_nfirstinmax, ima, rec, nfper, opt=opt, nsx=nsx, nsy=nsy


if not keyword_set(opt) then opt=' '

vs = size(ima)
Nx = vs[1]
Ny = vs[2]
if keyword_set(nsx) then ScaleX  = nsx else ScaleX  = fix(  alog(float(Nx) / 4. * 3.) / alog(2.))  
if keyword_set(nsy) then ScaleY  = nsy else ScaleY  = fix(  alog(float(Ny) / 4. * 3.) / alog(2.))  
 
 
DataFileName='xx_data.fits'
MR_File_Name='xx_nfirst.fits'
recfilename = 'xx_temp.fits'

writefits, DataFileName, ima
optx = ' -x ' + STRCOMPRESS(STRING(ScaleX), /REMOVE_ALL)
opty = ' -y ' + STRCOMPRESS(STRING(ScaleY), /REMOVE_ALL)
lastNx = Nx / 2^(ScaleX-3)
lastNy = Ny / 2^(ScaleY-3)
; print, lastNx, lastNy

com = 'mr_linecol ' + opt + optx + opty + ' ' + DataFileName + ' ' +  MR_File_Name
; print, com
spawn, com
wt = readfits(MR_File_Name, header)

n = long(n_elements(wt))

y = abs(wt)
y1 = y
y1(0:lastNx-1,0:lastNy-1) = 0
m = max(y1)
t = m - nfper*m
index = where(y lt t,c)
print, " LC => OPT = ",   opt + optx + opty, "   Total nbr of coeff = ", n
print, " Keep perc. P in max = ", nfper, " MAX = ", m, " Level = ", t, " ==> KEEP NFIRST COEF: NF = ", n - c
;   print, "MAX = ", m, " nfper in max = ", nfper, "  ==> Perc. of threshold  coeff = ", float(c) / n * 100. 
if c gt 0 then begin
    wt(index) = 0
   end 
  
writefits, MR_File_Name, wt, header

com = 'mr_linecol -r ' + opt + optx + opty + ' ' + MR_File_Name + ' ' +  recfilename
; print, com
spawn, com
rec = readfits(recfilename)
delete, DataFileName
delete, recfilename
delete, MR_File_Name
; load, rec
end

