;+
; NAME:
;        NFIRSTWT
;
; PURPOSE:
;	Computes a multiresolution transform of an image. Set to zero
;       the smallest coefficients and reconstruct the image.
;       The number of coefficient we want to preserve from the thresholding
;       can be given either by the keyword NBR or using the keyword  NFPER.
;       NFPER is the percentage of coefficients we want to keep.
;       Default transform is the orthogonal wavelet transform using five scales,
;       and a L2 normalization. The default value is 10%. 
;
; CALLING:
;
;      NFIRSTWT, Imag, Rec, Opt=Opt, nbr=nbr, nfper=nfper
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want transform
;    
;
; OUTPUTS:
;     Rec -- 2D: reconstructed image 
;
; KEYWORDS:
;      nbr -- int:  Number of coefficients we want to keep
;
;      nfper -- int: Percentage of coefficients we want to keep.
;                    Default is 10%.
;
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
;       (i.e. OWT with 5 scales, 10% of the coefficients are not thresholded). 
;               NFIRST, I, Rec
;
;       Ditto, but keep only 5%
;              NFIRST, I, Rec, nfper=5.
;
;       Ditto, but use an undecimated WT instead of the OWT
;              NFIRST, I, Rec, nfper=5., OPT='-t 24 -n 5 -L'
;
;-

pro nfirstwt, ima, rec, opt=opt, nbr=nbr, nfper=nfper

if not keyword_set(opt) then opt='-t14 -L -n5 -v'

MR_File_Name='xx_nfirst.mr'
recfilename = 'xx_temp.fits'

mr_transform, ima, wt, opt=opt, MR_File_Name=MR_File_Name
wt = readfits(MR_File_Name, header)

n = long(n_elements(wt))
if not keyword_set(nbr) then nbr = -1 $
else nfper =  float(nbr)  / float(n) * 100.

if not keyword_set(nfper) then nfper = -1
if nbr lt 0 and nfper lt 0 then nfper = 10.

nf = long(n * float(nfper) / 100.)
print, " WT, OPT = ", opt, "Total nbr of coeff = ", n
print, " Keep percentage P = ", nfper, " ==> KEEP NFIRST COEF: NF = ", nf

t = fltarr(n)
t(*) = wt
t1 = abs(t)
a = sort(t1)

t(a(0:n-1-nf)) = 0
wt(*)= t(*)
writefits, MR_File_Name, wt, header

com = 'mr_recons ' + MR_File_Name + ' ' +  recfilename
spawn, com
rec = readfits(recfilename)
delete, recfilename
delete, MR_File_Name
load, rec
end
