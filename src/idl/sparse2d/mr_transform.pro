;+
; NAME:
;        MR_TRANSFORM
;
; PURPOSE:
;	Computes a multiresolution transform of an image. If the
;       the keyword MR_File_Name is set, a file is created
;       which contains the multiresolution transform. A multiresolution file
;       has a .mr extension, if the parameter file name have not one, then
;       the extension is added to the file name. Result is store in DataTransf.
;       Dependiong on the options (see OPT), DataTransf can be a cube or an
;       image. If the keyword is not set, the created file is deleted.
;
; CALLING:
;
;      MR_TRANSFORM, Imag, DataTransf, MR_File_Name=MR_File_Name, OPT=Opt
;       
;
; INPUTS:
;     Imag -- 2D IDL array: image we want transform
;    
;     MR_File_Name: string containing the file multiresolution
;                   file name
;
; OUTPUTS:
;     DataTransf -- 2D or 3D array: multiresolution transform
;
; KEYWORDS:
;      MR_File_Name -- string: multiresolution file
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
;
; EXAMPLE:
;
;       Compute the multiresolution of the image I with default options
;       (i.e. a trou algorithm with 4 scales). The result is stored in 
;       the file result.mr
;               MR_TRANSFORM, I, Output, MR_File_Name='result.mr'
;
;       Compute the multiresolution of I by using the pyramidal median 
;       algorithm with 5 scales. The result is stored in the file
;       result_pyr_med.mr
;              MR_TRANSFORM, I, Output, MR_File_Name='result_pyr_med', $
;                            OPT='-t 10 -n 5'
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;       October, 1999 Header Update
;-

pro MR_TRANSFORM, Imag, DataTransf, MR_File_Name=MR_File_Name, OPT=Opt

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_transform'
        print, 'CALLING SEQUENCE: mr_transform, Imag,  DataTransf, MR_File_Name=MR_File_Name, OPT=Opt'
        goto, DONE
        end

vsize = size(Imag)
if vsize(0) NE 2 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr_transform, Imag,  DataTransf, MR_File_Name=MR_File_Name, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '
if keyword_set(MR_File_Name) then filename = MR_File_Name $
else filename = 'xx_temp.mr'
p = strpos(filename, '.mr')
if p LT 0 then filename = filename + '.mr'

NameImag='xx_imag.fits'
writefits, NameImag, Imag

com = 'mr_transform ' + Opt + ' ' + NameImag + ' ' +  filename
spawn, com
DataTransf = mr_read(filename)

delete, NameImag
if filename EQ 'xx_temp.mr' then delete, filename
DONE:

END
