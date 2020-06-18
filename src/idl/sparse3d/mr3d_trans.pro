;+
; NAME:
;        MR3D_TRANS
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
;      MR3D_TRANSFORM, Imag, DataTransf, MR_File_Name=MR_File_Name, OPT=Opt
;       
;
; INPUTS:
;     Cube -- 3D IDL array: cube we want transform
;    
;     MR_File_Name: string containing the file multiresolution
;                   file name
;
; OUTPUTS:
;     DataTransf -- 3D or 4D array: multiresolution transform
;
; KEYWORDS:
;      MR_File_Name -- string: multiresolution file
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;         [-t type_of_multiresolution_transform]
;              1: (bi-) orthogonal transform  
;              2: (bi-) orthogonal transform via lifting sheme 
;              3: A trous wavelet transform 
;              default is (bi-) orthogonal transform 
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
;           [-L]
;                Use a L2 normalization. Default is L1
;            [-i]
;              print statistical information about each band.
;
;
; EXTERNAL CALLS:
;       mr3d_trans (C++ program)
;
; EXAMPLE:
;
;       Compute the multiresolution of a cube I with default options
;       (i.e. a  (bi-) orthogonal transform  with 4 scales). The result is stored in 
;       the file result.mr
;               MR3D_TRANS, I, Output, MR_File_Name='result.mr'
;
;       Compute the multiresolution of I by using the a trous algorithm 
;       with 5 scales. The result is stored in the file
;       result_at.mr
;              MR3D_TRANS, I, Output, MR_File_Name='result_at', $
;                            OPT='-t 3 -n 5'
;
; HISTORY:
;	Written: Jean-Luc Starck  
;	September, 2005 File creation
;-

pro MR3D_TRANS, Cube, DataTransf, MR_File_Name=MR_File_Name, OPT=Opt

if N_PARAMS() LT 2 then begin 
        spawn, 'mr3d_trans'
        print, 'CALLING SEQUENCE: mr3d_trans, Cube,  DataTransf, MR_File_Name=MR_File_Name, OPT=Opt'
        goto, DONE
        end

vsize = size(Cube)
if vsize(0) NE 3 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr3d_trans, Cube, DataTransf, MR_File_Name=MR_File_Name, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '
if keyword_set(MR_File_Name) then filename = MR_File_Name $
else filename = 'xx_temp.mr'
p = strpos(filename, '.mr')
if p LT 0 then filename = filename + '.mr'

NameImag='xx_imag.fits'
writefits, NameImag, Cube

com = 'mr3d_trans ' + Opt + ' ' + NameImag + ' ' +  filename
spawn, com
DataTransf = readfits(filename)

delete, NameImag
if filename EQ 'xx_temp.mr' then delete, filename
DONE:

END
