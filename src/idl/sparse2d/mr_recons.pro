;+
; NAME:
;        MR_RECONS
;
; PURPOSE:
;	Reconstruct an image from its multiresolution transform. If the
;       the keyword MR_File_Name is set, a multiresolution transform   
;       must already exist. A multiresolution file
;       has a .mr extension, if the parameter file name have not one, then
;       the extension is added to the file name. If the MR_File_Name is not
;       set, then the Fits Header of the multiresolution file must be given.
;       Result is store in Imag.
;
; CALLING:
;
;      MR_RECONS, DataTransf, Imag, MR_File_Name=MR_File_Name, Header=Header
;       
;
; INPUTS:
;     DataTransf -- 2D or 3D array: multiresolution transform
;    
; OUTPUTS:
;     Imag -- 2D IDL array: image we want transform;
;
; KEYWORDS:
;      MR_File_Name -- string: multiresolution file
;
;      Header -- string: FITS Header
;
;
; EXTERNAL CALLS:
;       mr_recons (C++ program)
;
; EXAMPLE:
;
;       Compute the multiresolution of the image I with default options
;       (i.e. a trou algorithm with 4 scales). The result is stored in 
;       the file result.mr
;             MR_TRANSFORM, I, Output, MR_File_Name='result.mr'
;
;       Reconstruct the image from the transformation
;              MR_RECONS, Output, RecIma, MR_File_Name='result_pyr_med'
;
; HISTORY:
;	Written: Jean-Luc Starck 2003.
;	April, 2003 File creation
;-

pro MR_Recons, DataTransf, Imag, MR_File_Name=MR_File_Name, Header=Header

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_recons'
        print, 'CALLING SEQUENCE: mr_recons, DataTransf, Imag, MR_File_Name=MR_File_Name, Header=Header'
        goto, DONE
        end

if not keyword_set(MR_File_Name) and  not keyword_set(header) then begin
        print, 'CALLING SEQUENCE: mr_recons, DataTransf, Imag, MR_File_Name=MR_File_Name, Header=Header'
        print, '                 MR_File_Name or Header keyword must be set'
	goto, DONE
        end
	
	
vsize = size(DataTransf)
if vsize(0) LT 1 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr_recons, DataTransf, Imag, MR_File_Name=MR_File_Name, Header=Header'
        goto, DONE
        end

if keyword_set(MR_File_Name) then begin
   filename = MR_File_Name
   p = strpos(filename, '.mr')
   if p LT 0 then filename = filename + '.mr'
   w = readfits(filename, header)
   end else filename = 'xx_mr.mr'


NameImag='xx_imag.fits'

writefits, filename , DataTransf, header

com = 'mr_recons '  + filename + ' ' +  NameImag
spawn, com
Imag = readfits(NameImag)

delete, NameImag
; delete, filename
DONE:

END
