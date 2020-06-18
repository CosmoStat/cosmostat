;+
; NAME:
;        MR3D_RECONS
;
; PURPOSE:
;     	Reconstruct an image from its multiresolution transform. If the
;       the keyword MR_File_Name is set, a multiresolution transform   
;       must already exist. A multiresolution file
;       has a .mr extension, if the parameter file name have not one, then
;       the extension is added to the file name. If the MR_File_Name is not
;       set, then the Fits Header of the multiresolution file must be given.
;       Result is store in Cube.
;
;
; CALLING:
;
;      MR3D_RECONS,  DataTransf, Cube, MR_File_Name=MR_File_Name, Header=Header
;
;
; INPUTS:
;     DataTransf -- 3D or 4D array: multiresolution transform
;     
;     MR_File_Name: string containing the file multiresolution
;                   file name
;
; OUTPUTS:
; 	 Cube -- 3D IDL array: cube we want transform
;
; KEYWORDS:
;      MR_File_Name -- string: multiresolution file
;      Header -- string: FITS Header
;      
;
; EXTERNAL CALLS:
;       mr3d_recons (C++ program)
;
; EXAMPLE:
;       Compute the multiresolution of the Cube with default options
;       (i.e. a trou algorithm with 4 scales). The result is stored in 
;       the file result.mr
;             MR3D_TRANS, C, Output, OPT='-t3 -n4', header
;
;       Reconstruct the image from the transformation
;              MR3D_RECONS, Output, RecCube, header = header
; HISTORY:
;	Written: Sandrine Pires  
;	September, 2006 File creation
;-

pro MR3D_RECONS, DataTransf, Cube, MR_File_Name=MR_File_Name, Header=Header

if N_PARAMS() LT 2 then begin 
        spawn, 'mr3d_recons'
        print, 'CALLING SEQUENCE: mr3d_recons, DataTransf, Cube, MR_File_Name=MR_File_Name, Header=Header'
        goto, DONE
        end
if not keyword_set(MR_File_Name) and  not keyword_set(header) then begin
        print, 'CALLING SEQUENCE: mr3d_recons, DataTransf, Cube, MR_File_Name=MR_File_Name, Header=Header'
        print, 'MR_File_Name or Header keyword must be set'
	goto, DONE
        end
				
vsize = size(DataTransf)
if vsize(0) NE 4 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr3d_recons, DataTransf, Cube, MR_File_Name=MR_File_Name, Header=Header'
        goto, DONE
        end

if keyword_set(MR_File_Name) then begin
   filename = MR_File_Name
   p = strpos(filename, '.mr')
   if p LT 0 then filename = filename + '.mr'
   w = readfits(filename, header)
   end else filename = 'xx_mr.mr'


NameCube='xx_imag.fits'

writefits, filename , DataTransf, header

com = 'mr3d_recons '  + filename + ' ' +  NameCube
spawn, com
Cube = readfits(NameCube)

delete, NameCube
; delete, filename
DONE:

END
