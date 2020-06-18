;+
; NAME:
;        IM_PROG
;
; PURPOSE:
;	Computes the cosinus transform of an image.
;
; CALLING:
;
;      IM_DCT, ImagIn, ImagOut, DCT_File_Name=DCT_File_Name, OPT=Opt
;       
;
; INPUTS:
;     Imag -- 2D IDL array: image we want transform
;    
;     DCT_File_Name: string containing the DCT file name
;
; OUTPUTS:
;     DataTransf -- 2D: DCT transform
;
; KEYWORDS:
;      DCT_File_Name -- string: DCT file
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;          [-r]
;             inverse cosinus transform.
;
;          [-b BlockSize]
;             Block Size. Default is image size. 
;
;          [-O]
;             Block overlapping. Default is no. 
;
;          [-w]
;             Weight the data. Default is no. 
;
;          [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4
;
;
; EXTERNAL CALLS:
;       im_dct (C++ program)
;
; EXAMPLE:
;
;       Compute the DCT of the image. The result is stored in 
;       the file result.fits
;               IM_DCT, Ima, DCT, DCT_File_Name='result.fits'
;
;       Compute the inverse DCT transform. The result is stored in 
;       the file result.fits
;              MR_TRANSFORM, DCTIma, ImaRec,  OPT='-r'
;
; HISTORY:
;	Written: Jean-Luc Starck 2003.
;	April, 2003 File creation
;-

pro IM_IRADON, ImagIn, ImagOut,  OPT=Opt, rad_header=rad_header
ImagIn = float(ImagIn)
if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: IM_RADON, ImagIn, ImagOut, OPT=Opt'
        goto, DONE
        end

vsize = size(ImagIn)
if vsize(0) NE 2 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM_RADON, ImagIn, ImagOut, OPT=Opt'
         goto, DONE
        end
 
if not keyword_set(Opt) then Opt = ' '

filename = 'xx_temp.fits'
p = strpos(filename, '.fits')
if p LT 0 then filename = filename + '.fits'

NameImag='xx_imag.rad'
writefits, NameImag, ImagIn, rad_header

PROG = "im_radon -r "
com = PROG + ' ' + Opt + ' ' + NameImag + ' ' +  filename
spawn, com
ImagOut = readfits(filename)

delete, NameImag
if filename EQ 'xx_temp.rad' then delete, filename
DONE:

END


pro ttr, d, r, z

 n = dblarr(128,128)
 n [64, 64]= 1.
 d= n 
 im_radon, d, r, rad=h, opt=' -m 5 '
 r = r * 0
 r [64, 41]= 1.
 im_iradon, r, z, rad=h
load, z
end

