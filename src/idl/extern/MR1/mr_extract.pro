;+
; NAME:
;        MR_EXTRACT
;
; PURPOSE:
;	Extract a scale from a wavelet transform
;
; CALLING:
;
;      MR_Extract, Multiresolution_File_Name, ScaleImage, OPT=Opt
;       
;
; INPUTS:
;     Multiresolution_File_Name: string containing the file multiresolution
;                                file name
;
; OUTPUTS:
;     ScaleImage: Image of a scale of the multiresolution transform
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;
;        where options = 
;
;           [-b BandNumber] 
;                Extract a band.

;           [-s scale_number]
;                scale number to extract. default is 1.
;
;           [-x]
;                Write all bands separately as images with prefix 'band_j' 
;               (j being the scale number)
; 
;           [-B] 
;               Same as x option, but interpolate by block the scales.
;               This option is valid only if the chosen multiresolution 
;               transform is pyramidal (6,7,8,9,10,11,12). 
;
; EXTERNAL CALLS:
;       mr_extract (C++ program)
;
; EXAMPLE:
;
;       Extracts the first scale of the multiresolution file filetrans 
;               MR_EXTRACT, 'filetrans.mr', Imag_band1
;
;       Extract the scale 3 of  of the multiresolution file filetrans
;               MR_TRANSFORM, 'filetrans.mr', Imag_band3, OPT='-b 3'
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
pro MR_EXTRACT, MR_File_Name, Imag, OPT=Opt 

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_extract'
        print, 'CALL SEQUENCE: mr_extract MRFileName Image_Out OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '

OutName='xx_scale.fits'
com = 'mr_extract ' + Opt  + ' ' +  MR_File_Name + ' ' + OutName

spawn, com

Imag = rim(OutName)

delete, OutName
DONE:

END
