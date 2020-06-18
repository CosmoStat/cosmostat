Pro HEALPIX_Nested_Vectors, Resolution, Px, Py, Pz
;+
;  NAME
;     HEALPIX_Nested_Vectors
; PURPOSE:
;     Read Cartesian [X,Y,Z] direction vectors for HEALPIX pixel scheme 
; EXPLANATION:
;     IDL procedure to read Cartesian [X,Y,Z] direction vectors from (FITS) 
;     binary file for HEALPIX pixel scheme, nested pixel order.
; CALLING SEQUENCE:
;     HEALPIX_Nested_Vectors, resolution, Px, Py, Pz
; INPUTS:
;     Resolution - Scalar integer (3-9) giving the HEALPIX resolution
; OUTPUTS:
;     Px, Py, Pz, - vectors giving the  Cartesan direction vectors for the 
;         specified resolution.    The number of elements will be 768, 3072,
;         12288, 49152, 196608, 786432 or 3145728 as the resolution is 
;         increased from 3 to 9.         
; PROCEDURES USED:
;     CONCAT_DIR(), READFITS()
; REVISION HISTORY:
;      Written,  Al Kogut     May 26, 1998
;      (modified for FITS format 1 May 2000 JW)
;       se getenv to extract environment variables.  MRG, SSAI, 20 August 2002.
;      Use CONCAT_DIR() to concatenate directories  WL   May 2003
;-
;--------------------------------------------------------------------------

 if N_Params() LT 2 then begin
      print,'Syntax - HEALPIX_Nested_Vectors, resolution, Px, Py, Pz'
      return
 endif
 
; Find appropriate file in $MAP_REF/pixels
 In_Dir = concat_dir( strtrim(getenv('MAP_REF'), 2), 'pixels')
 Res_String = Strtrim( Resolution, 2)
 In_File = concat_dir( In_Dir,'healpix_vectors_nested_' + Res_String + '.fits')

 vector_array = readfits(In_File,/silent)
 Px = vector_array[*,0]
 Py = vector_array[*,1]
 Pz = vector_array[*,2]

Return
End
