PRO GET_HEAL_LUT, Table, PROJECT=Project_Choice, COORD=Coord_Choice, $
    SIZE=Size_Choice, SILENT=slnt
;+
; NAME:
;	GET_HEAL_LUT
; PURPOSE:
; 	Procedure to read and return a user-selectable look-up table.
; CALLING SEQUENCE:
;	GET_HEAL_LUT, Table, Coord = {1-9}, Project={1,2}, Size = {1,5},/SILENT
; OUTPUTS:
;	Table - The lookup table.
; OPTIONAL INPUT KEYWORDS:
;       The program will prompt for the values of COORD, PROJECT and SIZE
;           if they are not supplied as keywords.
;
;	COORD   - Scalar integer (1-9) giving the coordinate system 
;                 transformation:
;			                       Display
;			           Galactic    Celestial    Ecliptic
;			Native:
;			Galactic      1            2            3
;			Celestial     4            5            6
;			Ecliptic      7            8            9
;	PROJECT - Scalar integer (1-2) giving the type of projection:
;			1 - Mollweide
;			2 - Zenithal Equal Area
;	SIZE    - Scalar intger (1-5) giving the size of the output image 
;                 described by the table:
;			1 -- Small    (512 x 256)
;			2 -- Medium   (1024 x 512)
;			3 -- Large    (2048 x 1024)
;			4 -- X large  (4096 x 2048)
;			5 -- XX large (8192 x 4096, mollweide, native 
;			     coordinates only)
;	/SILENT  - If present and nonzero then all unnecessary output
;	          is suppressed.
; PROCEDURES USED:
;      CONCAT_DIR(), READFITS()
; MODIFICATION HISTORY:
;	Written by ?
;	SILENT keyword added.  Usage description/header fleshed out.
;	   MRG, RITSS, 16 November 2000.
;	Use getenv to extract environment variables.  MRG, SSAI, August 2002
;      Use CONCAT_DIR() to concatenate directories  WBL   May 2003
;-
on_error, 2
;
; Check arguments.
   If (n_params() LT 1) Then begin 
      print,'Syntax - GET_HEAL_LUT, Table, [Coord = {1-9}, Project={1,2},' + $
	' Size = {1,5},/SILENT] '
      return
   endif
;
mref = strtrim(getenv('MAP_REF'), 2)
Dir_String = concat_dir(mref, concat_dir('luts','heal'))

Project_Again:
IF( NOT KEYWORD_SET(Project_Choice) )THEN BEGIN
  PRINT,' Enter the type of projection you wish to employ:'
  PRINT,'    1 -- Mollweide'
  PRINT,'    2 -- Zenithal Equal Area'
  READ, Project_Choice
ENDIF
CASE Project_Choice OF
    1 : Project_String = 'mollweide'
    2 : Project_String = 'zea'
 ELSE : BEGIN
        PRINT,' Invalid value for Project_Choice, please select again...'
        Project_Choice=0
        GOTO, Project_Again
        END
ENDCASE

Coord_Again:
IF( NOT KEYWORD_SET(Coord_Choice) )THEN BEGIN
  PRINT,' Enter the number corresponding to the native cooordinates and display'
  PRINT,' coordinates of this map:'
  PRINT,'           Display:  Galactic    Celestial    Ecliptic'
  PRINT,'   Native:'
  PRINT,'   Galactic             1            2            3'
  PRINT,'   Celestial            4            5            6'
  PRINT,'   Ecliptic             7            8            9'
  READ, Coord_Choice
ENDIF
CASE Coord_Choice OF
    1 : Coord_String = '_'
    2 : Coord_String = '_g_to_c_'
    3 : Coord_String = '_g_to_e_'
    4 : Coord_String = '_c_to_g_'
    5 : Coord_String = '_'
    6 : Coord_String = '_c_to_e_'
    7 : Coord_String = '_e_to_g_'
    8 : Coord_String = '_e_to_c_'
    9 : Coord_String = '_'
 ELSE : BEGIN
        PRINT,' Invalid value for Coord_Choice, please select again...'
        Coord_Choice=0
        GOTO, Coord_Again
        END
ENDCASE

Size_Again:
IF( NOT KEYWORD_SET(Size_Choice) )THEN BEGIN
  PRINT,' Enter the size of the look-up table you wish to employ:
  PRINT,'   1 -- Small    (512 x 256)'
  PRINT,'   2 -- Medium   (1024 x 512)'
  PRINT,'   3 -- Large    (2048 x 1024)'
  PRINT,'   4 -- X large  (4096 x 2048)'
  PRINT,'   5 -- XX large (8192 x 4096, mollweide, native coordinates only)'
  READ, Size_Choice
ENDIF

if Size_Choice EQ 5 then $
    if (project_choice NE 1) or (Coord_String NE '_') then message, $
    'ERROR - XX large size available only for Mollweide native coordinates' 

CASE Size_Choice OF
    1 : Size_String = 'h10s'
    2 : Size_String = 'h10m'
    3 : Size_String = 'h10l'
    4 : Size_String = 'h10xl'
    5 : Size_String = 'h10xxl'
 ELSE : BEGIN
        PRINT,' Invalid value for Size_Choice, please select again...'
        Size_Choice=0
        GOTO, Size_Again
        END
ENDCASE

;Test whether gzip'ed versions of the LUT files are present

 LUT_File = Project_String+Coord_String+Size_String+'_lut.fits'
 fits_file = concat_dir(dir_string, lut_file) 
 openr,lun, fits_file,error=error,/get_lun
 if error NE 0 then  begin
        fits_file = fits_file + '.gz'
	openr,lun,fits_file,error=error,/get_lun
	if error NE 0 then $
	      message,'ERROR - Cannot locate look-up table ' + fits_file
	      
 endif
 free_lun,lun

; Load the look-up table
If NOT keyword_set(slnt) Then $
  message,' Reading the look-up table ' + lut_file,/INF
Table = readfits(fits_file,/SILENT)


RETURN
END
