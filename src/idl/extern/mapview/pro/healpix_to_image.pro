PRO HEALPIX_TO_IMAGE, T, Image, MASK=N_Obs, SCALE=Scale_Choice, $
    TRANSFER=Transfer_Choice, MIN=Tmin, MAX=Tmax, SILENT=slnt, _EXTRA=Extra 
;+
; NAME:
;       HEALPIX_TO_IMAGE
; PURPOSE:
;       To convert a 1-D map in 'healpix' pixel list format into a 2-D
;       byte-scaled image ready for display.
; CALLING SEQUENCE:
;       HEALPIX_TO_IMAGE, T, Image
; INPUTS:
;       T     - The map.
; OUTPUTS:
;       Image - The 2-D byte image.
; OPTIONAL INPUT KEYWORDS:
;       MASK     - A mask array indicating whether or not the unobserved
;                  pixels are painted with color 2 (see COMMENTS).
;       SCALE    - The display scaling option:
;                       1 - Auto-scaling
;                       2 - User-selectable scaling
;       TRANSFER - The color transfer option:
;                       1 - Linear: color ~ T
;                       2 - Logarithmic: color ~ LOG(T)
;                       3 - Histogram equalized
;       MIN      - The minimum scaling value.
;       MAX      - The maximum scaling value.
;       SILENT   - If present and nonzero, extraneous progress comments
;                  are suppressed.
;       _EXTRA   - Any remaining keywords are passed to "GET_HEAL_LUT";
;                  see that routine for details.
; COMMENTS:
;       The user is prompted for the necessary parameters not supplied
;       by keyword.
;
;       This procedure produces a 2-d image from a 1-d map in 'healpix'
;       pixel list format.  The user must supply the 1-d array and is
;       prompted for projection options from within this routine and 
;       from within GET_HEAL_LUT.  The keyword  mask governs whether 
;       or not unobserved pixels are painted with color 2  (typically 
;       grey). A sample calling sequence with mask set would be:
;
;               HEALPIX_TO_IMAGE, T, I, Mask=N_Obs
;
;       where N_Obs is in the same pixel list format as T, with unobserved 
;       pixels set to zero.  If Mask is not set, the map is assumed to be 
;       full sky.
;
;       Use REPROJ_HEALPIX to convert a HealPix image without any scaling
; PROCEDURES USED:
;       AVG(), GET_HEAL_LUT, GET_HEAL_RES(), LOADCT_012
; MODIFICATION HISTORY:
;       Written by ?
;       SILENT keyword added.  Usage description/header fleshed out.
;          MRG, RITSS, 16 November 2000.
;       Changed keyword_set test to n_elements tests for min and max.
;           RSH, SSAI, 16 Nov 2001.
;       Use STDDEV() instead of STDEV()   W. Landsman Feb. 2003
;       Avoid logarithm of negative numbers W. Landsman Mar. 2003
;       Fix incorrect call to STDDEV()  W. Landsman July 2004
;-
on_error, 2
;
; Check arguments.
If (n_params() LT 2) Then message, 'Syntax: HEALPIX_TO_IMAGE, T, Image'
;
; Check whether graphics is initialized under X windows
IF( !D.Name EQ 'X' AND !D.Window EQ -1 )THEN BEGIN
  ; Loading the color table also sets !D.N_Colors properly
  If (NOT keyword_set(slnt)) Then $
      PRINT,' Loading the default rainbow color table...'
  LOADCT_012, 13
ENDIF
N_Colors = MIN([252,!D.TABLE_SIZE-4])

Scale_Again:
IF( NOT KEYWORD_SET(Scale_Choice) )THEN BEGIN
  PRINT,' Select a display scaling option:'
  PRINT,'    1 -- Auto-scaling'
  PRINT,'    2 -- User-selectable scaling'
  READ, Scale_Choice
ENDIF
IF( Scale_Choice LT 1 OR Scale_Choice GT 2 )THEN BEGIN
  PRINT,' Invalid value for Scale_Choice, please select again...'
  Scale_Choice=0
  GOTO, Scale_Again
ENDIF

Transfer_Again:
IF( NOT KEYWORD_SET(Transfer_Choice) )THEN BEGIN
  PRINT,' Select a color transfer option:'
  PRINT,'    1 -- Linear: color ~ T'
  PRINT,'    2 -- Logarithmic: color ~ LOG(T)'
  PRINT,'    3 -- Histogram equalized'
  READ, Transfer_Choice
ENDIF
IF( Transfer_Choice LT 1 OR Transfer_Choice GT 3 )THEN BEGIN
  PRINT,' Invalid value for Transfer_Choice, please select again...'
  Transfer_Choice=0
  GOTO, Transfer_Again
ENDIF

; Obtain the look-up table
GET_HEAL_LUT, Table, SILENT=slnt, _EXTRA=Extra

; Obtain the pixel reduction factor from the size of the input array
N_Pixels = N_ELEMENTS(T)
Resolution = GET_HEAL_RES( N_Pixels )
IF( Resolution GT 10 OR Resolution LT 1 )THEN BEGIN
  PRINT,' Map resolution too high or low for the current look-up tables...'
  PRINT,' N_Pixels, Resolution:', N_Pixels, Resolution
  STOP
ENDIF ELSE BEGIN
  Factor = 4L^(10-Resolution)
ENDELSE

; Set a mask for determining the display range
Sky = LINDGEN(N_Pixels)
IF( KEYWORD_SET(N_Obs) )THEN BEGIN
  Sky = WHERE(N_Obs GT 0)
ENDIF

; Scale the maps
dflg = (NOT keyword_set(slnt)) OR       $
       (n_elements(Tmin) LE 0) OR (n_elements(Tmax) LE 0)

CASE Scale_Choice OF
 1 : BEGIN
      Tmin = MIN( T[Sky], MAX = Tmax)
      END
 2 : BEGIN
      If (dflg) Then Begin
	  x = moment( T[Sky], Maxmoment=2, sdev = rms)
          PRINT,' The minimum, median, and maximum values in the map are:'
          T_min = MIN( T[Sky], MAX = T_max)
          PRINT, T_min, MEDIAN(T[Sky]), T_max
          PRINT,' The mean and rms of the map are:'
          PRINT, x[0], Rms
      EndIf
      IF( n_elements(Tmin) LE 0 )THEN BEGIN
        PRINT,' Enter the minimum value to be displayed:'
        READ, Tmin
      ENDIF
      IF( n_elements(Tmax) LE 0 )THEN BEGIN
        PRINT,' Enter the maximum value to be displayed:'
        READ, Tmax
      ENDIF
     END
ENDCASE

CASE Transfer_Choice OF
 1 : Color = 3B + BYTSCL(T, Min=Tmin, Max=Tmax, Top=N_Colors )
 2 : BEGIN
       Log_T = alog10(T + (0.0001 - Tmin) )
      Log_Tmin   = min(log_t, max = log_Tmax)
      Color = 3B + BYTSCL(T, Min=Log_Tmin, Max=Log_Tmax, Top=N_Colors )
     END
 3 : Color = 3B + HIST_EQUAL( T, Min=Tmin, Max=Tmax, Top=N_Colors )
ENDCASE

; Paint the unobserved regions with color 2, if desired
IF( KEYWORD_SET(N_Obs) )THEN BEGIN
  If (NOT keyword_set(slnt)) Then $
      MESSAGE,' Masking the blank pixels...',/INF
  IF( MIN(N_Obs) EQ 0 )THEN Color[WHERE(N_Obs EQ 0)] = 2
ENDIF

; Load the map image
If (NOT keyword_set(slnt)) Then $
    MESSAGE,' Generating the image...',/INF
 Image = Color[ (Table > (-1))/Factor ]
 Image[WHERE( Table LT 0 )] = 0

RETURN
END
