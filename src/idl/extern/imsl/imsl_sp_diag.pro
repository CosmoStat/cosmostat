; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_diag.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Sp_diag, a, $
                  n_rows, $
                  symmetric = symmetric, $
                  band = band

@imsl_init.pro
   ON_ERROR, on_err_action

   IF KEYWORD_SET(symmetric) THEN a_sym = IMSL_1 ELSE a_sym = IMSL_0
   IF KEYWORD_SET(band) THEN a_band = IMSL_1 ELSE a_band = IMSL_0
   n = IMSL_LONG(n_rows(0))

   ; Band symmetric
   IF ((a_sym EQ 1) AND (a_band EQ 1)) THEN BEGIN
      RETURN, a(IMSL_2*n:*)
   END

   ;Coordinate, symmetric
   IF ((a_sym EQ 1) AND (a_band EQ 0)) THEN BEGIN
      RETURN, a(where(a.row EQ a.col)).val
   END

   ;Band general
   IF ((a_sym EQ 0) AND (a_band EQ 1)) THEN BEGIN
      RETURN, a(n:IMSL_2*n-1)
   END

   ;Coordinate general
   IF ((a_sym EQ 0) AND (a_band EQ 0)) THEN BEGIN
      RETURN, a(where(a.row EQ a.col)).val
   END
end
