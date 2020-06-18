; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_wilcoxon.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION IMSL_wilcoxon, x1, $            ;INPUT 1-D array: floating point
                   x2, $            ;INPUT 1-D array: floating point
                   Fuzz=fuzz, $     ;INPUT Scalar floating point
                   Double=double, $ ;INPUT Scalar ON/OFF flag
                   Stats=stats      ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are a simple arrays.
   ;
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN $
         message, "Incorrect number of arguments."
   size_x1 = IMSL_SIZE(x1)
   IF (size_x1(0) NE 1) THEN BEGIN
      message, "X1 must be a 1-D array."
   END
   IF (nargs EQ 2) THEN BEGIN
      size_x2 = IMSL_SIZE(x2)
      IF (size_x2(0) NE 1) THEN BEGIN
         message, "X2 must be a 1-D array."
      END
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x1(N_ELEMENTS(size_x1)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (nargs EQ 2) THEN BEGIN
      IF (size_x2(N_ELEMENTS(size_x2)-2) GT TYP_FLOAT) THEN type   =   TYP_DOUBLE
   END
   IF (KEYWORD_SET(DOUBLE) EQ true) THEN type   =   TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   nobs1 = size_x1(N_ELEMENTS(size_x1)-1)
   IF (nargs EQ 2) THEN nobs2 = size_x2(N_ELEMENTS(size_x2)-1)
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      IF (nargs EQ 2) THEN result = double(0.0) ELSE result = DBLARR(2)
      ;
      ; Floating point arguments and keywords
      x1_cvt = double(x1)
      IF (nargs EQ 2) THEN x2_cvt = double(x2)
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN $
        fuzz_cvt  =  (DOUBLE(fuzz))(0)
      ; Set the default value for fuzz in case of sign rank test.
      IF ((NOT KEYWORD_SET(fuzz)) AND (nargs EQ 1)) THEN fuzz_cvt = DOUBLE(0.0)
      ;
      ; Output keywords.
      IF (ARG_PRESENT(stats) EQ TRUE) THEN $
        stats_spc = dblarr(10)
   END ELSE BEGIN
      IF (nargs EQ 2) THEN result = float(0.0) ELSE result = FLTARR(2)
      ;
      ; Floating point arguments and keywords
      x1_cvt = float(x1)
      IF (nargs EQ 2) THEN x2_cvt = float(x2)
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN $
        fuzz_cvt = (float(fuzz))(0)
      ; Set the default value for fuzz in case of sign rank test.
      IF ((NOT KEYWORD_SET(fuzz)) AND (nargs EQ 1)) THEN fuzz_cvt = FLOAT(0.0)
      ;
      ; Output keywords.
      IF (ARG_PRESENT(stats) EQ TRUE) THEN $
        stats_spc = fltarr(10)
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_195, type, err_status, x1_cvt, x2_cvt, nobs1, $
                  nobs2, fuzz_cvt, stats_spc, result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(stats) EQ TRUE) THEN $
     stats=stats_spc
   ;
   ; Return.
   ;
   RETURN, result
END
