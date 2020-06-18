; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_kolmogorov1.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_kolmogorov1, f, $            ;INPUT Scalar STRING
                x, $                       ;INPUT 1-D array: floating point
                double = double, $         ;INPUT Scalar ON/OFF flag
                differences=differences, $ ;OUTPUT 1-D array: floating point
                nmissing=nmissing          ;OUTPUT Scalar LONG
    
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - X must be a 1D array
   ;                       
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   nobs = IMSL_N_ELEMENTS(x)
    ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_LONG(0)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = DBLARR(3)
      x_cvt = DOUBLE(x)
      IF (ARG_PRESENT(differences) EQ TRUE) THEN diffs_spc = DBLARR(3)
    END ELSE BEGIN
      result = FLTARR(3)
      x_cvt = FLOAT(x)
      IF (ARG_PRESENT(differences) EQ TRUE) THEN diffs_spc = FLTARR(3)
  END
   ;
   ; Call the system function.
   ;
  err_status = 0L
   MATHSTAT_269, type, err_status, $
                              f, $
                              x_cvt, $
                              nobs, $
                              diffs_spc, $
                              nmissing_spc, $
                              result
   
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc
   IF (ARG_PRESENT(differences) EQ TRUE) THEN differences = diffs_spc
 RETURN, result
END

